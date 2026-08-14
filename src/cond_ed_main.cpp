#define main ftlm_lanczos_driver_main
#include "main.cpp"
#undef main

namespace {

struct ConductivityParams {
  Params base;
  double omega_min = 0.0;
  double omega_max = 10.0;
  int omega_count = 501;
  double eta = 0.05;
  char direction = 'x';
  std::size_t max_block_dim = 1024;
  std::string dc_output;
};

struct EigenBlock {
  int particles = 0;
  std::size_t basis_dim = 0;
  std::vector<double> eigenvalues;
  std::vector<double> current_abs2;  // row + col * n in H eigenbasis
};

struct ThermoPoint {
  double density = 0.0;
  double charge_correlation = 0.0;
  double compressibility = 0.0;
  double log_partition = 0.0;
};

void print_cond_help(const char* argv0) {
  print_help(argv0);
  std::cout
      << "\nConductivity ED options:\n"
      << "  --omega-min X                  minimum frequency\n"
      << "  --omega-max X                  maximum frequency\n"
      << "  --omega-count N                number of frequency points\n"
      << "  --eta X                        Lorentzian broadening\n"
      << "  --conductivity-direction x|y   current direction\n"
      << "  --max-conductivity-block-dim N abort above this exact block size\n"
      << "  --dc-output PATH               optional omega=0 summary CSV\n";
}

ConductivityParams parse_cond_args(int argc, char** argv) {
  ConductivityParams params;
  std::vector<std::string> base_args;
  base_args.push_back(argv[0]);
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto next = [&](const std::string& name) {
      if (i + 1 >= argc) {
        die("Missing value for " + name);
      }
      return std::string(argv[++i]);
    };
    if (arg == "--help") {
      print_cond_help(argv[0]);
      std::exit(0);
    } else if (arg == "--omega-min") {
      params.omega_min = parse_value<double>(next(arg));
    } else if (arg == "--omega-max") {
      params.omega_max = parse_value<double>(next(arg));
    } else if (arg == "--omega-count") {
      params.omega_count = parse_value<int>(next(arg));
    } else if (arg == "--eta") {
      params.eta = parse_value<double>(next(arg));
    } else if (arg == "--conductivity-direction") {
      const std::string value = next(arg);
      if (value != "x" && value != "y") {
        die("--conductivity-direction must be x or y.");
      }
      params.direction = value[0];
    } else if (arg == "--max-conductivity-block-dim") {
      params.max_block_dim = parse_value<std::size_t>(next(arg));
    } else if (arg == "--dc-output") {
      params.dc_output = next(arg);
    } else {
      base_args.push_back(arg);
      if (arg != "--help") {
        // Options known to parse_args that require values.
        static const std::vector<std::string> value_args = {
            "--lx", "--ly", "--tx", "--ty", "--tp", "--phix", "--phiy",
            "--u", "--v", "--beta", "--beta-list", "--mu-min", "--mu-max",
            "--mu-count", "--samples", "--lanczos-steps", "--threads",
            "--exact-block-threshold", "--seed", "--max-sector-dim",
            "--checkpoint", "--trace-partition-csv", "--debug-block-nup",
            "--debug-block-ndown", "--debug-block-mx", "--debug-block-my",
            "--debug-block-csv", "--output"};
        if (std::find(value_args.begin(), value_args.end(), arg) != value_args.end()) {
          base_args.push_back(next(arg));
        }
      }
    }
  }

  std::vector<char*> argv_copy;
  argv_copy.reserve(base_args.size());
  for (std::string& arg : base_args) {
    argv_copy.push_back(arg.data());
  }
  params.base = parse_args(static_cast<int>(argv_copy.size()), argv_copy.data());

  if (params.omega_count < 2) {
    die("--omega-count must be at least 2.");
  }
  if (params.omega_max <= params.omega_min) {
    die("--omega-max must be greater than --omega-min.");
  }
  if (params.eta <= 0.0) {
    die("--eta must be positive.");
  }
  return params;
}

double direction_component(int dir, char direction) {
  const int dx[8] = {1, -1, 0, 0, 1, -1, 1, -1};
  const int dy[8] = {0, 0, 1, -1, 1, 1, -1, -1};
  return static_cast<double>((direction == 'x') ? dx[dir] : dy[dir]);
}

MomentumBlock build_current_block(
    const SectorBasis& sector,
    const Lattice& lattice,
    int mx,
    int my,
    char direction) {
  std::vector<std::pair<int, int>> translations;
  translations.reserve(static_cast<std::size_t>(lattice.sites));
  for (int dy = 0; dy < lattice.ly; ++dy) {
    for (int dx = 0; dx < lattice.lx; ++dx) {
      translations.push_back({dx, dy});
    }
  }

  std::vector<Complex> phase(translations.size(), Complex{0.0, 0.0});
  for (std::size_t i = 0; i < translations.size(); ++i) {
    const double arg =
        2.0 * kPi *
        (static_cast<double>(mx * translations[i].first) / static_cast<double>(lattice.lx) +
         static_cast<double>(my * translations[i].second) / static_cast<double>(lattice.ly));
    phase[i] = Complex(std::cos(arg), std::sin(arg));
  }

  std::vector<double> dnf(sector.parents.size(), 0.0);
  std::vector<int> active_index(sector.parents.size(), -1);
  int next_active = 0;
  for (std::size_t p = 0; p < sector.parents.size(); ++p) {
    Complex phsum = 0.0;
    const StateKey& rep = sector.parents[p].representative;
    for (std::size_t shift = 0; shift < translations.size(); ++shift) {
      int sign = 1;
      const StateKey translated = translate_state(
          rep,
          translations[shift].first,
          translations[shift].second,
          lattice.lx,
          lattice.ly,
          sign);
      if (translated == rep) {
        phsum += static_cast<double>(sign) * phase[shift];
      }
    }
    dnf[p] = std::abs(phsum) / std::sqrt(static_cast<double>(sector.parents[p].degeneracy));
    if (dnf[p] > kEps) {
      active_index[p] = next_active++;
    }
  }

  MomentumBlock block;
  block.mx = mx;
  block.my = my;
  block.particles = sector.particles;
  block.basis_dim = static_cast<std::size_t>(next_active);
  if (block.basis_dim == 0) {
    return block;
  }
  block.diagonal.assign(block.basis_dim, 0.0);
  block.row_ptr.reserve(block.basis_dim + 1);
  block.row_ptr.push_back(0);

  for (std::size_t p = 0; p < sector.parents.size(); ++p) {
    const int src = active_index[p];
    if (src < 0) {
      continue;
    }
    for (int edge = sector.hop_row_ptr[p]; edge < sector.hop_row_ptr[p + 1]; ++edge) {
      const CompactHop& hop = sector.hops[static_cast<std::size_t>(edge)];
      const double component = direction_component(hop.dir, direction);
      if (std::abs(component) < kEps) {
        continue;
      }
      const int tgt = active_index[static_cast<std::size_t>(hop.parent)];
      if (tgt < 0) {
        continue;
      }
      const Complex h_coeff =
          sector.hop_amplitudes[static_cast<std::size_t>(hop.dir)] *
          static_cast<double>(hop.sign) *
          phase[static_cast<std::size_t>(hop.shift)] *
          (dnf[static_cast<std::size_t>(hop.parent)] / dnf[p]);
      const Complex j_coeff = Complex{0.0, component} * h_coeff;
      block.col_idx.push_back(tgt);
      block.values.push_back(j_coeff);
    }
    block.row_ptr.push_back(static_cast<int>(block.col_idx.size()));
  }
  return block;
}

std::vector<Complex> dense_from_block(
    const MomentumBlock& block,
    bool include_diagonal,
    bool hermitize) {
  const int n = static_cast<int>(block.basis_dim);
  std::vector<Complex> dense(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    if (include_diagonal) {
      dense[static_cast<std::size_t>(i) + static_cast<std::size_t>(i) * n] =
          block.diagonal[static_cast<std::size_t>(i)];
    }
    for (int edge = block.row_ptr[static_cast<std::size_t>(i)];
         edge < block.row_ptr[static_cast<std::size_t>(i + 1)];
         ++edge) {
      const int j = block.col_idx[static_cast<std::size_t>(edge)];
      dense[static_cast<std::size_t>(i) + static_cast<std::size_t>(j) * n] +=
          block.values[static_cast<std::size_t>(edge)];
    }
  }
  if (hermitize) {
    for (int j = 0; j < n; ++j) {
      const std::size_t jj = static_cast<std::size_t>(j);
      dense[jj + jj * static_cast<std::size_t>(n)] =
          Complex(dense[jj + jj * static_cast<std::size_t>(n)].real(), 0.0);
      for (int i = 0; i < j; ++i) {
        const std::size_t idx_ij = static_cast<std::size_t>(i) + jj * static_cast<std::size_t>(n);
        const std::size_t idx_ji = jj + static_cast<std::size_t>(i) * static_cast<std::size_t>(n);
        const Complex avg = 0.5 * (dense[idx_ij] + std::conj(dense[idx_ji]));
        dense[idx_ij] = avg;
        dense[idx_ji] = std::conj(avg);
      }
    }
  }
  return dense;
}

void diagonalize_dense_hermitian(
    std::vector<Complex>& dense,
    std::vector<double>& eigenvalues) {
  int n = static_cast<int>(eigenvalues.size());
  if (n == 0) {
    return;
  }
  std::vector<double> rwork(static_cast<std::size_t>(std::max(1, 3 * n - 2)), 0.0);
  int lda = n;
  int lwork = -1;
  int info = 0;
  Complex work_query = 0.0;
  char jobz = 'V';
  char uplo = 'U';
  zheev_(&jobz, &uplo, &n, dense.data(), &lda, eigenvalues.data(),
         &work_query, &lwork, rwork.data(), &info);
  if (info != 0) {
    die("zheev workspace query failed with info=" + std::to_string(info));
  }
  lwork = std::max(1, static_cast<int>(std::ceil(work_query.real())));
  std::vector<Complex> work(static_cast<std::size_t>(lwork), 0.0);
  zheev_(&jobz, &uplo, &n, dense.data(), &lda, eigenvalues.data(),
         work.data(), &lwork, rwork.data(), &info);
  if (info != 0) {
    die("zheev diagonalization failed with info=" + std::to_string(info));
  }
}

std::vector<double> transform_current_abs2(
    const std::vector<Complex>& eigenvectors,
    const std::vector<Complex>& current_dense,
    int n) {
  std::vector<Complex> temp(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      Complex sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += current_dense[static_cast<std::size_t>(row) + static_cast<std::size_t>(k) * n] *
               eigenvectors[static_cast<std::size_t>(k) + static_cast<std::size_t>(col) * n];
      }
      temp[static_cast<std::size_t>(row) + static_cast<std::size_t>(col) * n] = sum;
    }
  }

  std::vector<double> abs2(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < n; ++row) {
      Complex sum = 0.0;
      for (int k = 0; k < n; ++k) {
        sum += std::conj(eigenvectors[static_cast<std::size_t>(k) + static_cast<std::size_t>(row) * n]) *
               temp[static_cast<std::size_t>(k) + static_cast<std::size_t>(col) * n];
      }
      abs2[static_cast<std::size_t>(row) + static_cast<std::size_t>(col) * n] = std::norm(sum);
    }
  }
  return abs2;
}

double logsumexp(const std::vector<double>& values) {
  const double top = *std::max_element(values.begin(), values.end());
  double sum = 0.0;
  for (double value : values) {
    sum += std::exp(value - top);
  }
  return top + std::log(sum);
}

ThermoPoint thermo_for(
    const std::vector<EigenBlock>& blocks,
    double beta,
    double mu,
    int sites) {
  std::vector<double> log_weights;
  log_weights.reserve(1024);
  for (const EigenBlock& block : blocks) {
    for (double eigenvalue : block.eigenvalues) {
      log_weights.push_back(-beta * (eigenvalue - mu * static_cast<double>(block.particles)));
    }
  }
  const double log_z = logsumexp(log_weights);
  double n_total = 0.0;
  double n2_total = 0.0;
  for (const EigenBlock& block : blocks) {
    const double particles = static_cast<double>(block.particles);
    for (double eigenvalue : block.eigenvalues) {
      const double weight = std::exp(-beta * (eigenvalue - mu * particles) - log_z);
      n_total += weight * particles;
      n2_total += weight * particles * particles;
    }
  }
  ThermoPoint point;
  point.log_partition = log_z;
  point.density = n_total / static_cast<double>(sites);
  point.charge_correlation =
      std::max(0.0, (n2_total - n_total * n_total) / static_cast<double>(sites));
  point.compressibility = beta * point.charge_correlation;
  return point;
}

double sigma_for_omega(
    const std::vector<EigenBlock>& blocks,
    double beta,
    double mu,
    double omega,
    double eta,
    double log_z,
    int sites) {
  double sigma = 0.0;
  constexpr double dE_eps = 1e-10;
  for (const EigenBlock& block : blocks) {
    const int n = static_cast<int>(block.basis_dim);
    const double particles = static_cast<double>(block.particles);
    for (int initial = 0; initial < n; ++initial) {
      const double ei = block.eigenvalues[static_cast<std::size_t>(initial)];
      const double log_wi = -beta * (ei - mu * particles) - log_z;
      const double wi = std::exp(log_wi);
      for (int final = 0; final < n; ++final) {
        const double ej = block.eigenvalues[static_cast<std::size_t>(final)];
        const double delta_e = ej - ei;
        if (delta_e < -dE_eps) {
          continue;
        }
        const double j2 =
            block.current_abs2[static_cast<std::size_t>(final) +
                               static_cast<std::size_t>(initial) * n];
        if (j2 <= 0.0) {
          continue;
        }
        double thermal_factor = 0.0;
        if (std::abs(delta_e) <= dE_eps) {
          thermal_factor = beta * wi;
        } else {
          const double wj =
              std::exp(-beta * (ej - mu * particles) - log_z);
          thermal_factor = (wi - wj) / delta_e;
        }
        const double broadened_pi_delta =
            eta / ((omega - delta_e) * (omega - delta_e) + eta * eta);
        sigma += thermal_factor * j2 * broadened_pi_delta;
      }
    }
  }
  return sigma / static_cast<double>(sites);
}

void write_conductivity_csv(
    const std::string& path,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    const std::vector<double>& omega_values,
    const std::vector<std::vector<ThermoPoint>>& thermo,
    const std::vector<std::vector<std::vector<double>>>& sigma) {
  std::ofstream out(path);
  if (!out) {
    die("Failed to open conductivity output file: " + path);
  }
  out << "beta,mu,omega,sigma,n,charge_correlation,compressibility,log_partition\n";
  out << std::setprecision(15);
  for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
    for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
      for (std::size_t iw = 0; iw < omega_values.size(); ++iw) {
        const ThermoPoint& point = thermo[ibeta][imu];
        out << beta_values[ibeta] << "," << mu_values[imu] << ","
            << omega_values[iw] << "," << sigma[ibeta][imu][iw] << ","
            << point.density << "," << point.charge_correlation << ","
            << point.compressibility << "," << point.log_partition << "\n";
      }
    }
  }
}

void write_dc_csv(
    const std::string& path,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    const std::vector<std::vector<ThermoPoint>>& thermo,
    const std::vector<std::vector<double>>& dc_sigma) {
  if (path.empty()) {
    return;
  }
  std::ofstream out(path);
  if (!out) {
    die("Failed to open DC output file: " + path);
  }
  out << "beta,mu,sigma_dc,n,charge_correlation,compressibility,log_partition\n";
  out << std::setprecision(15);
  for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
    for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
      const ThermoPoint& point = thermo[ibeta][imu];
      out << beta_values[ibeta] << "," << mu_values[imu] << ","
          << dc_sigma[ibeta][imu] << "," << point.density << ","
          << point.charge_correlation << "," << point.compressibility << ","
          << point.log_partition << "\n";
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const ConductivityParams params = parse_cond_args(argc, argv);
    const Params& base = params.base;
    const Lattice lattice = build_lattice(base.lx, base.ly);
    const std::vector<double> mu_values =
        linspace(base.mu_min, base.mu_max, base.mu_count);
    const std::vector<double> omega_values =
        linspace(params.omega_min, params.omega_max, params.omega_count);

    std::cout << "Running exact optical/DC conductivity on a "
              << base.lx << "x" << base.ly << " rectangular lattice ("
              << lattice.sites << " sites)\n";
    std::cout << "direction=" << params.direction
              << " eta=" << params.eta
              << " omega=[" << params.omega_min << "," << params.omega_max
              << "] count=" << params.omega_count
              << " max_block_dim=" << params.max_block_dim << "\n";

    std::vector<EigenBlock> blocks;
    for (int n_up = 0; n_up <= lattice.sites; ++n_up) {
      for (int n_down = 0; n_down <= lattice.sites; ++n_down) {
        const SectorBasis sector = build_sector_basis(lattice, n_up, n_down, base);
        if (sector.full_dim == 0) {
          continue;
        }
        if (sector.full_dim > base.max_sector_dim) {
          std::ostringstream msg;
          msg << "Sector (Nup=" << n_up << ", Ndown=" << n_down
              << ") has full dimension " << sector.full_dim
              << ", which exceeds --max-sector-dim=" << base.max_sector_dim;
          die(msg.str());
        }
        std::size_t active_sum = 0;
        for (int my = 0; my < lattice.ly; ++my) {
          for (int mx = 0; mx < lattice.lx; ++mx) {
            const MomentumBlock h_block = build_momentum_block(sector, lattice, mx, my);
            active_sum += h_block.basis_dim;
            if (h_block.basis_dim == 0) {
              continue;
            }
            if (h_block.basis_dim > params.max_block_dim) {
              std::ostringstream msg;
              msg << "Conductivity block (Nup=" << n_up << ", Ndown=" << n_down
                  << ", mx=" << mx << ", my=" << my << ") has dimension "
                  << h_block.basis_dim
                  << ", which exceeds --max-conductivity-block-dim="
                  << params.max_block_dim;
              die(msg.str());
            }

            const MomentumBlock j_block =
                build_current_block(sector, lattice, mx, my, params.direction);
            std::vector<Complex> h_dense = dense_from_block(h_block, true, true);
            std::vector<Complex> j_dense = dense_from_block(j_block, false, false);
            std::vector<double> eigenvalues(h_block.basis_dim, 0.0);
            diagonalize_dense_hermitian(h_dense, eigenvalues);

            EigenBlock block;
            block.particles = h_block.particles;
            block.basis_dim = h_block.basis_dim;
            block.eigenvalues = std::move(eigenvalues);
            block.current_abs2 = transform_current_abs2(
                h_dense, j_dense, static_cast<int>(h_block.basis_dim));
            blocks.push_back(std::move(block));
          }
        }
        std::cout << "  sector Nup=" << n_up << " Ndown=" << n_down
                  << " full_dim=" << sector.full_dim
                  << " parents=" << sector.parents.size()
                  << " active_k_dim_sum=" << active_sum << "\n";
      }
    }

    std::vector<std::vector<ThermoPoint>> thermo(
        base.beta_values.size(), std::vector<ThermoPoint>(mu_values.size()));
    std::vector<std::vector<std::vector<double>>> sigma(
        base.beta_values.size(),
        std::vector<std::vector<double>>(
            mu_values.size(), std::vector<double>(omega_values.size(), 0.0)));
    std::vector<std::vector<double>> dc_sigma(
        base.beta_values.size(), std::vector<double>(mu_values.size(), 0.0));

    for (std::size_t ibeta = 0; ibeta < base.beta_values.size(); ++ibeta) {
      const double beta = base.beta_values[ibeta];
      for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
        const double mu = mu_values[imu];
        thermo[ibeta][imu] = thermo_for(blocks, beta, mu, lattice.sites);
        for (std::size_t iw = 0; iw < omega_values.size(); ++iw) {
          sigma[ibeta][imu][iw] = sigma_for_omega(
              blocks,
              beta,
              mu,
              omega_values[iw],
              params.eta,
              thermo[ibeta][imu].log_partition,
              lattice.sites);
        }
        dc_sigma[ibeta][imu] = sigma_for_omega(
            blocks,
            beta,
            mu,
            0.0,
            params.eta,
            thermo[ibeta][imu].log_partition,
            lattice.sites);
        std::cout << "beta=" << std::setw(12) << beta
                  << "  mu=" << std::setw(12) << mu
                  << "  n=" << std::setw(12) << thermo[ibeta][imu].density
                  << "  sigma_dc=" << std::setw(12) << dc_sigma[ibeta][imu]
                  << "\n";
      }
    }

    write_conductivity_csv(
        base.output,
        base.beta_values,
        mu_values,
        omega_values,
        thermo,
        sigma);
    write_dc_csv(params.dc_output, base.beta_values, mu_values, thermo, dc_sigma);
    std::cout << "Wrote " << base.output << "\n";
    if (!params.dc_output.empty()) {
      std::cout << "Wrote " << params.dc_output << "\n";
    }
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }
}
