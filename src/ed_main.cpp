#define main ftlm_lanczos_driver_main
#include "main.cpp"
#undef main

namespace {

extern "C" void zheev_(
    char* jobz,
    char* uplo,
    int* n,
    Complex* a,
    int* lda,
    double* w,
    Complex* work,
    int* lwork,
    double* rwork,
    int* info);

std::vector<double> diagonalize_block_exact(const MomentumBlock& block) {
  int n = static_cast<int>(block.basis_dim);
  if (n == 0) {
    return {};
  }
  if (n == 1) {
    return {block.diagonal[0]};
  }

  std::vector<Complex> dense(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    dense[static_cast<std::size_t>(i) + static_cast<std::size_t>(i) * static_cast<std::size_t>(n)] =
        block.diagonal[static_cast<std::size_t>(i)];
    for (int edge = block.row_ptr[static_cast<std::size_t>(i)];
         edge < block.row_ptr[static_cast<std::size_t>(i + 1)];
         ++edge) {
      const int j = block.col_idx[static_cast<std::size_t>(edge)];
      dense[static_cast<std::size_t>(i) + static_cast<std::size_t>(j) * static_cast<std::size_t>(n)] +=
          block.values[static_cast<std::size_t>(edge)];
    }
  }

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

  std::vector<double> eigenvalues(static_cast<std::size_t>(n), 0.0);
  std::vector<double> rwork(static_cast<std::size_t>(std::max(1, 3 * n - 2)), 0.0);
  int lda = n;
  int lwork = -1;
  int info = 0;
  Complex work_query = 0.0;
  char jobz = 'N';
  char uplo = 'U';
  zheev_(&jobz, &uplo, &n, dense.data(), &lda, eigenvalues.data(),
         &work_query, &lwork, rwork.data(), &info);
  if (info != 0) {
    die("zheev workspace query failed with info=" + std::to_string(info));
  }

  lwork = std::max(1, static_cast<int>(std::real(work_query)));
  std::vector<Complex> work(static_cast<std::size_t>(lwork), 0.0);
  zheev_(&jobz, &uplo, &n, dense.data(), &lda, eigenvalues.data(),
         work.data(), &lwork, rwork.data(), &info);
  if (info != 0) {
    die("zheev diagonalization failed with info=" + std::to_string(info));
  }

  return eigenvalues;
}

void accumulate_exact_block(
    const std::vector<double>& eigenvalues,
    int particles,
    double beta,
    const std::vector<double>& mu_values,
    std::vector<double>& emin_by_mu,
    std::vector<double>& z_scaled,
    std::vector<double>& n_scaled) {
  for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
    const double mu = mu_values[imu];
    double block_min = std::numeric_limits<double>::infinity();
    for (double eigenvalue : eigenvalues) {
      block_min = std::min(block_min, eigenvalue - mu * static_cast<double>(particles));
    }
    if (!std::isfinite(block_min)) {
      continue;
    }

    if (!std::isfinite(emin_by_mu[imu])) {
      emin_by_mu[imu] = block_min;
    } else if (block_min < emin_by_mu[imu]) {
      const double scale = std::exp(-beta * (emin_by_mu[imu] - block_min));
      z_scaled[imu] *= scale;
      n_scaled[imu] *= scale;
      emin_by_mu[imu] = block_min;
    }

    for (double eigenvalue : eigenvalues) {
      const double shifted = eigenvalue - mu * static_cast<double>(particles) - emin_by_mu[imu];
      const double weight = std::exp(-beta * shifted);
      z_scaled[imu] += weight;
      n_scaled[imu] += static_cast<double>(particles) * weight;
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Params params = parse_args(argc, argv);
    const Lattice lattice = build_lattice(params.lx, params.ly);
    const std::vector<double> mu_values =
        linspace(params.mu_min, params.mu_max, params.mu_count);

    std::cout << "Exact diagonalization using the FTLM sector and k-block builder\n";
    std::cout << "Running ED on a " << params.lx << "x" << params.ly
              << " rectangular lattice (" << lattice.sites << " sites)\n";
    std::cout << "beta=" << params.beta << " U=" << params.u << " V=" << params.v
              << " tx=" << params.tx << " ty=" << params.ty << "\n";

    std::vector<double> densities(mu_values.size(), 0.0);
    std::vector<double> partitions(mu_values.size(), 0.0);
    std::vector<double> emin_by_mu(
        mu_values.size(), std::numeric_limits<double>::infinity());
    std::vector<double> z_scaled(mu_values.size(), 0.0);
    std::vector<double> n_scaled(mu_values.size(), 0.0);

    for (int n_up = 0; n_up <= lattice.sites; ++n_up) {
      for (int n_down = 0; n_down <= lattice.sites; ++n_down) {
        const SectorBasis sector = build_sector_basis(lattice, n_up, n_down, params);
        if (sector.full_dim == 0) {
          continue;
        }
        if (sector.full_dim > params.max_sector_dim) {
          std::ostringstream msg;
          msg << "Sector (Nup=" << n_up << ", Ndown=" << n_down << ") has full dimension "
              << sector.full_dim << ", which exceeds --max-sector-dim=" << params.max_sector_dim;
          die(msg.str());
        }

        std::size_t active_sum = 0;
        for (int my = 0; my < lattice.ly; ++my) {
          for (int mx = 0; mx < lattice.lx; ++mx) {
            const MomentumBlock block = build_momentum_block(sector, lattice, mx, my);
            active_sum += block.basis_dim;
            if (block.basis_dim == 0) {
              continue;
            }
            const std::vector<double> eigenvalues = diagonalize_block_exact(block);
            accumulate_exact_block(
                eigenvalues,
                block.particles,
                params.beta,
                mu_values,
                emin_by_mu,
                z_scaled,
                n_scaled);
          }
        }

        std::cout << "  sector Nup=" << n_up << " Ndown=" << n_down
                  << " full_dim=" << sector.full_dim
                  << " parents=" << sector.parents.size()
                  << " active_k_dim_sum=" << active_sum << "\n";
      }
    }

    for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
      partitions[imu] = z_scaled[imu];
      densities[imu] =
          n_scaled[imu] / (static_cast<double>(lattice.sites) * z_scaled[imu]);
      std::cout << "mu=" << std::setw(12) << mu_values[imu] << "  n=" << std::setw(12)
                << densities[imu] << "\n";
    }

    write_results(params.output, mu_values, densities, partitions);
    std::cout << "Wrote " << params.output << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }
}
