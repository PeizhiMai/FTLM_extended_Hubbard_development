#include <chrono>
#include <functional>
#include <unordered_set>

#define main ftlm_lanczos_driver_main
#include "main.cpp"
#undef main

namespace {

enum class ConductivityMemoryMode {
  Hybrid,
  StoreBases,
  Recompute,
};

struct ConductivityParams {
  Params base;
  double omega_min = 0.0;
  double omega_max = 10.0;
  int omega_count = 501;
  double eta = 0.05;
  char direction = 'x';
  std::string dc_output;
  std::string conductivity_checkpoint;
  ConductivityMemoryMode memory_mode = ConductivityMemoryMode::Hybrid;
  double basis_memory_mb = 512.0;
  double max_runtime_minutes = 0.0;
};

struct ExactConductivityBlock {
  int particles = 0;
  std::size_t basis_dim = 0;
  std::vector<double> eigenvalues;
  std::vector<double> current_abs2;  // final + initial * n in H eigenbasis
};

struct ThermoPoint {
  double density = 0.0;
  double charge_correlation = 0.0;
  double compressibility = 0.0;
  double log_partition = -std::numeric_limits<double>::infinity();
};

struct ConductivitySampleKey {
  BlockKey block;
  int sample = 0;

  bool operator==(const ConductivitySampleKey& other) const {
    return block == other.block && sample == other.sample;
  }
};

struct ConductivitySampleKeyHash {
  std::size_t operator()(const ConductivitySampleKey& key) const {
    std::size_t value = BlockKeyHash{}(key.block);
    value ^= std::hash<int>{}(key.sample) + 0x9e3779b9U + (value << 6U) + (value >> 2U);
    return value;
  }
};

struct ConductivitySampleRecord {
  ConductivitySampleKey key;
  int particles = 0;
  std::size_t basis_dim = 0;
  double trace_prefactor = 1.0;
  double current_norm = 0.0;
  std::vector<double> first_eigenvalues;
  std::vector<double> first_v0;
  std::vector<double> second_eigenvalues;
  std::vector<double> second_v0;
  std::vector<Complex> current_matrix;  // second_eigen + first_eigen * second_count
};

struct LanczosDecomposition {
  int steps = 0;
  std::vector<double> alpha;
  std::vector<double> beta;
  std::vector<double> eigenvalues;
  std::vector<double> first_components;
  std::vector<std::vector<double>> eigenvectors;  // Lanczos index + eigen index * steps
  std::vector<std::vector<Complex>> basis;
};

struct LogAccumulator {
  double log_sum = -std::numeric_limits<double>::infinity();

  void add_log(double value) {
    if (!std::isfinite(value)) {
      return;
    }
    if (!std::isfinite(log_sum)) {
      log_sum = value;
      return;
    }
    if (value > log_sum) {
      log_sum = value + std::log1p(std::exp(log_sum - value));
    } else {
      log_sum = log_sum + std::log1p(std::exp(value - log_sum));
    }
  }
};

const char* memory_mode_name(ConductivityMemoryMode mode) {
  switch (mode) {
    case ConductivityMemoryMode::Hybrid:
      return "hybrid";
    case ConductivityMemoryMode::StoreBases:
      return "store-bases";
    case ConductivityMemoryMode::Recompute:
      return "recompute";
  }
  return "unknown";
}

ConductivityMemoryMode parse_memory_mode(const std::string& value) {
  if (value == "hybrid") {
    return ConductivityMemoryMode::Hybrid;
  }
  if (value == "store-bases") {
    return ConductivityMemoryMode::StoreBases;
  }
  if (value == "recompute") {
    return ConductivityMemoryMode::Recompute;
  }
  die("--conductivity-memory-mode must be hybrid, store-bases, or recompute.");
}

void print_cond_help(const char* argv0) {
  print_help(argv0);
  std::cout
      << "\nFTLM conductivity options:\n"
      << "  --omega-min X                       minimum frequency\n"
      << "  --omega-max X                       maximum frequency\n"
      << "  --omega-count N                     number of frequency points\n"
      << "  --eta X                             Lorentzian broadening\n"
      << "  --conductivity-direction x|y        current direction\n"
      << "  --conductivity-checkpoint PATH      append/resume per-sample conductivity records\n"
      << "  --conductivity-memory-mode MODE     hybrid, store-bases, or recompute\n"
      << "  --conductivity-basis-memory-mb X    hybrid memory budget for storing two bases\n"
      << "  --max-runtime-minutes X             stop cleanly after completed sample records\n"
      << "  --dc-output PATH                    optional omega=0 summary CSV\n";
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
    } else if (arg == "--conductivity-checkpoint") {
      params.conductivity_checkpoint = next(arg);
    } else if (arg == "--conductivity-memory-mode") {
      params.memory_mode = parse_memory_mode(next(arg));
    } else if (arg == "--conductivity-basis-memory-mb") {
      params.basis_memory_mb = parse_value<double>(next(arg));
    } else if (arg == "--max-runtime-minutes") {
      params.max_runtime_minutes = parse_value<double>(next(arg));
    } else if (arg == "--dc-output") {
      params.dc_output = next(arg);
    } else {
      base_args.push_back(arg);
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

  std::vector<char*> argv_copy;
  argv_copy.reserve(base_args.size());
  for (std::string& arg : base_args) {
    argv_copy.push_back(arg.data());
  }
  params.base = parse_args(static_cast<int>(argv_copy.size()), argv_copy.data());

  if (!params.base.checkpoint.empty()) {
    die("Use --conductivity-checkpoint for conductivity runs; --checkpoint is reserved for thermodynamic FTLM.");
  }
  if (!params.base.trace_partition_csv.empty() || !params.base.debug_block_csv.empty() ||
      params.base.debug_block_nup >= 0 || params.base.debug_block_ndown >= 0 ||
      params.base.debug_block_mx >= 0 || params.base.debug_block_my >= 0) {
    die("Thermodynamic trace/debug options are not supported by ftlm_conductivity.");
  }
  if (params.omega_count < 2) {
    die("--omega-count must be at least 2.");
  }
  if (params.omega_max <= params.omega_min) {
    die("--omega-max must be greater than --omega-min.");
  }
  if (params.eta <= 0.0) {
    die("--eta must be positive.");
  }
  if (params.basis_memory_mb <= 0.0) {
    die("--conductivity-basis-memory-mb must be positive.");
  }
  if (params.max_runtime_minutes < 0.0) {
    die("--max-runtime-minutes must be non-negative.");
  }
  if (params.max_runtime_minutes > 0.0 && params.conductivity_checkpoint.empty()) {
    die("--max-runtime-minutes requires --conductivity-checkpoint so completed samples are kept.");
  }
  return params;
}

double direction_component(int dir, char direction) {
  const int dx[8] = {1, -1, 0, 0, 1, -1, 1, -1};
  const int dy[8] = {0, 0, 1, -1, 1, 1, -1, -1};
  if (dir < 0 || dir >= 8) {
    die("Invalid hopping direction in current operator.");
  }
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

void apply_operator(
    const MomentumBlock& op,
    const std::vector<Complex>& in,
    std::vector<Complex>& out) {
  std::fill(out.begin(), out.end(), Complex{0.0, 0.0});
  for (std::size_t i = 0; i < op.basis_dim; ++i) {
    out[i] += op.diagonal[i] * in[i];
    for (int edge = op.row_ptr[i]; edge < op.row_ptr[i + 1]; ++edge) {
      out[i] += op.values[static_cast<std::size_t>(edge)] *
                in[static_cast<std::size_t>(op.col_idx[static_cast<std::size_t>(edge)])];
    }
  }
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

ExactConductivityBlock diagonalize_exact_conductivity_block(
    const MomentumBlock& h_block,
    const MomentumBlock& j_block) {
  ExactConductivityBlock block;
  block.particles = h_block.particles;
  block.basis_dim = h_block.basis_dim;
  std::vector<Complex> h_dense = dense_from_block(h_block, true, true);
  std::vector<Complex> j_dense = dense_from_block(j_block, false, false);
  block.eigenvalues.assign(h_block.basis_dim, 0.0);
  diagonalize_dense_hermitian(h_dense, block.eigenvalues);
  block.current_abs2 = transform_current_abs2(
      h_dense, j_dense, static_cast<int>(h_block.basis_dim));
  return block;
}

LanczosDecomposition run_lanczos_decomposition(
    const MomentumBlock& block,
    const std::vector<Complex>& initial,
    int krylov_dim,
    bool store_basis) {
  LanczosDecomposition decomp;
  if (block.basis_dim == 0 || krylov_dim <= 0) {
    return decomp;
  }

  std::vector<Complex> q_prev(block.basis_dim, 0.0);
  std::vector<Complex> q_cur = initial;
  std::vector<Complex> q_next(block.basis_dim, 0.0);
  std::vector<Complex> hq(block.basis_dim, 0.0);
  decomp.alpha.reserve(static_cast<std::size_t>(krylov_dim));
  decomp.beta.reserve(static_cast<std::size_t>(krylov_dim));
  if (store_basis) {
    decomp.basis.reserve(static_cast<std::size_t>(krylov_dim));
  }

  double beta_prev = 0.0;
  for (int step = 0; step < krylov_dim; ++step) {
    if (store_basis) {
      decomp.basis.push_back(q_cur);
    }
    apply_hamiltonian(block, q_cur, hq);
    const double alpha = dot(q_cur, hq).real();
    for (std::size_t i = 0; i < block.basis_dim; ++i) {
      q_next[i] = hq[i] - alpha * q_cur[i] - beta_prev * q_prev[i];
    }
    decomp.alpha.push_back(alpha);
    ++decomp.steps;
    const double beta = std::sqrt(complex_norm2(q_next));
    if (beta < 1e-12 || step + 1 == krylov_dim) {
      break;
    }
    decomp.beta.push_back(beta);
    q_prev.swap(q_cur);
    q_cur.swap(q_next);
    for (Complex& value : q_cur) {
      value /= beta;
    }
    beta_prev = beta;
  }

  std::vector<std::vector<double>> tri(
      static_cast<std::size_t>(decomp.steps),
      std::vector<double>(static_cast<std::size_t>(decomp.steps), 0.0));
  for (int i = 0; i < decomp.steps; ++i) {
    tri[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] =
        decomp.alpha[static_cast<std::size_t>(i)];
    if (i + 1 < decomp.steps) {
      tri[static_cast<std::size_t>(i)][static_cast<std::size_t>(i + 1)] =
          decomp.beta[static_cast<std::size_t>(i)];
      tri[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(i)] =
          decomp.beta[static_cast<std::size_t>(i)];
    }
  }

  std::vector<double> raw_eigenvalues;
  std::vector<std::vector<double>> raw_eigenvectors;
  jacobi_diagonalize(tri, raw_eigenvalues, raw_eigenvectors);

  std::vector<int> order(decomp.steps);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
    return raw_eigenvalues[static_cast<std::size_t>(lhs)] <
           raw_eigenvalues[static_cast<std::size_t>(rhs)];
  });

  decomp.eigenvalues.resize(static_cast<std::size_t>(decomp.steps));
  decomp.first_components.resize(static_cast<std::size_t>(decomp.steps));
  decomp.eigenvectors.assign(
      static_cast<std::size_t>(decomp.steps),
      std::vector<double>(static_cast<std::size_t>(decomp.steps), 0.0));
  for (int sorted = 0; sorted < decomp.steps; ++sorted) {
    const int original = order[static_cast<std::size_t>(sorted)];
    decomp.eigenvalues[static_cast<std::size_t>(sorted)] =
        raw_eigenvalues[static_cast<std::size_t>(original)];
    decomp.first_components[static_cast<std::size_t>(sorted)] =
        raw_eigenvectors[0][static_cast<std::size_t>(original)];
    for (int row = 0; row < decomp.steps; ++row) {
      decomp.eigenvectors[static_cast<std::size_t>(row)][static_cast<std::size_t>(sorted)] =
          raw_eigenvectors[static_cast<std::size_t>(row)][static_cast<std::size_t>(original)];
    }
  }
  return decomp;
}

void for_each_lanczos_vector(
    const MomentumBlock& block,
    const LanczosDecomposition& decomp,
    const std::vector<Complex>& initial,
    const std::function<void(int, const std::vector<Complex>&)>& callback) {
  if (decomp.steps == 0) {
    return;
  }
  std::vector<Complex> q_prev(block.basis_dim, 0.0);
  std::vector<Complex> q_cur = initial;
  std::vector<Complex> q_next(block.basis_dim, 0.0);
  std::vector<Complex> hq(block.basis_dim, 0.0);
  double beta_prev = 0.0;
  for (int step = 0; step < decomp.steps; ++step) {
    callback(step, q_cur);
    if (step + 1 == decomp.steps) {
      break;
    }
    apply_hamiltonian(block, q_cur, hq);
    for (std::size_t i = 0; i < block.basis_dim; ++i) {
      q_next[i] = hq[i] - decomp.alpha[static_cast<std::size_t>(step)] * q_cur[i] -
                  beta_prev * q_prev[i];
    }
    const double beta = decomp.beta[static_cast<std::size_t>(step)];
    q_prev.swap(q_cur);
    q_cur.swap(q_next);
    for (Complex& value : q_cur) {
      value /= beta;
    }
    beta_prev = beta;
  }
}

std::vector<Complex> project_current_store_bases(
    const MomentumBlock& current_block,
    const LanczosDecomposition& first,
    const LanczosDecomposition& second) {
  const int ma = first.steps;
  const int mb = second.steps;
  std::vector<Complex> krylov_matrix(static_cast<std::size_t>(mb) * ma, 0.0);
  std::vector<Complex> current_q(current_block.basis_dim, 0.0);
  for (int ia = 0; ia < ma; ++ia) {
    apply_operator(current_block, first.basis[static_cast<std::size_t>(ia)], current_q);
    for (int jb = 0; jb < mb; ++jb) {
      krylov_matrix[static_cast<std::size_t>(jb) + static_cast<std::size_t>(ia) * mb] =
          dot(second.basis[static_cast<std::size_t>(jb)], current_q);
    }
  }

  std::vector<Complex> temp(static_cast<std::size_t>(mb) * ma, 0.0);
  for (int lb = 0; lb < mb; ++lb) {
    for (int ia = 0; ia < ma; ++ia) {
      Complex sum = 0.0;
      for (int jb = 0; jb < mb; ++jb) {
        sum += second.eigenvectors[static_cast<std::size_t>(jb)][static_cast<std::size_t>(lb)] *
               krylov_matrix[static_cast<std::size_t>(jb) + static_cast<std::size_t>(ia) * mb];
      }
      temp[static_cast<std::size_t>(lb) + static_cast<std::size_t>(ia) * mb] = sum;
    }
  }

  std::vector<Complex> projected(static_cast<std::size_t>(mb) * ma, 0.0);
  for (int lb = 0; lb < mb; ++lb) {
    for (int la = 0; la < ma; ++la) {
      Complex sum = 0.0;
      for (int ia = 0; ia < ma; ++ia) {
        sum += temp[static_cast<std::size_t>(lb) + static_cast<std::size_t>(ia) * mb] *
               first.eigenvectors[static_cast<std::size_t>(ia)][static_cast<std::size_t>(la)];
      }
      projected[static_cast<std::size_t>(lb) + static_cast<std::size_t>(la) * mb] = sum;
    }
  }
  return projected;
}

std::vector<Complex> project_current_recompute_bases(
    const MomentumBlock& h_block,
    const MomentumBlock& current_block,
    const LanczosDecomposition& first,
    const LanczosDecomposition& second,
    const std::vector<Complex>& first_initial,
    const std::vector<Complex>& second_initial) {
  const int ma = first.steps;
  const int mb = second.steps;
  std::vector<Complex> projected(static_cast<std::size_t>(mb) * ma, 0.0);
  std::vector<Complex> current_q(current_block.basis_dim, 0.0);

  for_each_lanczos_vector(
      h_block,
      first,
      first_initial,
      [&](int ia, const std::vector<Complex>& qa) {
        apply_operator(current_block, qa, current_q);
        for_each_lanczos_vector(
            h_block,
            second,
            second_initial,
            [&](int jb, const std::vector<Complex>& qb) {
              const Complex kji = dot(qb, current_q);
              if (std::norm(kji) < 1e-32) {
                return;
              }
              for (int lb = 0; lb < mb; ++lb) {
                const double ub = second.eigenvectors[static_cast<std::size_t>(jb)][static_cast<std::size_t>(lb)];
                if (std::abs(ub) < 1e-15) {
                  continue;
                }
                for (int la = 0; la < ma; ++la) {
                  const double ua = first.eigenvectors[static_cast<std::size_t>(ia)][static_cast<std::size_t>(la)];
                  if (std::abs(ua) < 1e-15) {
                    continue;
                  }
                  projected[static_cast<std::size_t>(lb) + static_cast<std::size_t>(la) * mb] +=
                      ub * ua * kji;
                }
              }
            });
      });
  return projected;
}

std::vector<Complex> random_initial_vector(std::size_t dim, std::uint64_t seed) {
  std::mt19937_64 rng(mix_u64(seed));
  std::normal_distribution<double> normal(0.0, 1.0);
  std::vector<Complex> vec(dim, 0.0);
  for (Complex& value : vec) {
    value = Complex(normal(rng), normal(rng));
  }
  normalize(vec);
  return vec;
}

double estimated_store_basis_mb(std::size_t basis_dim, int krylov_dim) {
  const double vectors = 2.0 * static_cast<double>(krylov_dim) + 6.0;
  return vectors * static_cast<double>(basis_dim) * static_cast<double>(sizeof(Complex)) /
         (1024.0 * 1024.0);
}

bool should_store_bases(
    const ConductivityParams& params,
    std::size_t basis_dim,
    int krylov_dim) {
  if (params.memory_mode == ConductivityMemoryMode::StoreBases) {
    return true;
  }
  if (params.memory_mode == ConductivityMemoryMode::Recompute) {
    return false;
  }
  return estimated_store_basis_mb(basis_dim, krylov_dim) <= params.basis_memory_mb;
}

ConductivitySampleRecord run_conductivity_sample(
    const MomentumBlock& h_block,
    const MomentumBlock& current_block,
    const Params& base,
    const ConductivitySampleKey& sample_key,
    bool store_bases) {
  ConductivitySampleRecord record;
  record.key = sample_key;
  record.particles = h_block.particles;
  record.basis_dim = h_block.basis_dim;
  record.trace_prefactor = static_cast<double>(h_block.basis_dim) /
                           static_cast<double>(base.samples);

  const int krylov_dim = std::min<int>(base.lanczos_steps, static_cast<int>(h_block.basis_dim));
  const std::uint64_t block_seed = make_block_seed(
      base,
      sample_key.block.n_up,
      sample_key.block.n_down,
      sample_key.block.mx,
      sample_key.block.my);
  const std::vector<Complex> r = random_initial_vector(
      h_block.basis_dim,
      block_seed ^ static_cast<std::uint64_t>(sample_key.sample));

  LanczosDecomposition first = run_lanczos_decomposition(h_block, r, krylov_dim, store_bases);
  record.first_eigenvalues = first.eigenvalues;
  record.first_v0 = first.first_components;

  std::vector<Complex> jr(h_block.basis_dim, 0.0);
  apply_operator(current_block, r, jr);
  record.current_norm = std::sqrt(complex_norm2(jr));
  if (record.current_norm < 1e-14) {
    return record;
  }
  for (Complex& value : jr) {
    value /= record.current_norm;
  }

  LanczosDecomposition second = run_lanczos_decomposition(h_block, jr, krylov_dim, store_bases);
  record.second_eigenvalues = second.eigenvalues;
  record.second_v0 = second.first_components;
  if (store_bases) {
    record.current_matrix = project_current_store_bases(current_block, first, second);
  } else {
    record.current_matrix = project_current_recompute_bases(
        h_block, current_block, first, second, r, jr);
  }
  return record;
}

std::string beta_list_text(const std::vector<double>& values) {
  std::ostringstream out;
  out << std::setprecision(17);
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i > 0) {
      out << ",";
    }
    out << values[i];
  }
  return out.str();
}

std::string conductivity_checkpoint_metadata_text(const ConductivityParams& params) {
  const Params& base = params.base;
  std::ostringstream out;
  out << std::setprecision(17)
      << "lx=" << base.lx << "\n"
      << "ly=" << base.ly << "\n"
      << "tx=" << base.tx << "\n"
      << "ty=" << base.ty << "\n"
      << "tp=" << base.tp << "\n"
      << "phix=" << base.phix << "\n"
      << "phiy=" << base.phiy << "\n"
      << "u=" << base.u << "\n"
      << "v=" << base.v << "\n"
      << "samples=" << base.samples << "\n"
      << "lanczos_steps=" << base.lanczos_steps << "\n"
      << "exact_block_threshold=" << base.exact_block_threshold << "\n"
      << "seed=" << base.seed << "\n"
      << "direction=" << params.direction << "\n";
  return out.str();
}

void validate_or_write_conductivity_checkpoint_metadata(const ConductivityParams& params) {
  if (params.conductivity_checkpoint.empty()) {
    return;
  }
  const std::string path = params.conductivity_checkpoint + ".meta";
  const std::string expected = conductivity_checkpoint_metadata_text(params);
  {
    std::ifstream in(path);
    if (in) {
      const std::string existing(
          (std::istreambuf_iterator<char>(in)),
          std::istreambuf_iterator<char>());
      if (existing != expected) {
        die("Conductivity checkpoint metadata does not match current run: " + path);
      }
      return;
    }
  }
  std::ofstream out(path);
  if (!out) {
    die("Failed to write conductivity checkpoint metadata: " + path);
  }
  out << expected;
}

void ensure_conductivity_checkpoint_header(const std::string& path) {
  if (path.empty()) {
    return;
  }
  {
    std::ifstream in(path, std::ios::binary);
    if (in && in.peek() != std::ifstream::traits_type::eof()) {
      return;
    }
  }
  std::ofstream out(path, std::ios::binary);
  if (!out) {
    die("Failed to create conductivity checkpoint file: " + path);
  }
  const char magic[8] = {'F', 'T', 'L', 'M', 'C', 'C', '1', '\n'};
  out.write(magic, sizeof(magic));
}

template <typename T>
void write_vector_binary(std::ostream& out, const std::vector<T>& values) {
  if (!values.empty()) {
    out.write(reinterpret_cast<const char*>(values.data()),
              static_cast<std::streamsize>(sizeof(T) * values.size()));
  }
}

template <typename T>
bool read_vector_binary(std::istream& in, std::vector<T>& values) {
  if (values.empty()) {
    return true;
  }
  return static_cast<bool>(in.read(
      reinterpret_cast<char*>(values.data()),
      static_cast<std::streamsize>(sizeof(T) * values.size())));
}

void append_conductivity_checkpoint_record(
    const std::string& path,
    const ConductivitySampleRecord& record) {
  if (path.empty()) {
    return;
  }
  ensure_conductivity_checkpoint_header(path);
  std::ofstream out(path, std::ios::binary | std::ios::app);
  if (!out) {
    die("Failed to append conductivity checkpoint file: " + path);
  }
  const char record_magic[8] = {'C', 'C', 'S', 'R', 'E', 'C', '1', '\n'};
  out.write(record_magic, sizeof(record_magic));
  write_binary(out, record.key.block.n_up);
  write_binary(out, record.key.block.n_down);
  write_binary(out, record.key.block.mx);
  write_binary(out, record.key.block.my);
  write_binary(out, record.key.sample);
  write_binary(out, record.particles);
  const std::uint64_t basis_dim = static_cast<std::uint64_t>(record.basis_dim);
  write_binary(out, basis_dim);
  write_binary(out, record.trace_prefactor);
  write_binary(out, record.current_norm);
  const std::uint64_t first_count = static_cast<std::uint64_t>(record.first_eigenvalues.size());
  const std::uint64_t second_count = static_cast<std::uint64_t>(record.second_eigenvalues.size());
  const std::uint64_t matrix_count = static_cast<std::uint64_t>(record.current_matrix.size());
  write_binary(out, first_count);
  write_binary(out, second_count);
  write_binary(out, matrix_count);
  write_vector_binary(out, record.first_eigenvalues);
  write_vector_binary(out, record.first_v0);
  write_vector_binary(out, record.second_eigenvalues);
  write_vector_binary(out, record.second_v0);
  write_vector_binary(out, record.current_matrix);
  out.flush();
}

bool read_conductivity_record_body(std::istream& in, ConductivitySampleRecord& record) {
  std::uint64_t basis_dim = 0;
  std::uint64_t first_count = 0;
  std::uint64_t second_count = 0;
  std::uint64_t matrix_count = 0;
  if (!read_binary(in, record.key.block.n_up) ||
      !read_binary(in, record.key.block.n_down) ||
      !read_binary(in, record.key.block.mx) ||
      !read_binary(in, record.key.block.my) ||
      !read_binary(in, record.key.sample) ||
      !read_binary(in, record.particles) ||
      !read_binary(in, basis_dim) ||
      !read_binary(in, record.trace_prefactor) ||
      !read_binary(in, record.current_norm) ||
      !read_binary(in, first_count) ||
      !read_binary(in, second_count) ||
      !read_binary(in, matrix_count)) {
    return false;
  }
  if (first_count > 1000000ULL || second_count > 1000000ULL ||
      matrix_count > 100000000ULL) {
    die("Conductivity checkpoint record is implausibly large.");
  }
  if (second_count > 0 && matrix_count != first_count * second_count) {
    die("Conductivity checkpoint record has inconsistent matrix dimensions.");
  }
  if (second_count == 0 && matrix_count != 0) {
    die("Conductivity checkpoint record has a matrix but no second Lanczos space.");
  }
  record.basis_dim = static_cast<std::size_t>(basis_dim);
  record.first_eigenvalues.resize(static_cast<std::size_t>(first_count));
  record.first_v0.resize(static_cast<std::size_t>(first_count));
  record.second_eigenvalues.resize(static_cast<std::size_t>(second_count));
  record.second_v0.resize(static_cast<std::size_t>(second_count));
  record.current_matrix.resize(static_cast<std::size_t>(matrix_count));
  return read_vector_binary(in, record.first_eigenvalues) &&
         read_vector_binary(in, record.first_v0) &&
         read_vector_binary(in, record.second_eigenvalues) &&
         read_vector_binary(in, record.second_v0) &&
         read_vector_binary(in, record.current_matrix);
}

void for_each_checkpoint_record(
    const std::string& path,
    const std::function<void(const ConductivitySampleRecord&)>& callback) {
  if (path.empty()) {
    return;
  }
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return;
  }
  char magic[8] = {};
  if (!in.read(magic, sizeof(magic))) {
    return;
  }
  const char expected_magic[8] = {'F', 'T', 'L', 'M', 'C', 'C', '1', '\n'};
  if (std::memcmp(magic, expected_magic, sizeof(magic)) != 0) {
    die("Conductivity checkpoint has an unrecognized format: " + path);
  }

  const char expected_record[8] = {'C', 'C', 'S', 'R', 'E', 'C', '1', '\n'};
  while (true) {
    char record_magic[8] = {};
    if (!in.read(record_magic, sizeof(record_magic))) {
      break;
    }
    if (std::memcmp(record_magic, expected_record, sizeof(record_magic)) != 0) {
      break;
    }
    ConductivitySampleRecord record;
    if (!read_conductivity_record_body(in, record)) {
      break;
    }
    callback(record);
  }
}

std::unordered_set<ConductivitySampleKey, ConductivitySampleKeyHash>
read_completed_conductivity_samples(const std::string& path) {
  std::unordered_set<ConductivitySampleKey, ConductivitySampleKeyHash> completed;
  for_each_checkpoint_record(path, [&](const ConductivitySampleRecord& record) {
    completed.insert(record.key);
  });
  return completed;
}

void for_each_all_sample_record(
    const std::string& checkpoint,
    const std::vector<ConductivitySampleRecord>& memory_records,
    const std::function<void(const ConductivitySampleRecord&)>& callback) {
  if (!checkpoint.empty()) {
    for_each_checkpoint_record(checkpoint, callback);
  } else {
    for (const ConductivitySampleRecord& record : memory_records) {
      callback(record);
    }
  }
}

void add_exact_thermo_logs(
    const std::vector<ExactConductivityBlock>& exact_blocks,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    std::vector<std::vector<LogAccumulator>>& logz) {
  for (const ExactConductivityBlock& block : exact_blocks) {
    const double particles = static_cast<double>(block.particles);
    for (double energy : block.eigenvalues) {
      for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
        const double beta = beta_values[ibeta];
        for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
          const double mu = mu_values[imu];
          logz[ibeta][imu].add_log(-beta * (energy - mu * particles));
        }
      }
    }
  }
}

void add_record_thermo_logs(
    const ConductivitySampleRecord& record,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    std::vector<std::vector<LogAccumulator>>& logz) {
  const double particles = static_cast<double>(record.particles);
  for (std::size_t state = 0; state < record.first_eigenvalues.size(); ++state) {
    const double prefactor = record.trace_prefactor * record.first_v0[state] * record.first_v0[state];
    if (prefactor <= 0.0) {
      continue;
    }
    const double log_prefactor = std::log(prefactor);
    const double energy = record.first_eigenvalues[state];
    for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
      const double beta = beta_values[ibeta];
      for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
        const double mu = mu_values[imu];
        logz[ibeta][imu].add_log(log_prefactor - beta * (energy - mu * particles));
      }
    }
  }
}

std::vector<std::vector<ThermoPoint>> compute_thermo_points(
    const std::vector<ExactConductivityBlock>& exact_blocks,
    const std::string& checkpoint,
    const std::vector<ConductivitySampleRecord>& memory_records,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    int sites) {
  std::vector<std::vector<LogAccumulator>> logz(
      beta_values.size(), std::vector<LogAccumulator>(mu_values.size()));
  add_exact_thermo_logs(exact_blocks, beta_values, mu_values, logz);
  for_each_all_sample_record(
      checkpoint,
      memory_records,
      [&](const ConductivitySampleRecord& record) {
        add_record_thermo_logs(record, beta_values, mu_values, logz);
      });

  std::vector<std::vector<ThermoPoint>> thermo(
      beta_values.size(), std::vector<ThermoPoint>(mu_values.size()));
  std::vector<std::vector<double>> n_total(
      beta_values.size(), std::vector<double>(mu_values.size(), 0.0));
  std::vector<std::vector<double>> n2_total(
      beta_values.size(), std::vector<double>(mu_values.size(), 0.0));

  for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
    for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
      if (!std::isfinite(logz[ibeta][imu].log_sum)) {
        die("No thermodynamic weight accumulated for beta/mu point.");
      }
      thermo[ibeta][imu].log_partition = logz[ibeta][imu].log_sum;
    }
  }

  auto add_weight = [&](double energy, double prefactor, int particles) {
    if (prefactor <= 0.0) {
      return;
    }
    const double log_prefactor = std::log(prefactor);
    const double p = static_cast<double>(particles);
    for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
      const double beta = beta_values[ibeta];
      for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
        const double mu = mu_values[imu];
        const double logw = log_prefactor - beta * (energy - mu * p) -
                            thermo[ibeta][imu].log_partition;
        const double weight = std::exp(logw);
        n_total[ibeta][imu] += weight * p;
        n2_total[ibeta][imu] += weight * p * p;
      }
    }
  };

  for (const ExactConductivityBlock& block : exact_blocks) {
    for (double energy : block.eigenvalues) {
      add_weight(energy, 1.0, block.particles);
    }
  }
  for_each_all_sample_record(
      checkpoint,
      memory_records,
      [&](const ConductivitySampleRecord& record) {
        for (std::size_t state = 0; state < record.first_eigenvalues.size(); ++state) {
          add_weight(
              record.first_eigenvalues[state],
              record.trace_prefactor * record.first_v0[state] * record.first_v0[state],
              record.particles);
        }
      });

  for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
    const double beta = beta_values[ibeta];
    for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
      ThermoPoint& point = thermo[ibeta][imu];
      point.density = n_total[ibeta][imu] / static_cast<double>(sites);
      point.charge_correlation = std::max(
          0.0,
          (n2_total[ibeta][imu] - n_total[ibeta][imu] * n_total[ibeta][imu]) /
              static_cast<double>(sites));
      point.compressibility = beta * point.charge_correlation;
    }
  }
  return thermo;
}

void add_lorentzian_contribution(
    double amplitude,
    double delta_e,
    double thermal_factor,
    const std::vector<double>& omega_values,
    double eta,
    int sites,
    std::vector<double>& sigma_row,
    double& dc_sigma) {
  const double scale = thermal_factor * amplitude / static_cast<double>(sites);
  if (std::abs(scale) < 1e-300) {
    return;
  }
  for (std::size_t iw = 0; iw < omega_values.size(); ++iw) {
    const double diff = omega_values[iw] - delta_e;
    const double broadened_pi_delta = eta / (diff * diff + eta * eta);
    sigma_row[iw] += scale * broadened_pi_delta;
  }
  const double dc_diff = -delta_e;
  dc_sigma += scale * eta / (dc_diff * dc_diff + eta * eta);
}

void accumulate_exact_sigma(
    const ExactConductivityBlock& block,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    const std::vector<double>& omega_values,
    const std::vector<std::vector<ThermoPoint>>& thermo,
    double eta,
    int sites,
    std::vector<std::vector<std::vector<double>>>& sigma,
    std::vector<std::vector<double>>& dc_sigma) {
  constexpr double dE_eps = 1e-10;
  const int n = static_cast<int>(block.basis_dim);
  const double particles = static_cast<double>(block.particles);
  for (int initial = 0; initial < n; ++initial) {
    const double ei = block.eigenvalues[static_cast<std::size_t>(initial)];
    for (int final = 0; final < n; ++final) {
      const double ej = block.eigenvalues[static_cast<std::size_t>(final)];
      const double delta_e = ej - ei;
      if (delta_e < -dE_eps) {
        continue;
      }
      const double j2 = block.current_abs2[static_cast<std::size_t>(final) +
                                           static_cast<std::size_t>(initial) * n];
      if (j2 <= 0.0) {
        continue;
      }
      for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
        const double beta = beta_values[ibeta];
        for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
          const double mu = mu_values[imu];
          const double log_z = thermo[ibeta][imu].log_partition;
          const double wi = std::exp(-beta * (ei - mu * particles) - log_z);
          double thermal_factor = 0.0;
          if (std::abs(delta_e) <= dE_eps) {
            thermal_factor = beta * wi;
          } else {
            const double wj = std::exp(-beta * (ej - mu * particles) - log_z);
            thermal_factor = (wi - wj) / delta_e;
          }
          add_lorentzian_contribution(
              j2,
              delta_e,
              thermal_factor,
              omega_values,
              eta,
              sites,
              sigma[ibeta][imu],
              dc_sigma[ibeta][imu]);
        }
      }
    }
  }
}

void accumulate_sample_sigma(
    const ConductivitySampleRecord& record,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    const std::vector<double>& omega_values,
    const std::vector<std::vector<ThermoPoint>>& thermo,
    double eta,
    int sites,
    std::vector<std::vector<std::vector<double>>>& sigma,
    std::vector<std::vector<double>>& dc_sigma) {
  constexpr double dE_eps = 1e-10;
  const std::size_t ma = record.first_eigenvalues.size();
  const std::size_t mb = record.second_eigenvalues.size();
  if (ma == 0 || mb == 0 || record.current_matrix.empty() || record.current_norm <= 0.0) {
    return;
  }
  const double particles = static_cast<double>(record.particles);
  for (std::size_t ia = 0; ia < ma; ++ia) {
    const double ea = record.first_eigenvalues[ia];
    const double va = record.first_v0[ia];
    for (std::size_t ib = 0; ib < mb; ++ib) {
      const double eb = record.second_eigenvalues[ib];
      const double delta_e = eb - ea;
      if (delta_e < -dE_eps) {
        continue;
      }
      const Complex matrix_element = record.current_matrix[ib + ia * mb];
      const double amplitude = record.trace_prefactor * record.current_norm *
          (va * record.second_v0[ib] * std::conj(matrix_element)).real();
      if (std::abs(amplitude) < 1e-18) {
        continue;
      }
      for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
        const double beta = beta_values[ibeta];
        for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
          const double mu = mu_values[imu];
          const double log_z = thermo[ibeta][imu].log_partition;
          const double wa = std::exp(-beta * (ea - mu * particles) - log_z);
          double thermal_factor = 0.0;
          if (std::abs(delta_e) <= dE_eps) {
            thermal_factor = beta * wa;
          } else {
            const double wb = std::exp(-beta * (eb - mu * particles) - log_z);
            thermal_factor = (wa - wb) / delta_e;
          }
          add_lorentzian_contribution(
              amplitude,
              delta_e,
              thermal_factor,
              omega_values,
              eta,
              sites,
              sigma[ibeta][imu],
              dc_sigma[ibeta][imu]);
        }
      }
    }
  }
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
    const std::vector<double> mu_values = linspace(base.mu_min, base.mu_max, base.mu_count);
    const std::vector<double> omega_values =
        linspace(params.omega_min, params.omega_max, params.omega_count);

    validate_or_write_conductivity_checkpoint_metadata(params);
    ensure_conductivity_checkpoint_header(params.conductivity_checkpoint);
    auto completed_samples = read_completed_conductivity_samples(params.conductivity_checkpoint);

    std::cout << "Running two-Lanczos FTLM optical/DC conductivity on a "
              << base.lx << "x" << base.ly << " rectangular lattice ("
              << lattice.sites << " sites)\n";
    std::cout << "beta=" << beta_list_text(base.beta_values)
              << " U=" << base.u << " V=" << base.v
              << " tx=" << base.tx << " ty=" << base.ty << " tp=" << base.tp
              << " phix=" << base.phix << " phiy=" << base.phiy
              << " direction=" << params.direction
              << " eta=" << params.eta
              << " omega=[" << params.omega_min << "," << params.omega_max
              << "] count=" << params.omega_count
              << " samples=" << base.samples
              << " lanczos_steps=" << base.lanczos_steps
              << " exact_block_threshold=" << base.exact_block_threshold
              << " memory_mode=" << memory_mode_name(params.memory_mode)
              << " basis_memory_mb=" << params.basis_memory_mb << "\n";
    if (!params.conductivity_checkpoint.empty()) {
      std::cout << "conductivity_checkpoint=" << params.conductivity_checkpoint
                << " completed_samples=" << completed_samples.size() << "\n";
    }

    const auto start = std::chrono::steady_clock::now();
    const auto max_runtime = std::chrono::duration<double, std::ratio<60>>(
        params.max_runtime_minutes);
    bool stopped_for_runtime = false;
    std::vector<ExactConductivityBlock> exact_blocks;
    std::vector<ConductivitySampleRecord> memory_records;

    for (int n_up = 0; n_up <= lattice.sites && !stopped_for_runtime; ++n_up) {
      for (int n_down = 0; n_down <= lattice.sites && !stopped_for_runtime; ++n_down) {
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
        for (int my = 0; my < lattice.ly && !stopped_for_runtime; ++my) {
          for (int mx = 0; mx < lattice.lx && !stopped_for_runtime; ++mx) {
            const MomentumBlock h_block = build_momentum_block(sector, lattice, mx, my);
            active_sum += h_block.basis_dim;
            if (h_block.basis_dim == 0) {
              continue;
            }
            const BlockKey block_key{n_up, n_down, mx, my};
            const MomentumBlock j_block =
                build_current_block(sector, lattice, mx, my, params.direction);
            const bool use_exact = h_block.basis_dim <= base.exact_block_threshold;
            if (use_exact) {
              exact_blocks.push_back(diagonalize_exact_conductivity_block(h_block, j_block));
              continue;
            }

            const int krylov_dim = std::min<int>(base.lanczos_steps, static_cast<int>(h_block.basis_dim));
            const bool store_bases = should_store_bases(params, h_block.basis_dim, krylov_dim);
            std::cout << "  FTLM block Nup=" << n_up
                      << " Ndown=" << n_down
                      << " mx=" << mx
                      << " my=" << my
                      << " basis_dim=" << h_block.basis_dim
                      << " krylov_dim=" << krylov_dim
                      << " basis_strategy=" << (store_bases ? "store-bases" : "recompute")
                      << " est_store_mb=" << estimated_store_basis_mb(h_block.basis_dim, krylov_dim)
                      << "\n";

            for (int sample = 0; sample < base.samples; ++sample) {
              const ConductivitySampleKey sample_key{block_key, sample};
              if (completed_samples.find(sample_key) != completed_samples.end()) {
                continue;
              }
              ConductivitySampleRecord record = run_conductivity_sample(
                  h_block, j_block, base, sample_key, store_bases);
              append_conductivity_checkpoint_record(params.conductivity_checkpoint, record);
              if (params.conductivity_checkpoint.empty()) {
                memory_records.push_back(std::move(record));
              } else {
                completed_samples.insert(sample_key);
              }
              if (params.max_runtime_minutes > 0.0) {
                const auto elapsed = std::chrono::steady_clock::now() - start;
                if (elapsed >= max_runtime) {
                  stopped_for_runtime = true;
                  break;
                }
              }
            }
          }
        }
        std::cout << "  sector Nup=" << n_up << " Ndown=" << n_down
                  << " full_dim=" << sector.full_dim
                  << " parents=" << sector.parents.size()
                  << " active_k_dim_sum=" << active_sum << "\n";
      }
    }

    if (stopped_for_runtime) {
      std::cout << "Reached --max-runtime-minutes after a completed sample; checkpoint flushed. "
                << "Rerun the same command to resume and write final CSVs.\n";
      return 0;
    }

    const std::vector<std::vector<ThermoPoint>> thermo = compute_thermo_points(
        exact_blocks,
        params.conductivity_checkpoint,
        memory_records,
        base.beta_values,
        mu_values,
        lattice.sites);

    std::vector<std::vector<std::vector<double>>> sigma(
        base.beta_values.size(),
        std::vector<std::vector<double>>(
            mu_values.size(), std::vector<double>(omega_values.size(), 0.0)));
    std::vector<std::vector<double>> dc_sigma(
        base.beta_values.size(), std::vector<double>(mu_values.size(), 0.0));

    for (const ExactConductivityBlock& block : exact_blocks) {
      accumulate_exact_sigma(
          block,
          base.beta_values,
          mu_values,
          omega_values,
          thermo,
          params.eta,
          lattice.sites,
          sigma,
          dc_sigma);
    }
    for_each_all_sample_record(
        params.conductivity_checkpoint,
        memory_records,
        [&](const ConductivitySampleRecord& record) {
          accumulate_sample_sigma(
              record,
              base.beta_values,
              mu_values,
              omega_values,
              thermo,
              params.eta,
              lattice.sites,
              sigma,
              dc_sigma);
        });

    write_conductivity_csv(base.output, base.beta_values, mu_values, omega_values, thermo, sigma);
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
