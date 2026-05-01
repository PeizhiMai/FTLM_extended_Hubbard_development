#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace {

using Complex = std::complex<double>;

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kEps = 1e-12;

struct Params {
  int lx = 2;
  int ly = 2;
  double tx = 1.0;
  double ty = 1.0;
  double tp = 0.0;
  double phix = 1.0;
  double phiy = 1.0;
  double u = 8.0;
  double v = 0.0;
  double beta = 2.0;
  std::vector<double> beta_values;
  double mu_min = -4.0;
  double mu_max = 10.0;
  int mu_count = 61;
  int samples = 5;
  int lanczos_steps = 80;
  int threads = 0;
  std::size_t exact_block_threshold = 256;
  std::uint64_t seed = 12345;
  std::size_t max_sector_dim = 2000000;
  std::string output = "n_vs_mu.csv";
  std::string checkpoint;
  std::string trace_partition_csv;
  std::string debug_block_csv;
  int debug_block_nup = -1;
  int debug_block_ndown = -1;
  int debug_block_mx = -1;
  int debug_block_my = -1;
};

struct StateKey {
  std::uint64_t up = 0;
  std::uint64_t down = 0;

  bool operator==(const StateKey& other) const {
    return up == other.up && down == other.down;
  }
};

struct StateKeyHash {
  std::size_t operator()(const StateKey& key) const {
    const std::size_t h1 = std::hash<std::uint64_t>{}(key.up);
    const std::size_t h2 = std::hash<std::uint64_t>{}(key.down);
    return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6U) + (h1 >> 2U));
  }
};

struct Relation {
  int parent = -1;
  int shift = -1;
  int sign = 1;
};

struct CompactHop {
  int parent = -1;
  int shift = -1;
  int sign = 1;
  Complex amplitude = 0.0;
};

struct ParentData {
  StateKey representative;
  int degeneracy = 0;
  double diagonal = 0.0;
  std::vector<CompactHop> hops;
};

struct SectorBasis {
  int n_up = 0;
  int n_down = 0;
  int particles = 0;
  std::size_t full_dim = 0;
  std::vector<ParentData> parents;
};

struct Edge {
  int target = -1;
  Complex amplitude = 0.0;
};

struct MomentumBlock {
  int mx = 0;
  int my = 0;
  int particles = 0;
  std::size_t basis_dim = 0;
  std::vector<double> diagonal;
  std::vector<int> row_ptr;
  std::vector<int> col_idx;
  std::vector<Complex> values;
};

struct LanczosSpectrum {
  std::vector<double> eigenvalues;
  std::vector<double> sample_weights;
  std::size_t basis_dim = 0;
  int particles = 0;
};

struct BlockThermo {
  double z = 0.0;
  double n_total = 0.0;
  double n2_total = 0.0;
  double min_shifted_energy = std::numeric_limits<double>::infinity();
};

struct StoredSpectrum {
  std::vector<double> eigenvalues;
  std::vector<double> sample_weights;
  double trace_prefactor = 1.0;
  std::size_t basis_dim = 0;
  int particles = 0;
  bool exact = false;
};

struct BlockKey {
  int n_up = 0;
  int n_down = 0;
  int mx = 0;
  int my = 0;

  bool operator==(const BlockKey& other) const {
    return n_up == other.n_up && n_down == other.n_down &&
           mx == other.mx && my == other.my;
  }
};

struct BlockKeyHash {
  std::size_t operator()(const BlockKey& key) const {
    std::size_t value = std::hash<int>{}(key.n_up);
    value ^= std::hash<int>{}(key.n_down) + 0x9e3779b9U + (value << 6U) + (value >> 2U);
    value ^= std::hash<int>{}(key.mx) + 0x9e3779b9U + (value << 6U) + (value >> 2U);
    value ^= std::hash<int>{}(key.my) + 0x9e3779b9U + (value << 6U) + (value >> 2U);
    return value;
  }
};

struct Lattice {
  int lx = 0;
  int ly = 0;
  int sites = 0;
  std::vector<std::array<int, 8>> neighbors;
  std::vector<std::pair<int, int>> unique_bonds;
};

[[noreturn]] void die(const std::string& message) {
  throw std::runtime_error(message);
}

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

template <typename T>
T parse_value(const std::string& text) {
  std::istringstream in(text);
  T value{};
  in >> value;
  if (!in || !in.eof()) {
    die("Failed to parse value: " + text);
  }
  return value;
}

std::vector<double> parse_double_list(const std::string& text) {
  std::vector<double> values;
  std::istringstream in(text);
  std::string item;
  while (std::getline(in, item, ',')) {
    if (item.empty()) {
      die("Empty entry in comma-separated value list: " + text);
    }
    values.push_back(parse_value<double>(item));
  }
  if (values.empty()) {
    die("Expected at least one value in comma-separated list: " + text);
  }
  return values;
}

void print_help(const char* argv0) {
  std::cout
      << "Usage: " << argv0 << " [options]\n"
      << "  --lx N                lattice size in x\n"
      << "  --ly N                lattice size in y\n"
      << "  --tx X                hopping along x bonds\n"
      << "  --ty Y                hopping along y bonds\n"
      << "  --tp X                second-neighbor diagonal hopping\n"
      << "  --phix X              homogeneous x twist in units of 2*pi\n"
      << "  --phiy Y              homogeneous y twist in units of 2*pi\n"
      << "  --u U                 onsite interaction\n"
      << "  --v V                 nearest-neighbor interaction\n"
      << "  --beta B              inverse temperature\n"
      << "  --beta-list LIST      comma-separated inverse temperatures\n"
      << "  --mu-min X            minimum chemical potential\n"
      << "  --mu-max X            maximum chemical potential\n"
      << "  --mu-count N          number of mu points\n"
      << "  --samples N           FTLM random vectors per k block\n"
      << "  --lanczos-steps N     Lanczos steps per random vector\n"
      << "  --threads N           OpenMP threads (0 uses runtime default)\n"
      << "  --exact-block-threshold N\n"
      << "                        use exact diagonalization for k blocks with dimension <= N\n"
      << "  --seed N              random seed\n"
      << "  --max-sector-dim N    abort if a full sector exceeds this size\n"
      << "  --checkpoint PATH     append/resume compact per-block spectra\n"
      << "  --trace-partition-csv PATH\n"
      << "                        write per-block FTLM vs ED thermodynamic traces\n"
      << "  --debug-block-nup N   selected block debug Nup\n"
      << "  --debug-block-ndown N selected block debug Ndown\n"
      << "  --debug-block-mx N    selected block debug mx\n"
      << "  --debug-block-my N    selected block debug my\n"
      << "  --debug-block-csv PATH\n"
      << "                        write selected-block FTLM vs ED comparison\n"
      << "  --output PATH         CSV output path\n"
      << "  --help                show this message\n";
}

Params parse_args(int argc, char** argv) {
  Params p;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto next = [&](const std::string& name) {
      if (i + 1 >= argc) {
        die("Missing value for " + name);
      }
      return std::string(argv[++i]);
    };

    if (arg == "--help") {
      print_help(argv[0]);
      std::exit(0);
    } else if (arg == "--lx") {
      p.lx = parse_value<int>(next(arg));
    } else if (arg == "--ly") {
      p.ly = parse_value<int>(next(arg));
    } else if (arg == "--tx") {
      p.tx = parse_value<double>(next(arg));
    } else if (arg == "--ty") {
      p.ty = parse_value<double>(next(arg));
    } else if (arg == "--tp") {
      p.tp = parse_value<double>(next(arg));
    } else if (arg == "--phix") {
      p.phix = parse_value<double>(next(arg));
    } else if (arg == "--phiy") {
      p.phiy = parse_value<double>(next(arg));
    } else if (arg == "--u") {
      p.u = parse_value<double>(next(arg));
    } else if (arg == "--v") {
      p.v = parse_value<double>(next(arg));
    } else if (arg == "--beta") {
      p.beta = parse_value<double>(next(arg));
    } else if (arg == "--beta-list") {
      p.beta_values = parse_double_list(next(arg));
    } else if (arg == "--mu-min") {
      p.mu_min = parse_value<double>(next(arg));
    } else if (arg == "--mu-max") {
      p.mu_max = parse_value<double>(next(arg));
    } else if (arg == "--mu-count") {
      p.mu_count = parse_value<int>(next(arg));
    } else if (arg == "--samples") {
      p.samples = parse_value<int>(next(arg));
    } else if (arg == "--lanczos-steps") {
      p.lanczos_steps = parse_value<int>(next(arg));
    } else if (arg == "--threads") {
      p.threads = parse_value<int>(next(arg));
    } else if (arg == "--exact-block-threshold") {
      p.exact_block_threshold = parse_value<std::size_t>(next(arg));
    } else if (arg == "--seed") {
      p.seed = parse_value<std::uint64_t>(next(arg));
    } else if (arg == "--max-sector-dim") {
      p.max_sector_dim = parse_value<std::size_t>(next(arg));
    } else if (arg == "--checkpoint") {
      p.checkpoint = next(arg);
    } else if (arg == "--trace-partition-csv") {
      p.trace_partition_csv = next(arg);
    } else if (arg == "--debug-block-nup") {
      p.debug_block_nup = parse_value<int>(next(arg));
    } else if (arg == "--debug-block-ndown") {
      p.debug_block_ndown = parse_value<int>(next(arg));
    } else if (arg == "--debug-block-mx") {
      p.debug_block_mx = parse_value<int>(next(arg));
    } else if (arg == "--debug-block-my") {
      p.debug_block_my = parse_value<int>(next(arg));
    } else if (arg == "--debug-block-csv") {
      p.debug_block_csv = next(arg);
    } else if (arg == "--output") {
      p.output = next(arg);
    } else {
      die("Unknown argument: " + arg);
    }
  }

  if (p.lx <= 0 || p.ly <= 0) {
    die("Lattice dimensions must be positive.");
  }
  if (p.lx * p.ly > 63) {
    die("This implementation supports at most 63 sites.");
  }
  if (p.beta_values.empty()) {
    p.beta_values.push_back(p.beta);
  } else {
    p.beta = p.beta_values.front();
  }
  for (double beta : p.beta_values) {
    if (beta <= 0.0) {
      die("All beta values must be positive.");
    }
  }
  if (p.samples <= 0 || p.lanczos_steps <= 0) {
    die("--samples and --lanczos-steps must be positive.");
  }
  if (p.threads < 0) {
    die("--threads must be non-negative.");
  }
  if (p.mu_count < 2) {
    die("--mu-count must be at least 2.");
  }
  const bool partial_block_debug =
      (p.debug_block_nup >= 0) || (p.debug_block_ndown >= 0) ||
      (p.debug_block_mx >= 0) || (p.debug_block_my >= 0) ||
      !p.debug_block_csv.empty();
  const bool full_block_debug =
      (p.debug_block_nup >= 0) && (p.debug_block_ndown >= 0) &&
      (p.debug_block_mx >= 0) && (p.debug_block_my >= 0);
  if (partial_block_debug && !full_block_debug) {
    die("Selected-block debug requires --debug-block-nup/--debug-block-ndown/--debug-block-mx/--debug-block-my together.");
  }
  return p;
}

std::string checkpoint_metadata_text(const Params& params) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "lx=" << params.lx << "\n"
      << "ly=" << params.ly << "\n"
      << "tx=" << params.tx << "\n"
      << "ty=" << params.ty << "\n"
      << "tp=" << params.tp << "\n"
      << "phix=" << params.phix << "\n"
      << "phiy=" << params.phiy << "\n"
      << "u=" << params.u << "\n"
      << "v=" << params.v << "\n"
      << "samples=" << params.samples << "\n"
      << "lanczos_steps=" << params.lanczos_steps << "\n"
      << "exact_block_threshold=" << params.exact_block_threshold << "\n"
      << "seed=" << params.seed << "\n";
  return out.str();
}

void validate_or_write_checkpoint_metadata(const Params& params) {
  if (params.checkpoint.empty()) {
    return;
  }
  const std::string path = params.checkpoint + ".meta";
  const std::string expected = checkpoint_metadata_text(params);
  {
    std::ifstream in(path);
    if (in) {
      const std::string existing(
          (std::istreambuf_iterator<char>(in)),
          std::istreambuf_iterator<char>());
      if (existing != expected) {
        die("Checkpoint metadata does not match current run: " + path);
      }
      return;
    }
  }
  std::ofstream out(path);
  if (!out) {
    die("Failed to write checkpoint metadata: " + path);
  }
  out << expected;
}

int popcount64(std::uint64_t value) {
  return static_cast<int>(__builtin_popcountll(value));
}

std::vector<int> occupied_sites(std::uint64_t bits, int sites) {
  std::vector<int> occ;
  occ.reserve(static_cast<std::size_t>(popcount64(bits)));
  for (int i = 0; i < sites; ++i) {
    if ((bits >> i) & 1ULL) {
      occ.push_back(i);
    }
  }
  return occ;
}

Lattice build_lattice(int lx, int ly) {
  Lattice lattice;
  lattice.lx = lx;
  lattice.ly = ly;
  lattice.sites = lx * ly;
  lattice.neighbors.resize(static_cast<std::size_t>(lattice.sites));

  auto index = [lx](int x, int y) { return y * lx + x; };
  for (int y = 0; y < ly; ++y) {
    for (int x = 0; x < lx; ++x) {
      const int s = index(x, y);
      lattice.neighbors[s] = {
          index((x + 1) % lx, y),
          index((x - 1 + lx) % lx, y),
          index(x, (y + 1) % ly),
          index(x, (y - 1 + ly) % ly),
          index((x + 1) % lx, (y + 1) % ly),
          index((x - 1 + lx) % lx, (y + 1) % ly),
          index((x + 1) % lx, (y - 1 + ly) % ly),
          index((x - 1 + lx) % lx, (y - 1 + ly) % ly),
      };
      lattice.unique_bonds.push_back({s, index((x + 1) % lx, y)});
      lattice.unique_bonds.push_back({s, index(x, (y + 1) % ly)});
    }
  }

  std::sort(lattice.unique_bonds.begin(), lattice.unique_bonds.end());
  lattice.unique_bonds.erase(
      std::unique(lattice.unique_bonds.begin(), lattice.unique_bonds.end()),
      lattice.unique_bonds.end());
  return lattice;
}

std::vector<std::uint64_t> generate_combinations(int sites, int particles) {
  std::vector<std::uint64_t> states;
  if (particles < 0 || particles > sites) {
    return states;
  }
  if (particles == 0) {
    states.push_back(0);
    return states;
  }

  std::uint64_t state = (std::uint64_t{1} << particles) - 1ULL;
  const std::uint64_t limit = (std::uint64_t{1} << sites);
  while (state < limit) {
    states.push_back(state);
    const std::uint64_t c = state & (~state + 1ULL);
    const std::uint64_t r = state + c;
    if (r == 0) {
      break;
    }
    state = (((r ^ state) >> 2U) / c) | r;
  }
  return states;
}

double diagonal_energy(
    const StateKey& state,
    const Lattice& lattice,
    double u,
    double v) {
  const int doublons = popcount64(state.up & state.down);
  int nn_density_sum = 0;
  for (const auto& bond : lattice.unique_bonds) {
    const int ni = static_cast<int>((state.up >> bond.first) & 1ULL) +
                   static_cast<int>((state.down >> bond.first) & 1ULL);
    const int nj = static_cast<int>((state.up >> bond.second) & 1ULL) +
                   static_cast<int>((state.down >> bond.second) & 1ULL);
    nn_density_sum += ni * nj;
  }
  return u * static_cast<double>(doublons) + v * static_cast<double>(nn_density_sum);
}

double fermion_sign_hop(std::uint64_t bits, int from, int to) {
  const int lo = std::min(from, to);
  const int hi = std::max(from, to);
  if (hi - lo <= 1) {
    return 1.0;
  }
  const std::uint64_t width = static_cast<std::uint64_t>(hi - lo - 1);
  const std::uint64_t mask = ((std::uint64_t{1} << width) - 1ULL) << (lo + 1);
  return (popcount64(bits & mask) & 1) ? -1.0 : 1.0;
}

std::uint64_t mix_u64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

std::uint64_t make_block_seed(
    const Params& params,
    int n_up,
    int n_down,
    int mx,
    int my) {
  std::uint64_t seed = mix_u64(params.seed);
  seed ^= mix_u64(static_cast<std::uint64_t>(n_up + 1));
  seed ^= mix_u64(static_cast<std::uint64_t>(n_down + 17));
  seed ^= mix_u64(static_cast<std::uint64_t>(mx + 257));
  seed ^= mix_u64(static_cast<std::uint64_t>(my + 65537));
  return mix_u64(seed);
}

int requested_omp_threads(const Params& params) {
#if defined(_OPENMP)
  return (params.threads > 0) ? params.threads : omp_get_max_threads();
#else
  (void)params;
  return 1;
#endif
}

Complex twist_phase_for_dir(const Params& params, const Lattice& lattice, int dir) {
  const double theta_x = 2.0 * kPi * params.phix / static_cast<double>(lattice.lx);
  const double theta_y = 2.0 * kPi * params.phiy / static_cast<double>(lattice.ly);
  double angle = 0.0;
  switch (dir) {
    case 0:
      angle = theta_x;
      break;
    case 1:
      angle = -theta_x;
      break;
    case 2:
      angle = theta_y;
      break;
    case 3:
      angle = -theta_y;
      break;
    case 4:
      angle = theta_x + theta_y;
      break;
    case 5:
      angle = -theta_x + theta_y;
      break;
    case 6:
      angle = theta_x - theta_y;
      break;
    case 7:
      angle = -theta_x - theta_y;
      break;
    default:
      die("Invalid hopping direction for twist phase.");
  }
  return Complex(std::cos(angle), std::sin(angle));
}

std::pair<std::uint64_t, int> translate_spin(
    std::uint64_t bits,
    int dx,
    int dy,
    int lx,
    int ly) {
  std::vector<int> moved;
  moved.reserve(static_cast<std::size_t>(popcount64(bits)));
  const int sites = lx * ly;
  for (int s = 0; s < sites; ++s) {
    if (!((bits >> s) & 1ULL)) {
      continue;
    }
    const int x = s % lx;
    const int y = s / lx;
    const int nx = (x + dx) % lx;
    const int ny = (y + dy) % ly;
    moved.push_back(ny * lx + nx);
  }

  int inversions = 0;
  for (std::size_t i = 0; i < moved.size(); ++i) {
    for (std::size_t j = i + 1; j < moved.size(); ++j) {
      if (moved[i] > moved[j]) {
        ++inversions;
      }
    }
  }

  std::uint64_t translated = 0;
  for (int site : moved) {
    translated |= (std::uint64_t{1} << site);
  }
  return {translated, (inversions & 1) ? -1 : 1};
}

StateKey translate_state(const StateKey& state, int dx, int dy, int lx, int ly, int& sign) {
  const auto up = translate_spin(state.up, dx, dy, lx, ly);
  const auto down = translate_spin(state.down, dx, dy, lx, ly);
  sign = up.second * down.second;
  return {up.first, down.first};
}

SectorBasis build_sector_basis(const Lattice& lattice, int n_up, int n_down, const Params& params) {
  SectorBasis basis;
  basis.n_up = n_up;
  basis.n_down = n_down;
  basis.particles = n_up + n_down;

  const auto up_states = generate_combinations(lattice.sites, n_up);
  const auto down_states = generate_combinations(lattice.sites, n_down);
  basis.full_dim = up_states.size() * down_states.size();

  std::vector<std::pair<int, int>> translations;
  translations.reserve(static_cast<std::size_t>(lattice.sites));
  for (int dy = 0; dy < lattice.ly; ++dy) {
    for (int dx = 0; dx < lattice.lx; ++dx) {
      translations.push_back({dx, dy});
    }
  }

  std::unordered_map<StateKey, Relation, StateKeyHash> relations;
  relations.reserve(basis.full_dim);

  for (const auto up : up_states) {
    for (const auto down : down_states) {
      const StateKey state{up, down};
      if (relations.find(state) != relations.end()) {
        continue;
      }

      const int parent_index = static_cast<int>(basis.parents.size());
      ParentData parent;
      parent.representative = state;
      parent.diagonal = diagonal_energy(state, lattice, params.u, params.v);

      relations.emplace(state, Relation{parent_index, 0, 1});
      int degeneracy = 0;
      for (int shift = 0; shift < static_cast<int>(translations.size()); ++shift) {
        int sign = 1;
        const StateKey translated = translate_state(
            state,
            translations[static_cast<std::size_t>(shift)].first,
            translations[static_cast<std::size_t>(shift)].second,
            lattice.lx,
            lattice.ly,
            sign);
        if (translated == state) {
          ++degeneracy;
        } else {
          relations.emplace(translated, Relation{parent_index, shift, sign});
        }
      }
      parent.degeneracy = degeneracy;
      basis.parents.push_back(std::move(parent));
    }
  }

  for (ParentData& parent : basis.parents) {
    const StateKey state = parent.representative;
    const auto up_occ = occupied_sites(state.up, lattice.sites);
    for (int from : up_occ) {
      for (int dir = 0; dir < 8; ++dir) {
        const int to = lattice.neighbors[from][dir];
        if ((state.up >> to) & 1ULL) {
          continue;
        }
        const double hop = (dir < 2) ? params.tx : ((dir < 4) ? params.ty : params.tp);
        const Complex twisted_hop =
            -hop * fermion_sign_hop(state.up, from, to) *
            twist_phase_for_dir(params, lattice, dir);
        const StateKey target{
            state.up ^ (std::uint64_t{1} << from) ^ (std::uint64_t{1} << to),
            state.down,
        };
        const Relation rel = relations.at(target);
        parent.hops.push_back({rel.parent, rel.shift, rel.sign, twisted_hop});
      }
    }

    const auto down_occ = occupied_sites(state.down, lattice.sites);
    for (int from : down_occ) {
      for (int dir = 0; dir < 8; ++dir) {
        const int to = lattice.neighbors[from][dir];
        if ((state.down >> to) & 1ULL) {
          continue;
        }
        const double hop = (dir < 2) ? params.tx : ((dir < 4) ? params.ty : params.tp);
        const Complex twisted_hop =
            -hop * fermion_sign_hop(state.down, from, to) *
            twist_phase_for_dir(params, lattice, dir);
        const StateKey target{
            state.up,
            state.down ^ (std::uint64_t{1} << from) ^ (std::uint64_t{1} << to),
        };
        const Relation rel = relations.at(target);
        parent.hops.push_back({rel.parent, rel.shift, rel.sign, twisted_hop});
      }
    }
  }

  return basis;
}

MomentumBlock build_momentum_block(
    const SectorBasis& sector,
    const Lattice& lattice,
    int mx,
    int my) {
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
    block.diagonal[static_cast<std::size_t>(src)] = sector.parents[p].diagonal;
    for (const CompactHop& hop : sector.parents[p].hops) {
      const int tgt = active_index[static_cast<std::size_t>(hop.parent)];
      if (tgt < 0) {
        continue;
      }
      const Complex coeff =
          hop.amplitude * static_cast<double>(hop.sign) *
          phase[static_cast<std::size_t>(hop.shift)] *
          (dnf[static_cast<std::size_t>(hop.parent)] / dnf[p]);
      block.col_idx.push_back(tgt);
      block.values.push_back(coeff);
    }
    block.row_ptr.push_back(static_cast<int>(block.col_idx.size()));
  }

  return block;
}

double complex_norm2(const std::vector<Complex>& vec) {
  double sum = 0.0;
  for (const Complex& value : vec) {
    sum += std::norm(value);
  }
  return sum;
}

Complex dot(const std::vector<Complex>& a, const std::vector<Complex>& b) {
  Complex sum = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    sum += std::conj(a[i]) * b[i];
  }
  return sum;
}

void normalize(std::vector<Complex>& vec) {
  const double norm = std::sqrt(complex_norm2(vec));
  if (norm <= 0.0) {
    die("Attempted to normalize a zero vector.");
  }
  for (Complex& value : vec) {
    value /= norm;
  }
}

void apply_hamiltonian(
    const MomentumBlock& block,
    const std::vector<Complex>& in,
    std::vector<Complex>& out) {
  std::fill(out.begin(), out.end(), Complex{0.0, 0.0});
  for (std::size_t i = 0; i < block.basis_dim; ++i) {
    out[i] += block.diagonal[i] * in[i];
    for (int edge = block.row_ptr[i]; edge < block.row_ptr[i + 1]; ++edge) {
      out[i] += block.values[static_cast<std::size_t>(edge)] *
                in[static_cast<std::size_t>(block.col_idx[static_cast<std::size_t>(edge)])];
    }
  }
}

std::vector<double> diagonalize_block_exact(const MomentumBlock& block) {
  int n = static_cast<int>(block.basis_dim);
  if (n == 0) {
    return {};
  }
  if (n == 1) {
    double eigenvalue = block.diagonal[0];
    for (int edge = block.row_ptr[0]; edge < block.row_ptr[1]; ++edge) {
      if (block.col_idx[static_cast<std::size_t>(edge)] == 0) {
        eigenvalue += block.values[static_cast<std::size_t>(edge)].real();
      }
    }
    return {eigenvalue};
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

void jacobi_diagonalize(
    std::vector<std::vector<double>>& matrix,
    std::vector<double>& eigenvalues,
    std::vector<std::vector<double>>& eigenvectors) {
  const int n = static_cast<int>(matrix.size());
  if (n == 0) {
    eigenvalues.clear();
    eigenvectors.clear();
    return;
  }
  if (n == 1) {
    eigenvalues = {matrix[0][0]};
    eigenvectors = {{1.0}};
    return;
  }

  eigenvectors.assign(static_cast<std::size_t>(n), std::vector<double>(static_cast<std::size_t>(n), 0.0));
  for (int i = 0; i < n; ++i) {
    eigenvectors[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 1.0;
  }

  for (int iter = 0; iter < 50 * n * n; ++iter) {
    int p = 0;
    int q = 1;
    double max_offdiag = 0.0;
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        const double value = std::abs(matrix[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
        if (value > max_offdiag) {
          max_offdiag = value;
          p = i;
          q = j;
        }
      }
    }
    if (max_offdiag < 1e-12) {
      break;
    }

    const double app = matrix[static_cast<std::size_t>(p)][static_cast<std::size_t>(p)];
    const double aqq = matrix[static_cast<std::size_t>(q)][static_cast<std::size_t>(q)];
    const double apq = matrix[static_cast<std::size_t>(p)][static_cast<std::size_t>(q)];
    const double tau = (aqq - app) / (2.0 * apq);
    const double t =
        (tau >= 0.0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
    const double c = 1.0 / std::sqrt(1.0 + t * t);
    const double s = t * c;

    for (int k = 0; k < n; ++k) {
      if (k == p || k == q) {
        continue;
      }
      const double mkp = matrix[static_cast<std::size_t>(k)][static_cast<std::size_t>(p)];
      const double mkq = matrix[static_cast<std::size_t>(k)][static_cast<std::size_t>(q)];
      matrix[static_cast<std::size_t>(k)][static_cast<std::size_t>(p)] = c * mkp - s * mkq;
      matrix[static_cast<std::size_t>(p)][static_cast<std::size_t>(k)] =
          matrix[static_cast<std::size_t>(k)][static_cast<std::size_t>(p)];
      matrix[static_cast<std::size_t>(k)][static_cast<std::size_t>(q)] = s * mkp + c * mkq;
      matrix[static_cast<std::size_t>(q)][static_cast<std::size_t>(k)] =
          matrix[static_cast<std::size_t>(k)][static_cast<std::size_t>(q)];
    }

    matrix[static_cast<std::size_t>(p)][static_cast<std::size_t>(p)] =
        c * c * app - 2.0 * s * c * apq + s * s * aqq;
    matrix[static_cast<std::size_t>(q)][static_cast<std::size_t>(q)] =
        s * s * app + 2.0 * s * c * apq + c * c * aqq;
    matrix[static_cast<std::size_t>(p)][static_cast<std::size_t>(q)] = 0.0;
    matrix[static_cast<std::size_t>(q)][static_cast<std::size_t>(p)] = 0.0;

    for (int k = 0; k < n; ++k) {
      const double vkp = eigenvectors[static_cast<std::size_t>(k)][static_cast<std::size_t>(p)];
      const double vkq = eigenvectors[static_cast<std::size_t>(k)][static_cast<std::size_t>(q)];
      eigenvectors[static_cast<std::size_t>(k)][static_cast<std::size_t>(p)] = c * vkp - s * vkq;
      eigenvectors[static_cast<std::size_t>(k)][static_cast<std::size_t>(q)] = s * vkp + c * vkq;
    }
  }

  eigenvalues.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    eigenvalues[static_cast<std::size_t>(i)] = matrix[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
  }
}

LanczosSpectrum run_ftlm_block(
    const MomentumBlock& block,
    const Params& params,
    std::uint64_t block_seed) {
  LanczosSpectrum spectrum;
  spectrum.basis_dim = block.basis_dim;
  spectrum.particles = block.particles;

  const int krylov_dim =
      std::min<int>(params.lanczos_steps, static_cast<int>(block.basis_dim));
  std::vector<std::vector<double>> eigenvalues_by_sample(
      static_cast<std::size_t>(params.samples));
  std::vector<std::vector<double>> sample_weights_by_sample(
      static_cast<std::size_t>(params.samples));

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) num_threads(params.threads > 0 ? params.threads : omp_get_max_threads())
#endif
  for (int sample = 0; sample < params.samples; ++sample) {
    std::mt19937_64 rng(mix_u64(block_seed ^ static_cast<std::uint64_t>(sample)));
    std::normal_distribution<double> normal(0.0, 1.0);
    std::vector<Complex> q_prev(block.basis_dim, 0.0);
    std::vector<Complex> q_cur(block.basis_dim, 0.0);
    std::vector<Complex> q_next(block.basis_dim, 0.0);
    std::vector<Complex> hq(block.basis_dim, 0.0);

    for (Complex& value : q_cur) {
      value = Complex(normal(rng), normal(rng));
    }
    normalize(q_cur);
    std::fill(q_prev.begin(), q_prev.end(), Complex{0.0, 0.0});

    std::vector<double> alpha;
    std::vector<double> beta;
    alpha.reserve(static_cast<std::size_t>(krylov_dim));
    beta.reserve(static_cast<std::size_t>(krylov_dim));

    int actual_steps = 0;
    double beta_prev = 0.0;
    for (int step = 0; step < krylov_dim; ++step) {
      apply_hamiltonian(block, q_cur, hq);
      const Complex a_complex = dot(q_cur, hq);
      const double a = a_complex.real();
      for (std::size_t i = 0; i < block.basis_dim; ++i) {
        q_next[i] = hq[i] - a * q_cur[i] - beta_prev * q_prev[i];
      }

      alpha.push_back(a);
      const double b = std::sqrt(complex_norm2(q_next));
      ++actual_steps;
      if (b < 1e-12 || step + 1 == krylov_dim) {
        break;
      }

      beta.push_back(b);
      q_prev.swap(q_cur);
      q_cur.swap(q_next);
      for (Complex& value : q_cur) {
        value /= b;
      }
      beta_prev = b;
    }

    std::vector<std::vector<double>> tri(
        static_cast<std::size_t>(actual_steps),
        std::vector<double>(static_cast<std::size_t>(actual_steps), 0.0));
    for (int i = 0; i < actual_steps; ++i) {
      tri[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = alpha[static_cast<std::size_t>(i)];
      if (i + 1 < actual_steps) {
        tri[static_cast<std::size_t>(i)][static_cast<std::size_t>(i + 1)] =
            beta[static_cast<std::size_t>(i)];
        tri[static_cast<std::size_t>(i + 1)][static_cast<std::size_t>(i)] =
            beta[static_cast<std::size_t>(i)];
      }
    }

    std::vector<double> eigenvalues;
    std::vector<std::vector<double>> eigenvectors;
    jacobi_diagonalize(tri, eigenvalues, eigenvectors);

    std::vector<int> order(actual_steps);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
      return eigenvalues[static_cast<std::size_t>(lhs)] <
             eigenvalues[static_cast<std::size_t>(rhs)];
    });

    std::vector<double>& sample_eigenvalues =
        eigenvalues_by_sample[static_cast<std::size_t>(sample)];
    std::vector<double>& sample_weights =
        sample_weights_by_sample[static_cast<std::size_t>(sample)];
    sample_eigenvalues.reserve(static_cast<std::size_t>(actual_steps));
    sample_weights.reserve(static_cast<std::size_t>(actual_steps));
    for (int sorted = 0; sorted < actual_steps; ++sorted) {
      const int original = order[static_cast<std::size_t>(sorted)];
      const double v1 = eigenvectors[0][static_cast<std::size_t>(original)];
      sample_eigenvalues.push_back(eigenvalues[static_cast<std::size_t>(original)]);
      sample_weights.push_back((v1 * v1) / static_cast<double>(params.samples));
    }
  }

  for (int sample = 0; sample < params.samples; ++sample) {
    const std::vector<double>& sample_eigenvalues =
        eigenvalues_by_sample[static_cast<std::size_t>(sample)];
    const std::vector<double>& sample_weights =
        sample_weights_by_sample[static_cast<std::size_t>(sample)];
    spectrum.eigenvalues.insert(
        spectrum.eigenvalues.end(),
        sample_eigenvalues.begin(),
        sample_eigenvalues.end());
    spectrum.sample_weights.insert(
        spectrum.sample_weights.end(),
        sample_weights.begin(),
        sample_weights.end());
  }

  return spectrum;
}

BlockThermo raw_block_thermo_from_ftlm(
    const LanczosSpectrum& spectrum,
    double beta,
    double mu) {
  BlockThermo thermo;
  thermo.min_shifted_energy = std::numeric_limits<double>::infinity();
  for (double eigenvalue : spectrum.eigenvalues) {
    thermo.min_shifted_energy = std::min(
        thermo.min_shifted_energy,
        eigenvalue - mu * static_cast<double>(spectrum.particles));
  }
  if (!std::isfinite(thermo.min_shifted_energy)) {
    return thermo;
  }

  for (std::size_t i = 0; i < spectrum.eigenvalues.size(); ++i) {
    const double shifted =
        spectrum.eigenvalues[i] - mu * static_cast<double>(spectrum.particles) -
        thermo.min_shifted_energy;
    const double weight =
        static_cast<double>(spectrum.basis_dim) * spectrum.sample_weights[i] *
        std::exp(-beta * shifted);
    const double particles = static_cast<double>(spectrum.particles);
    thermo.z += weight;
    thermo.n_total += particles * weight;
    thermo.n2_total += particles * particles * weight;
  }
  return thermo;
}

BlockThermo raw_block_thermo_from_exact(
    const std::vector<double>& eigenvalues,
    int particles,
    double beta,
    double mu) {
  BlockThermo thermo;
  for (double eigenvalue : eigenvalues) {
    thermo.min_shifted_energy = std::min(
        thermo.min_shifted_energy,
        eigenvalue - mu * static_cast<double>(particles));
  }
  if (!std::isfinite(thermo.min_shifted_energy)) {
    return thermo;
  }

  for (double eigenvalue : eigenvalues) {
    const double shifted =
        eigenvalue - mu * static_cast<double>(particles) - thermo.min_shifted_energy;
    const double weight = std::exp(-beta * shifted);
    const double particle_count = static_cast<double>(particles);
    thermo.z += weight;
    thermo.n_total += particle_count * weight;
    thermo.n2_total += particle_count * particle_count * weight;
  }
  return thermo;
}

void accumulate_block_thermo(
    const std::vector<double>& eigenvalues,
    const std::vector<double>* sample_weights,
    double trace_prefactor,
    int particles,
    double beta,
    const std::vector<double>& mu_values,
    std::vector<double>& emin_by_mu,
    std::vector<double>& z_scaled,
    std::vector<double>& n_scaled,
    std::vector<double>& n2_scaled) {
  for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
    const double mu = mu_values[imu];
    double block_min = std::numeric_limits<double>::infinity();
    for (double eigenvalue : eigenvalues) {
      block_min = std::min(
          block_min,
          eigenvalue - mu * static_cast<double>(particles));
    }
    if (!std::isfinite(block_min)) {
      continue;
    }

    if (!std::isfinite(emin_by_mu[imu])) {
      emin_by_mu[imu] = block_min;
    } else if (block_min < emin_by_mu[imu]) {
      const double scale =
          std::exp(-beta * (emin_by_mu[imu] - block_min));
      z_scaled[imu] *= scale;
      n_scaled[imu] *= scale;
      n2_scaled[imu] *= scale;
      emin_by_mu[imu] = block_min;
    }

    const double particle_count = static_cast<double>(particles);
    for (std::size_t i = 0; i < eigenvalues.size(); ++i) {
      const double shifted =
          eigenvalues[i] - mu * particle_count - emin_by_mu[imu];
      const double sample_weight =
          (sample_weights == nullptr) ? 1.0 : (*sample_weights)[i];
      const double weight =
          trace_prefactor * sample_weight * std::exp(-beta * shifted);
      z_scaled[imu] += weight;
      n_scaled[imu] += particle_count * weight;
      n2_scaled[imu] += particle_count * particle_count * weight;
    }
  }
}

bool block_matches_debug(
    const Params& params,
    int n_up,
    int n_down,
    int mx,
    int my) {
  return n_up == params.debug_block_nup &&
         n_down == params.debug_block_ndown &&
         mx == params.debug_block_mx &&
         my == params.debug_block_my;
}

void write_selected_block_debug(
    const std::string& path,
    int n_up,
    int n_down,
    int mx,
    int my,
    const LanczosSpectrum& spectrum,
    const std::vector<double>& exact_eigenvalues,
    double beta,
    const std::vector<double>& mu_values) {
  std::ofstream out(path);
  if (!out) {
    die("Failed to open selected block debug file: " + path);
  }
  out << "nup,ndown,mx,my,basis_dim,mu,z_ftlm,z_ed,z_ratio,"
      << "n_ftlm_total,n_ed_total,n_diff,n2_ftlm_total,n2_ed_total,n2_diff,"
      << "min_ftlm,min_ed,sum_sample_weights\n";
  out << std::setprecision(15);
  const double sum_sample_weights =
      std::accumulate(spectrum.sample_weights.begin(), spectrum.sample_weights.end(), 0.0);
  for (double mu : mu_values) {
    const BlockThermo ftlm = raw_block_thermo_from_ftlm(spectrum, beta, mu);
    const BlockThermo exact = raw_block_thermo_from_exact(exact_eigenvalues, spectrum.particles, beta, mu);
    const double z_ratio = (exact.z > 0.0) ? (ftlm.z / exact.z) : 0.0;
    out << n_up << "," << n_down << "," << mx << "," << my << ","
        << spectrum.basis_dim << "," << mu << ","
        << ftlm.z << "," << exact.z << "," << z_ratio << ","
        << ftlm.n_total << "," << exact.n_total << ","
        << (ftlm.n_total - exact.n_total) << ","
        << ftlm.n2_total << "," << exact.n2_total << ","
        << (ftlm.n2_total - exact.n2_total) << ","
        << ftlm.min_shifted_energy << "," << exact.min_shifted_energy << ","
        << sum_sample_weights << "\n";
  }
}

template <typename T>
void write_binary(std::ostream& out, const T& value) {
  out.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template <typename T>
bool read_binary(std::istream& in, T& value) {
  return static_cast<bool>(in.read(reinterpret_cast<char*>(&value), sizeof(T)));
}

void ensure_checkpoint_header(const std::string& path) {
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
    die("Failed to create checkpoint file: " + path);
  }
  const char magic[8] = {'F', 'T', 'L', 'M', 'C', 'P', '1', '\n'};
  out.write(magic, sizeof(magic));
}

std::unordered_map<BlockKey, StoredSpectrum, BlockKeyHash> read_checkpoint(
    const std::string& path) {
  std::unordered_map<BlockKey, StoredSpectrum, BlockKeyHash> spectra;
  if (path.empty()) {
    return spectra;
  }
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return spectra;
  }

  char magic[8] = {};
  if (!in.read(magic, sizeof(magic))) {
    return spectra;
  }
  const char expected_magic[8] = {'F', 'T', 'L', 'M', 'C', 'P', '1', '\n'};
  if (std::memcmp(magic, expected_magic, sizeof(magic)) != 0) {
    die("Checkpoint has an unrecognized format: " + path);
  }

  const char expected_record[8] = {'B', 'L', 'K', 'R', 'E', 'C', '1', '\n'};
  while (true) {
    char record_magic[8] = {};
    if (!in.read(record_magic, sizeof(record_magic))) {
      break;
    }
    if (std::memcmp(record_magic, expected_record, sizeof(record_magic)) != 0) {
      break;
    }

    BlockKey key;
    StoredSpectrum spectrum;
    std::uint8_t exact_flag = 0;
    std::uint64_t basis_dim = 0;
    std::uint64_t eigen_count = 0;
    std::uint64_t weight_count = 0;
    if (!read_binary(in, key.n_up) ||
        !read_binary(in, key.n_down) ||
        !read_binary(in, key.mx) ||
        !read_binary(in, key.my) ||
        !read_binary(in, spectrum.particles) ||
        !read_binary(in, basis_dim) ||
        !read_binary(in, spectrum.trace_prefactor) ||
        !read_binary(in, exact_flag) ||
        !read_binary(in, eigen_count) ||
        !read_binary(in, weight_count)) {
      break;
    }
    if (eigen_count > 100000000ULL || weight_count > 100000000ULL) {
      die("Checkpoint record is implausibly large: " + path);
    }
    spectrum.basis_dim = static_cast<std::size_t>(basis_dim);
    spectrum.exact = exact_flag != 0;
    spectrum.eigenvalues.resize(static_cast<std::size_t>(eigen_count));
    spectrum.sample_weights.resize(static_cast<std::size_t>(weight_count));
    if (!in.read(reinterpret_cast<char*>(spectrum.eigenvalues.data()),
                 static_cast<std::streamsize>(sizeof(double) * spectrum.eigenvalues.size())) ||
        !in.read(reinterpret_cast<char*>(spectrum.sample_weights.data()),
                 static_cast<std::streamsize>(sizeof(double) * spectrum.sample_weights.size()))) {
      break;
    }
    spectra[key] = std::move(spectrum);
  }
  return spectra;
}

void append_checkpoint_record(
    const std::string& path,
    const BlockKey& key,
    const StoredSpectrum& spectrum) {
  if (path.empty()) {
    return;
  }
  ensure_checkpoint_header(path);
  std::ofstream out(path, std::ios::binary | std::ios::app);
  if (!out) {
    die("Failed to append checkpoint file: " + path);
  }
  const char record_magic[8] = {'B', 'L', 'K', 'R', 'E', 'C', '1', '\n'};
  out.write(record_magic, sizeof(record_magic));
  write_binary(out, key.n_up);
  write_binary(out, key.n_down);
  write_binary(out, key.mx);
  write_binary(out, key.my);
  write_binary(out, spectrum.particles);
  const std::uint64_t basis_dim = static_cast<std::uint64_t>(spectrum.basis_dim);
  write_binary(out, basis_dim);
  write_binary(out, spectrum.trace_prefactor);
  const std::uint8_t exact_flag = spectrum.exact ? 1 : 0;
  write_binary(out, exact_flag);
  const std::uint64_t eigen_count =
      static_cast<std::uint64_t>(spectrum.eigenvalues.size());
  const std::uint64_t weight_count =
      static_cast<std::uint64_t>(spectrum.sample_weights.size());
  write_binary(out, eigen_count);
  write_binary(out, weight_count);
  out.write(reinterpret_cast<const char*>(spectrum.eigenvalues.data()),
            static_cast<std::streamsize>(sizeof(double) * spectrum.eigenvalues.size()));
  out.write(reinterpret_cast<const char*>(spectrum.sample_weights.data()),
            static_cast<std::streamsize>(sizeof(double) * spectrum.sample_weights.size()));
  out.flush();
}

std::vector<double> linspace(double start, double stop, int count) {
  std::vector<double> values(static_cast<std::size_t>(count));
  const double step = (stop - start) / static_cast<double>(count - 1);
  for (int i = 0; i < count; ++i) {
    values[static_cast<std::size_t>(i)] = start + step * static_cast<double>(i);
  }
  return values;
}

void write_results(
    const std::string& path,
    const std::vector<double>& mu_values,
    const std::vector<double>& densities,
    const std::vector<double>& charge_correlations,
    const std::vector<double>& compressibilities,
    const std::vector<double>& partitions) {
  std::ofstream out(path);
  if (!out) {
    die("Failed to open output file: " + path);
  }
  out << "mu,n,charge_correlation,compressibility,partition_like\n";
  out << std::setprecision(15);
  for (std::size_t i = 0; i < mu_values.size(); ++i) {
    out << mu_values[i] << "," << densities[i] << ","
        << charge_correlations[i] << "," << compressibilities[i] << ","
        << partitions[i] << "\n";
  }
}

void write_multi_beta_results(
    const std::string& path,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    const std::vector<std::vector<double>>& densities,
    const std::vector<std::vector<double>>& charge_correlations,
    const std::vector<std::vector<double>>& compressibilities,
    const std::vector<std::vector<double>>& partitions) {
  std::ofstream out(path);
  if (!out) {
    die("Failed to open output file: " + path);
  }
  out << "beta,mu,n,charge_correlation,compressibility,partition_like\n";
  out << std::setprecision(15);
  for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
    for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
      out << beta_values[ibeta] << "," << mu_values[imu] << ","
          << densities[ibeta][imu] << ","
          << charge_correlations[ibeta][imu] << ","
          << compressibilities[ibeta][imu] << ","
          << partitions[ibeta][imu] << "\n";
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Params params = parse_args(argc, argv);
    const Lattice lattice = build_lattice(params.lx, params.ly);
    const bool multi_beta = params.beta_values.size() > 1;
    const std::vector<double> mu_values =
        linspace(params.mu_min, params.mu_max, params.mu_count);
    const bool trace_blocks = !params.trace_partition_csv.empty();
    const bool debug_selected_block =
        params.debug_block_nup >= 0 && params.debug_block_ndown >= 0 &&
        params.debug_block_mx >= 0 && params.debug_block_my >= 0;
    if (!params.checkpoint.empty() && (trace_blocks || debug_selected_block)) {
      die("--checkpoint is not currently compatible with trace/debug block output.");
    }
    validate_or_write_checkpoint_metadata(params);
    ensure_checkpoint_header(params.checkpoint);
    auto checkpoint_spectra = read_checkpoint(params.checkpoint);

    std::cout << "Reference used: "
              << "reference_ftlm_hub_cond/fltm_hub_cond/hubTri2Dcond_omp.f"
              << " and cond_spect_omp.f\n";
    std::cout << "Running FTLM with translation-symmetry k blocks on a "
              << params.lx << "x" << params.ly << " rectangular lattice ("
              << lattice.sites << " sites)\n";
    std::cout << "beta=";
    for (std::size_t i = 0; i < params.beta_values.size(); ++i) {
      if (i > 0) {
        std::cout << ",";
      }
      std::cout << params.beta_values[i];
    }
    std::cout << " U=" << params.u << " V=" << params.v
              << " tx=" << params.tx << " ty=" << params.ty << " tp=" << params.tp
              << " phix=" << params.phix << " phiy=" << params.phiy
              << " threads=" << requested_omp_threads(params)
              << " exact_block_threshold=" << params.exact_block_threshold << "\n";
    if (!params.checkpoint.empty()) {
      std::cout << "checkpoint=" << params.checkpoint
                << " loaded_blocks=" << checkpoint_spectra.size() << "\n";
    }
    std::vector<StoredSpectrum> stored_spectra;
    std::vector<double> single_emin_by_mu(
        mu_values.size(), std::numeric_limits<double>::infinity());
    std::vector<double> single_z_scaled(mu_values.size(), 0.0);
    std::vector<double> single_n_scaled(mu_values.size(), 0.0);
    std::vector<double> single_n2_scaled(mu_values.size(), 0.0);
    std::ofstream trace_out;
    if (trace_blocks) {
      trace_out.open(params.trace_partition_csv);
      if (!trace_out) {
        die("Failed to open trace CSV: " + params.trace_partition_csv);
      }
      trace_out << "nup,ndown,mx,my,particles,basis_dim,mu,"
                << "method,z_used,z_ed,z_ratio,n_used_total,n_ed_total,n_total_diff,"
                << "n2_used_total,n2_ed_total,n2_total_diff,"
                << "min_ftlm,min_ed,sum_sample_weights\n";
      trace_out << std::setprecision(15);
    }
    bool wrote_selected_block_debug = false;

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

            const BlockKey block_key{n_up, n_down, mx, my};
            const auto checkpoint_it = checkpoint_spectra.find(block_key);
            if (checkpoint_it != checkpoint_spectra.end()) {
              const StoredSpectrum& stored = checkpoint_it->second;
              if (stored.basis_dim != block.basis_dim ||
                  stored.particles != block.particles) {
                die("Checkpoint block metadata does not match current run.");
              }
              if (multi_beta) {
                stored_spectra.push_back(stored);
              } else {
                const std::vector<double>* sample_weights =
                    stored.sample_weights.empty() ? nullptr : &stored.sample_weights;
                accumulate_block_thermo(
                    stored.eigenvalues,
                    sample_weights,
                    stored.trace_prefactor,
                    stored.particles,
                    params.beta,
                    mu_values,
                    single_emin_by_mu,
                    single_z_scaled,
                    single_n_scaled,
                    single_n2_scaled);
              }
              continue;
            }

            const bool use_exact_block = block.basis_dim <= params.exact_block_threshold;
            const bool need_exact_block =
                use_exact_block || trace_blocks || (debug_selected_block &&
                                 block_matches_debug(params, n_up, n_down, mx, my));
            LanczosSpectrum spectrum;
            if (!use_exact_block) {
              spectrum =
                  run_ftlm_block(block, params, make_block_seed(params, n_up, n_down, mx, my));
            }
            std::vector<double> exact_eigenvalues;
            if (need_exact_block) {
              exact_eigenvalues = diagonalize_block_exact(block);
            }

            if (trace_blocks) {
              const double sum_sample_weights = use_exact_block
                  ? 1.0
                  : std::accumulate(
                        spectrum.sample_weights.begin(),
                        spectrum.sample_weights.end(),
                        0.0);
              for (double mu : mu_values) {
                const BlockThermo used = use_exact_block
                    ? raw_block_thermo_from_exact(exact_eigenvalues, block.particles, params.beta, mu)
                    : raw_block_thermo_from_ftlm(spectrum, params.beta, mu);
                const BlockThermo exact =
                    raw_block_thermo_from_exact(exact_eigenvalues, block.particles, params.beta, mu);
                const double z_ratio = (exact.z > 0.0) ? (used.z / exact.z) : 0.0;
                trace_out << n_up << "," << n_down << "," << mx << "," << my << ","
                          << block.particles << "," << block.basis_dim << ","
                          << mu << "," << (use_exact_block ? "exact" : "ftlm") << ","
                          << used.z << "," << exact.z << "," << z_ratio
                          << "," << used.n_total << "," << exact.n_total << ","
                          << (used.n_total - exact.n_total) << ","
                          << used.n2_total << "," << exact.n2_total << ","
                          << (used.n2_total - exact.n2_total) << ","
                          << used.min_shifted_energy << "," << exact.min_shifted_energy << ","
                          << sum_sample_weights << "\n";
              }
            }

            if (debug_selected_block && block_matches_debug(params, n_up, n_down, mx, my)) {
              const double sum_sample_weights = use_exact_block
                  ? 1.0
                  : std::accumulate(
                        spectrum.sample_weights.begin(),
                        spectrum.sample_weights.end(),
                        0.0);
              std::cout << "Selected block debug Nup=" << n_up
                        << " Ndown=" << n_down
                        << " mx=" << mx
                        << " my=" << my
                        << " basis_dim=" << block.basis_dim
                        << " method=" << (use_exact_block ? "exact" : "ftlm")
                        << " sum_sample_weights=" << sum_sample_weights << "\n";
              for (double mu : mu_values) {
                const BlockThermo used = use_exact_block
                    ? raw_block_thermo_from_exact(exact_eigenvalues, block.particles, params.beta, mu)
                    : raw_block_thermo_from_ftlm(spectrum, params.beta, mu);
                const BlockThermo exact =
                    raw_block_thermo_from_exact(exact_eigenvalues, block.particles, params.beta, mu);
                const double z_ratio = (exact.z > 0.0) ? (used.z / exact.z) : 0.0;
                std::cout << "  debug mu=" << std::setw(12) << mu
                          << " z_used=" << std::setw(14) << used.z
                          << " z_ed=" << std::setw(14) << exact.z
                          << " ratio=" << std::setw(14) << z_ratio
                          << " n_diff=" << std::setw(14) << (used.n_total - exact.n_total)
                          << "\n";
              }
              if (!params.debug_block_csv.empty() && !use_exact_block) {
                write_selected_block_debug(
                    params.debug_block_csv,
                    n_up,
                    n_down,
                    mx,
                    my,
                    spectrum,
                    exact_eigenvalues,
                    params.beta,
                    mu_values);
                wrote_selected_block_debug = true;
              } else if (!params.debug_block_csv.empty()) {
                write_selected_block_debug(
                    params.debug_block_csv,
                    n_up,
                    n_down,
                    mx,
                    my,
                    LanczosSpectrum{{}, {}, block.basis_dim, block.particles},
                    exact_eigenvalues,
                    params.beta,
                    mu_values);
                wrote_selected_block_debug = true;
              }
            }

            if (multi_beta || !params.checkpoint.empty()) {
              StoredSpectrum stored;
              stored.particles = block.particles;
              stored.basis_dim = block.basis_dim;
              stored.exact = use_exact_block;
              if (use_exact_block) {
                stored.eigenvalues = std::move(exact_eigenvalues);
                stored.trace_prefactor = 1.0;
              } else {
                stored.eigenvalues = std::move(spectrum.eigenvalues);
                stored.sample_weights = std::move(spectrum.sample_weights);
                stored.trace_prefactor = static_cast<double>(stored.basis_dim);
              }
              append_checkpoint_record(params.checkpoint, block_key, stored);
              if (multi_beta) {
                stored_spectra.push_back(std::move(stored));
              } else {
                const std::vector<double>* sample_weights =
                    stored.sample_weights.empty() ? nullptr : &stored.sample_weights;
                accumulate_block_thermo(
                    stored.eigenvalues,
                    sample_weights,
                    stored.trace_prefactor,
                    stored.particles,
                    params.beta,
                    mu_values,
                    single_emin_by_mu,
                    single_z_scaled,
                    single_n_scaled,
                    single_n2_scaled);
              }
            } else if (use_exact_block) {
              accumulate_block_thermo(
                  exact_eigenvalues,
                  nullptr,
                  1.0,
                  block.particles,
                  params.beta,
                  mu_values,
                  single_emin_by_mu,
                  single_z_scaled,
                  single_n_scaled,
                  single_n2_scaled);
            } else {
              accumulate_block_thermo(
                  spectrum.eigenvalues,
                  &spectrum.sample_weights,
                  static_cast<double>(spectrum.basis_dim),
                  spectrum.particles,
                  params.beta,
                  mu_values,
                  single_emin_by_mu,
                  single_z_scaled,
                  single_n_scaled,
                  single_n2_scaled);
            }
          }
        }
        std::cout << "  sector Nup=" << n_up << " Ndown=" << n_down
                  << " full_dim=" << sector.full_dim
                  << " parents=" << sector.parents.size()
                  << " active_k_dim_sum=" << active_sum << "\n";
      }
    }

    if (debug_selected_block && !params.debug_block_csv.empty() && !wrote_selected_block_debug) {
      die("Selected block debug target was not found among active FTLM blocks.");
    }

    if (multi_beta) {
      std::vector<std::vector<double>> all_densities(
          params.beta_values.size(), std::vector<double>(mu_values.size(), 0.0));
      std::vector<std::vector<double>> all_charge_correlations(
          params.beta_values.size(), std::vector<double>(mu_values.size(), 0.0));
      std::vector<std::vector<double>> all_compressibilities(
          params.beta_values.size(), std::vector<double>(mu_values.size(), 0.0));
      std::vector<std::vector<double>> all_partitions(
          params.beta_values.size(), std::vector<double>(mu_values.size(), 0.0));

      for (std::size_t ibeta = 0; ibeta < params.beta_values.size(); ++ibeta) {
        const double beta = params.beta_values[ibeta];
        std::vector<double> emin_by_mu(
            mu_values.size(), std::numeric_limits<double>::infinity());
        std::vector<double> z_scaled(mu_values.size(), 0.0);
        std::vector<double> n_scaled(mu_values.size(), 0.0);
        std::vector<double> n2_scaled(mu_values.size(), 0.0);

        for (const StoredSpectrum& stored : stored_spectra) {
          const std::vector<double>* sample_weights =
              stored.sample_weights.empty() ? nullptr : &stored.sample_weights;
          accumulate_block_thermo(
              stored.eigenvalues,
              sample_weights,
              stored.trace_prefactor,
              stored.particles,
              beta,
              mu_values,
              emin_by_mu,
              z_scaled,
              n_scaled,
              n2_scaled);
        }

        for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
          all_partitions[ibeta][imu] = z_scaled[imu];
          const double sites = static_cast<double>(lattice.sites);
          const double mean_particles = n_scaled[imu] / z_scaled[imu];
          const double mean_particles_squared = n2_scaled[imu] / z_scaled[imu];
          const double charge_correlation =
              (mean_particles_squared - mean_particles * mean_particles) / sites;
          all_densities[ibeta][imu] = mean_particles / sites;
          all_charge_correlations[ibeta][imu] = std::max(0.0, charge_correlation);
          all_compressibilities[ibeta][imu] =
              beta * all_charge_correlations[ibeta][imu];
          std::cout << "beta=" << std::setw(12) << beta
                    << "  mu=" << std::setw(12) << mu_values[imu]
                    << "  n=" << std::setw(12) << all_densities[ibeta][imu]
                    << "  kappa=" << std::setw(12)
                    << all_compressibilities[ibeta][imu] << "\n";
        }
      }
      write_multi_beta_results(
          params.output,
          params.beta_values,
          mu_values,
          all_densities,
          all_charge_correlations,
          all_compressibilities,
          all_partitions);
    } else {
      std::vector<double> densities(mu_values.size(), 0.0);
      std::vector<double> charge_correlations(mu_values.size(), 0.0);
      std::vector<double> compressibilities(mu_values.size(), 0.0);
      std::vector<double> partitions(mu_values.size(), 0.0);
      for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
        partitions[imu] = single_z_scaled[imu];
        const double sites = static_cast<double>(lattice.sites);
        const double mean_particles = single_n_scaled[imu] / single_z_scaled[imu];
        const double mean_particles_squared =
            single_n2_scaled[imu] / single_z_scaled[imu];
        const double charge_correlation =
            (mean_particles_squared - mean_particles * mean_particles) / sites;
        densities[imu] = mean_particles / sites;
        charge_correlations[imu] = std::max(0.0, charge_correlation);
        compressibilities[imu] = params.beta * charge_correlations[imu];
        std::cout << "beta=" << std::setw(12) << params.beta
                  << "  mu=" << std::setw(12) << mu_values[imu]
                  << "  n=" << std::setw(12) << densities[imu]
                  << "  kappa=" << std::setw(12)
                  << compressibilities[imu] << "\n";
      }
      write_results(
          params.output,
          mu_values,
          densities,
          charge_correlations,
          compressibilities,
          partitions);
    }
    std::cout << "Wrote " << params.output << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }
}
