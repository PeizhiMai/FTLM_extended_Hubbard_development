#include "thermo_checkpoint.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <csignal>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <ctime>
#include <unordered_map>
#include <utility>
#include <vector>

#include <sys/resource.h>
#include <unistd.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace {

using Complex = std::complex<double>;
using BlockKey = ftlm_checkpoint::BlockKey;
using BlockKeyHash = ftlm_checkpoint::BlockKeyHash;
using SectorKey = ftlm_checkpoint::SectorKey;

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
  double max_runtime_minutes = 0.0;
  std::string progress_jsonl;
  std::string progress_state;
  int progress_interval_seconds = 60;
  bool only_block = false;
  int only_block_nup = -1;
  int only_block_ndown = -1;
  int only_block_mx = -1;
  int only_block_my = -1;
  std::string twist_id;
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
  std::int8_t sign = 1;
  std::uint8_t dir = 0;
};
static_assert(sizeof(CompactHop) <= 12, "CompactHop must remain memory-compact for 4x4 sectors.");

struct ParentData {
  StateKey representative;
  int degeneracy = 0;
  double diagonal = 0.0;
};

struct SectorBasis {
  int n_up = 0;
  int n_down = 0;
  int particles = 0;
  std::size_t full_dim = 0;
  std::vector<ParentData> parents;
  std::vector<int> hop_row_ptr;
  std::vector<CompactHop> hops;
  std::array<Complex, 8> hop_amplitudes{};
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

std::array<int, 4> parse_block_key(const std::string& text) {
  std::array<int, 4> values{};
  std::istringstream in(text);
  std::string item;
  for (int i = 0; i < 4; ++i) {
    if (!std::getline(in, item, ',') || item.empty()) {
      die("--only-block expects NUP,NDOWN,MX,MY");
    }
    values[static_cast<std::size_t>(i)] = parse_value<int>(item);
  }
  if (std::getline(in, item, ',')) die("--only-block expects exactly four integers.");
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
      << "  --samples N           target total FTLM random vectors per k block\n"
      << "  --lanczos-steps N     Lanczos steps per random vector\n"
      << "  --threads N           OpenMP threads (0 uses runtime default)\n"
      << "  --exact-block-threshold N\n"
      << "                        use exact diagonalization for k blocks with dimension <= N\n"
      << "  --seed N              random seed\n"
      << "  --max-sector-dim N    abort if a full sector exceeds this size\n"
      << "  --checkpoint PATH     append/resume compact spectra (new files use extensible v2)\n"
      << "  --max-runtime-minutes X\n"
      << "                        stop cleanly after active samples finish\n"
      << "  --progress-jsonl PATH append machine-readable progress events\n"
      << "  --progress-state PATH atomically replace current progress JSON\n"
      << "  --progress-interval-seconds N\n"
      << "                        heartbeat interval (default 60)\n"
      << "  --only-block NUP,NDOWN,MX,MY\n"
      << "                        run only one block for a resource probe\n"
      << "  --twist-id TEXT       label included in progress output\n"
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
    } else if (arg == "--max-runtime-minutes") {
      p.max_runtime_minutes = parse_value<double>(next(arg));
    } else if (arg == "--progress-jsonl") {
      p.progress_jsonl = next(arg);
    } else if (arg == "--progress-state") {
      p.progress_state = next(arg);
    } else if (arg == "--progress-interval-seconds") {
      p.progress_interval_seconds = parse_value<int>(next(arg));
    } else if (arg == "--only-block") {
      const auto key = parse_block_key(next(arg));
      p.only_block = true;
      p.only_block_nup = key[0];
      p.only_block_ndown = key[1];
      p.only_block_mx = key[2];
      p.only_block_my = key[3];
    } else if (arg == "--twist-id") {
      p.twist_id = next(arg);
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
  if (p.max_runtime_minutes < 0.0) {
    die("--max-runtime-minutes must be non-negative.");
  }
  if (p.progress_interval_seconds <= 0) {
    die("--progress-interval-seconds must be positive.");
  }
  if (p.only_block &&
      (p.only_block_nup < 0 || p.only_block_nup > p.lx * p.ly ||
       p.only_block_ndown < 0 || p.only_block_ndown > p.lx * p.ly ||
       p.only_block_mx < 0 || p.only_block_mx >= p.lx ||
       p.only_block_my < 0 || p.only_block_my >= p.ly)) {
    die("--only-block is outside the requested lattice sectors/momenta.");
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

std::string checkpoint_metadata_text_v1(const Params& params) {
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

std::string checkpoint_metadata_text_v2(const Params& params) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "format=2\n"
      << "rng_version=" << ftlm_checkpoint::kRngVersion << "\n"
      << "rng_algorithm=" << ftlm_checkpoint::kRngAlgorithm << "\n"
      << "seed_derivation=" << ftlm_checkpoint::kSeedDerivation << "\n"
      << "lx=" << params.lx << "\n"
      << "ly=" << params.ly << "\n"
      << "tx=" << params.tx << "\n"
      << "ty=" << params.ty << "\n"
      << "tp=" << params.tp << "\n"
      << "phix=" << params.phix << "\n"
      << "phiy=" << params.phiy << "\n"
      << "u=" << params.u << "\n"
      << "v=" << params.v << "\n"
      << "lanczos_steps=" << params.lanczos_steps << "\n"
      << "exact_block_threshold=" << params.exact_block_threshold << "\n"
      << "seed=" << params.seed << "\n";
  return out.str();
}

ftlm_checkpoint::Format validate_or_write_checkpoint_metadata(const Params& params) {
  if (params.checkpoint.empty()) {
    return ftlm_checkpoint::Format::kMissing;
  }
  ftlm_checkpoint::Format format = ftlm_checkpoint::detect_format(params.checkpoint);
  if (format == ftlm_checkpoint::Format::kMissing ||
      format == ftlm_checkpoint::Format::kEmpty) {
    ftlm_checkpoint::initialize_v2_checkpoint(params.checkpoint);
    format = ftlm_checkpoint::Format::kV2;
  }
  if (format == ftlm_checkpoint::Format::kUnknown) {
    die("Checkpoint has an unrecognized format: " + params.checkpoint);
  }
  const std::string path = params.checkpoint + ".meta";
  const std::string expected = format == ftlm_checkpoint::Format::kV1
      ? checkpoint_metadata_text_v1(params)
      : checkpoint_metadata_text_v2(params);
  {
    std::ifstream in(path);
    if (in) {
      const std::string existing(
          (std::istreambuf_iterator<char>(in)),
          std::istreambuf_iterator<char>());
      if (existing != expected) {
        die("Checkpoint metadata does not match current run: " + path);
      }
      return format;
    }
  }
  std::ofstream out(path);
  if (!out) {
    die("Failed to write checkpoint metadata: " + path);
  }
  out << expected;
  return format;
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

// std::normal_distribution is intentionally avoided here: its transform is
// implementation-defined and can change when a compiler or standard library
// is upgraded.  This explicit 53-bit Box-Muller transform, together with the
// RNG/seed metadata above, makes every permanent sample ID reproducible.
class DeterministicNormal {
 public:
  explicit DeterministicNormal(std::uint64_t seed) : engine_(seed) {}

  double next() {
    if (has_spare_) {
      has_spare_ = false;
      return spare_;
    }
    const double u1 = uniform_open_01();
    const double u2 = uniform_open_01();
    const double radius = std::sqrt(-2.0 * std::log(u1));
    const double angle = 2.0 * kPi * u2;
    spare_ = radius * std::sin(angle);
    has_spare_ = true;
    return radius * std::cos(angle);
  }

 private:
  double uniform_open_01() {
    constexpr double inverse_two_to_53 = 1.0 / 9007199254740992.0;
    const std::uint64_t mantissa = engine_() >> 11U;
    return (static_cast<double>(mantissa) + 0.5) * inverse_two_to_53;
  }

  std::mt19937_64 engine_;
  bool has_spare_ = false;
  double spare_ = 0.0;
};

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

  int active_directions = 0;
  for (int dir = 0; dir < 8; ++dir) {
    const double hop = (dir < 2) ? params.tx : ((dir < 4) ? params.ty : params.tp);
    if (std::abs(hop) > kEps) {
      basis.hop_amplitudes[static_cast<std::size_t>(dir)] =
          -hop * twist_phase_for_dir(params, lattice, dir);
      ++active_directions;
    }
  }
  basis.hop_row_ptr.reserve(basis.parents.size() + 1);
  basis.hop_row_ptr.push_back(0);
  if (lattice.sites > 1 && !basis.parents.empty()) {
    const double expected_per_parent = static_cast<double>(active_directions) *
        (static_cast<double>(n_up * (lattice.sites - n_up)) +
         static_cast<double>(n_down * (lattice.sites - n_down))) /
        static_cast<double>(lattice.sites - 1);
    const long double reserve_count = static_cast<long double>(basis.parents.size()) *
        (1.15L * expected_per_parent + 2.0L);
    if (reserve_count < static_cast<long double>(std::numeric_limits<int>::max())) {
      basis.hops.reserve(static_cast<std::size_t>(reserve_count));
    }
  }

  for (const ParentData& parent : basis.parents) {
    const StateKey state = parent.representative;
    const auto up_occ = occupied_sites(state.up, lattice.sites);
    for (int from : up_occ) {
      for (int dir = 0; dir < 8; ++dir) {
        const int to = lattice.neighbors[from][dir];
        if ((state.up >> to) & 1ULL) {
          continue;
        }
        const double hop = (dir < 2) ? params.tx : ((dir < 4) ? params.ty : params.tp);
        if (std::abs(hop) <= kEps) {
          continue;
        }
        const StateKey target{
            state.up ^ (std::uint64_t{1} << from) ^ (std::uint64_t{1} << to),
            state.down,
        };
        const Relation rel = relations.at(target);
        const int sign = rel.sign *
            (fermion_sign_hop(state.up, from, to) > 0.0 ? 1 : -1);
        basis.hops.push_back({
            rel.parent, rel.shift, static_cast<std::int8_t>(sign),
            static_cast<std::uint8_t>(dir)});
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
        if (std::abs(hop) <= kEps) {
          continue;
        }
        const StateKey target{
            state.up,
            state.down ^ (std::uint64_t{1} << from) ^ (std::uint64_t{1} << to),
        };
        const Relation rel = relations.at(target);
        const int sign = rel.sign *
            (fermion_sign_hop(state.down, from, to) > 0.0 ? 1 : -1);
        basis.hops.push_back({
            rel.parent, rel.shift, static_cast<std::int8_t>(sign),
            static_cast<std::uint8_t>(dir)});
      }
    }
    if (basis.hops.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      die("Sector compact-hop table exceeds 32-bit indexing.");
    }
    basis.hop_row_ptr.push_back(static_cast<int>(basis.hops.size()));
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
    for (int edge = sector.hop_row_ptr[p]; edge < sector.hop_row_ptr[p + 1]; ++edge) {
      const CompactHop& hop = sector.hops[static_cast<std::size_t>(edge)];
      const int tgt = active_index[static_cast<std::size_t>(hop.parent)];
      if (tgt < 0) {
        continue;
      }
      const Complex coeff =
          sector.hop_amplitudes[static_cast<std::size_t>(hop.dir)] *
          static_cast<double>(hop.sign) *
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

struct LanczosSample {
  std::vector<double> eigenvalues;
  std::vector<double> overlaps;
};

LanczosSample run_ftlm_sample(
    const MomentumBlock& block,
    const Params& params,
    std::uint64_t block_seed,
    int sample_id,
    const std::function<void(int)>& step_callback = {}) {
  const int krylov_dim =
      std::min<int>(params.lanczos_steps, static_cast<int>(block.basis_dim));
  DeterministicNormal normal(
      mix_u64(block_seed ^ static_cast<std::uint64_t>(sample_id)));
  std::vector<Complex> q_prev(block.basis_dim, 0.0);
  std::vector<Complex> q_cur(block.basis_dim, 0.0);
  std::vector<Complex> q_next(block.basis_dim, 0.0);

  for (Complex& value : q_cur) value = Complex(normal.next(), normal.next());
  normalize(q_cur);

  std::vector<double> alpha;
  std::vector<double> beta;
  alpha.reserve(static_cast<std::size_t>(krylov_dim));
  beta.reserve(static_cast<std::size_t>(krylov_dim));

  int actual_steps = 0;
  double beta_prev = 0.0;
  for (int step = 0; step < krylov_dim; ++step) {
    apply_hamiltonian(block, q_cur, q_next);
    const double a = dot(q_cur, q_next).real();
    for (std::size_t i = 0; i < block.basis_dim; ++i) {
      q_next[i] -= a * q_cur[i] + beta_prev * q_prev[i];
    }
    alpha.push_back(a);
    const double b = std::sqrt(complex_norm2(q_next));
    ++actual_steps;
    if (step_callback) step_callback(actual_steps);
    if (b < 1e-12 || step + 1 == krylov_dim) break;
    beta.push_back(b);
    q_prev.swap(q_cur);
    q_cur.swap(q_next);
    for (Complex& value : q_cur) value /= b;
    beta_prev = b;
  }

  std::vector<std::vector<double>> tri(
      static_cast<std::size_t>(actual_steps),
      std::vector<double>(static_cast<std::size_t>(actual_steps), 0.0));
  for (int i = 0; i < actual_steps; ++i) {
    tri[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] =
        alpha[static_cast<std::size_t>(i)];
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

  LanczosSample result;
  result.eigenvalues.reserve(static_cast<std::size_t>(actual_steps));
  result.overlaps.reserve(static_cast<std::size_t>(actual_steps));
  for (int sorted = 0; sorted < actual_steps; ++sorted) {
    const int original = order[static_cast<std::size_t>(sorted)];
    const double v1 = eigenvectors[0][static_cast<std::size_t>(original)];
    result.eigenvalues.push_back(eigenvalues[static_cast<std::size_t>(original)]);
    result.overlaps.push_back(v1 * v1);
  }
  return result;
}

LanczosSpectrum run_ftlm_block(
    const MomentumBlock& block,
    const Params& params,
    std::uint64_t block_seed) {
  LanczosSpectrum spectrum;
  spectrum.basis_dim = block.basis_dim;
  spectrum.particles = block.particles;
  std::vector<LanczosSample> samples(static_cast<std::size_t>(params.samples));
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) num_threads(params.threads > 0 ? params.threads : omp_get_max_threads())
#endif
  for (int sample = 0; sample < params.samples; ++sample) {
    samples[static_cast<std::size_t>(sample)] =
        run_ftlm_sample(block, params, block_seed, sample);
  }
  for (const LanczosSample& sample : samples) {
    spectrum.eigenvalues.insert(
        spectrum.eigenvalues.end(), sample.eigenvalues.begin(), sample.eigenvalues.end());
    for (double overlap : sample.overlaps) {
      spectrum.sample_weights.push_back(overlap / static_cast<double>(params.samples));
    }
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
    const std::vector<double>& partitions,
    const std::vector<double>& log_partitions) {
  std::ofstream out(path);
  if (!out) {
    die("Failed to open output file: " + path);
  }
  out << "mu,n,charge_correlation,compressibility,partition_like,log_partition\n";
  out << std::setprecision(15);
  for (std::size_t i = 0; i < mu_values.size(); ++i) {
    out << mu_values[i] << "," << densities[i] << ","
        << charge_correlations[i] << "," << compressibilities[i] << ","
        << partitions[i] << "," << log_partitions[i] << "\n";
  }
}

void write_multi_beta_results(
    const std::string& path,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    const std::vector<std::vector<double>>& densities,
    const std::vector<std::vector<double>>& charge_correlations,
    const std::vector<std::vector<double>>& compressibilities,
    const std::vector<std::vector<double>>& partitions,
    const std::vector<std::vector<double>>& log_partitions) {
  std::ofstream out(path);
  if (!out) {
    die("Failed to open output file: " + path);
  }
  out << "beta,mu,n,charge_correlation,compressibility,partition_like,log_partition\n";
  out << std::setprecision(15);
  for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
    for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
      out << beta_values[ibeta] << "," << mu_values[imu] << ","
          << densities[ibeta][imu] << ","
          << charge_correlations[ibeta][imu] << ","
          << compressibilities[ibeta][imu] << ","
          << partitions[ibeta][imu] << ","
          << log_partitions[ibeta][imu] << "\n";
    }
  }
}

std::atomic<bool> g_stop_requested{false};
std::atomic<bool> g_signal_requested{false};

void stop_signal_handler(int) {
  g_signal_requested.store(true, std::memory_order_relaxed);
  g_stop_requested.store(true, std::memory_order_relaxed);
}

double peak_rss_gb() {
  struct rusage usage {};
  if (getrusage(RUSAGE_SELF, &usage) != 0) return 0.0;
#if defined(__APPLE__)
  return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0 * 1024.0);
#else
  return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#endif
}

std::string iso_timestamp() {
  const std::time_t now = std::time(nullptr);
  std::tm tm{};
#if defined(_WIN32)
  gmtime_s(&tm, &now);
#else
  gmtime_r(&now, &tm);
#endif
  std::ostringstream out;
  out << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
  return out.str();
}

std::string environment_value(const char* name) {
  const char* value = std::getenv(name);
  return value == nullptr ? std::string{} : std::string(value);
}

class ProgressReporter {
 public:
  ProgressReporter(
      const Params& params,
      const std::vector<ftlm_checkpoint::ManifestEntry>& manifest,
      ftlm_checkpoint::CheckpointData& checkpoint,
      std::mutex& checkpoint_mutex)
      : params_(params),
        manifest_(manifest),
        checkpoint_(checkpoint),
        checkpoint_mutex_(checkpoint_mutex),
        start_(std::chrono::steady_clock::now()) {
    {
      std::lock_guard<std::mutex> lock(checkpoint_mutex_);
      initial_status_ = ftlm_checkpoint::completion_status(
          checkpoint_, manifest_, params_.samples);
    }
    monitor_ = std::thread([this]() { monitor_loop(); });
  }

  ~ProgressReporter() {
    stop();
  }

  void set_block(const BlockKey& key, std::uint64_t dimension) {
    {
      std::lock_guard<std::mutex> lock(state_mutex_);
      current_block_ = key;
      current_dimension_ = dimension;
      has_current_block_ = true;
      active_steps_.clear();
    }
    emit("BLOCK_START", -1);
  }

  void clear_block() {
    // Emit while the block is still present in the snapshot so the durable
    // completion event identifies exactly which unit just finished.
    emit("BLOCK_COMPLETE", -1);
    {
      std::lock_guard<std::mutex> lock(state_mutex_);
      active_steps_.clear();
      active_started_.clear();
      has_current_block_ = false;
      current_dimension_ = 0;
    }
  }

  void start_sample(int sample_id) {
    std::lock_guard<std::mutex> lock(state_mutex_);
    active_steps_[sample_id] = 0;
    active_started_[sample_id] = std::chrono::steady_clock::now();
  }

  void update_sample_step(int sample_id, int step) {
    std::lock_guard<std::mutex> lock(state_mutex_);
    active_steps_[sample_id] = step;
  }

  void finish_sample(int sample_id) {
    {
      std::lock_guard<std::mutex> lock(state_mutex_);
      active_steps_.erase(sample_id);
      const auto started = active_started_.find(sample_id);
      last_sample_duration_seconds_ = started == active_started_.end()
          ? -1.0
          : std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started->second).count();
      if (started != active_started_.end()) active_started_.erase(started);
      ++new_samples_;
      if (has_current_block_) {
        has_last_durable_sample_ = true;
        last_durable_block_ = current_block_;
        last_durable_sample_id_ = sample_id;
      }
    }
    emit("CHECKPOINTED", sample_id);
  }

  void event(const std::string& name) {
    emit(name, -1);
  }

  bool runtime_expired() const {
    if (params_.max_runtime_minutes <= 0.0) return false;
    return elapsed_seconds() >= params_.max_runtime_minutes * 60.0;
  }

  double elapsed_seconds() const {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_).count();
  }

  void stop() {
    bool expected = false;
    if (!stopped_.compare_exchange_strong(expected, true)) return;
    monitor_cv_.notify_all();
    if (monitor_.joinable()) monitor_.join();
  }

 private:
  struct StateSnapshot {
    bool has_block = false;
    BlockKey block;
    std::uint64_t dimension = 0;
    std::map<int, int> active_steps;
    std::uint64_t new_samples = 0;
    bool has_last_durable_sample = false;
    BlockKey last_durable_block;
    int last_durable_sample_id = -1;
    double last_sample_duration_seconds = -1.0;
  };

  StateSnapshot state_snapshot() const {
    std::lock_guard<std::mutex> lock(state_mutex_);
    return {
        has_current_block_, current_block_, current_dimension_, active_steps_, new_samples_,
        has_last_durable_sample_, last_durable_block_, last_durable_sample_id_,
        last_sample_duration_seconds_};
  }

  void monitor_loop() {
    auto last_heartbeat = std::chrono::steady_clock::now();
    bool stop_reported = false;
    while (!stopped_.load()) {
      {
        std::unique_lock<std::mutex> lock(monitor_wait_mutex_);
        monitor_cv_.wait_for(
            lock, std::chrono::seconds(1), [this]() { return stopped_.load(); });
      }
      // Do not overwrite a terminal RUN_COMPLETE/RUN_INCOMPLETE state with a
      // heartbeat while the owner is joining this monitor thread.
      if (stopped_.load()) break;
      if (runtime_expired()) {
        g_stop_requested.store(true);
        if (!stop_reported) {
          emit("STOP_REQUESTED_RUNTIME", -1);
          stop_reported = true;
        }
      } else if (g_signal_requested.load() && !stop_reported) {
        emit("STOP_REQUESTED_SIGNAL", -1);
        stop_reported = true;
      }
      const auto now = std::chrono::steady_clock::now();
      if (std::chrono::duration_cast<std::chrono::seconds>(now - last_heartbeat).count() >=
          params_.progress_interval_seconds) {
        emit("HEARTBEAT", -1);
        last_heartbeat = now;
      }
    }
  }

  void emit(const std::string& event_name, int sample_id) {
    ftlm_checkpoint::CompletionStatus status;
    {
      std::lock_guard<std::mutex> lock(checkpoint_mutex_);
      status = ftlm_checkpoint::completion_status(checkpoint_, manifest_, params_.samples);
    }
    const StateSnapshot state = state_snapshot();
    long double active_weight = 0.0L;
    if (state.has_block) {
      for (const auto& item : state.active_steps) {
        active_weight += static_cast<long double>(state.dimension) * item.second;
      }
    }
    const double live_fraction = status.total_weight > 0.0L
        ? static_cast<double>((status.durable_weight + active_weight) / status.total_weight)
        : status.weighted_fraction;
    const double elapsed = elapsed_seconds();
    const double progress_delta = live_fraction - initial_status_.weighted_fraction;
    const double eta_hours = progress_delta > 1e-12
        ? elapsed * std::max(0.0, 1.0 - live_fraction) / progress_delta / 3600.0
        : -1.0;
    std::uintmax_t checkpoint_bytes = 0;
    if (!params_.checkpoint.empty() && std::filesystem::exists(params_.checkpoint)) {
      checkpoint_bytes = std::filesystem::file_size(params_.checkpoint);
    }
    std::ostringstream active;
    bool first = true;
    for (const auto& item : state.active_steps) {
      if (!first) active << ';';
      first = false;
      active << item.first << ':' << item.second << '/' << params_.lanczos_steps;
    }
    if (first) active << '-';

    std::ostringstream human;
    human << std::fixed << std::setprecision(2)
          << "FTLM_PROGRESS event=" << event_name
          << " timestamp=" << iso_timestamp()
          << " twist=" << (params_.twist_id.empty() ? "-" : params_.twist_id)
          << " target_R=" << params_.samples
          << " elapsed_s=" << elapsed
          << " min_complete_R=" << status.minimum_complete_samples
          << " samples=" << status.durable_sample_records << '/'
          << status.expected_sample_records
          << " blocks=" << status.completed_blocks << '/' << status.expected_blocks
          << " weighted_progress=" << 100.0 * live_fraction
          << " checkpoint_mb=" << static_cast<double>(checkpoint_bytes) / (1024.0 * 1024.0)
          << " max_rss_gb=" << peak_rss_gb()
          << " reused_samples=" << initial_status_.durable_sample_records
          << " new_samples=" << state.new_samples;
    if (eta_hours >= 0.0) human << " eta_hours=" << eta_hours;
    if (state.has_block) {
      human << " sector=" << state.block.n_up << ',' << state.block.n_down
            << " block=" << state.block.mx << ',' << state.block.my
            << " dim=" << state.dimension
            << " active_samples=" << active.str();
    }
    if (sample_id >= 0) human << " sample=" << sample_id;
    if (sample_id >= 0 && state.last_sample_duration_seconds >= 0.0) {
      human << " sample_seconds=" << state.last_sample_duration_seconds;
    }

    std::ostringstream json;
    json << std::setprecision(15)
         << '{'
         << "\"event\":\"" << event_name << "\","
         << "\"timestamp\":\"" << iso_timestamp() << "\","
         << "\"twist\":\"" << (params_.twist_id.empty() ? "" : params_.twist_id) << "\","
         << "\"target_R\":" << params_.samples << ','
         << "\"elapsed_seconds\":" << elapsed << ','
         << "\"minimum_complete_R\":" << status.minimum_complete_samples << ','
         << "\"durable_sample_records\":" << status.durable_sample_records << ','
         << "\"expected_sample_records\":" << status.expected_sample_records << ','
         << "\"completed_blocks\":" << status.completed_blocks << ','
         << "\"expected_blocks\":" << status.expected_blocks << ','
         << "\"weighted_progress\":" << live_fraction << ','
         << "\"checkpoint_bytes\":" << checkpoint_bytes << ','
         << "\"max_rss_gb\":" << peak_rss_gb() << ','
         << "\"reused_samples\":" << initial_status_.durable_sample_records << ','
         << "\"new_samples\":" << state.new_samples << ','
         << "\"eta_hours\":" << eta_hours << ','
         << "\"slurm_job_id\":\"" << environment_value("SLURM_JOB_ID") << "\","
         << "\"node\":\"" << environment_value("SLURMD_NODENAME") << "\","
         << "\"active_samples\":\"" << active.str() << "\","
         << "\"last_sample_duration_seconds\":"
         << state.last_sample_duration_seconds << ','
         << "\"current_block\":";
    if (state.has_block) {
      json << '{'
           << "\"n_up\":" << state.block.n_up << ','
           << "\"n_down\":" << state.block.n_down << ','
           << "\"mx\":" << state.block.mx << ','
           << "\"my\":" << state.block.my << ','
           << "\"basis_dim\":" << state.dimension << '}';
    } else {
      json << "null";
    }
    json << ",\"checkpointed_sample\":";
    if (sample_id >= 0) json << sample_id; else json << "null";
    json << ",\"last_durable_sample\":";
    if (state.has_last_durable_sample) {
      json << '{'
           << "\"n_up\":" << state.last_durable_block.n_up << ','
           << "\"n_down\":" << state.last_durable_block.n_down << ','
           << "\"mx\":" << state.last_durable_block.mx << ','
           << "\"my\":" << state.last_durable_block.my << ','
           << "\"sample_id\":" << state.last_durable_sample_id << '}';
    } else {
      json << "null";
    }
    json << '}';

    std::lock_guard<std::mutex> output_lock(output_mutex_);
    std::cout << human.str() << std::endl;
    if (!params_.progress_jsonl.empty()) {
      std::ofstream out(params_.progress_jsonl, std::ios::app);
      if (out) out << json.str() << '\n';
    }
    if (!params_.progress_state.empty()) {
      const std::string temporary =
          params_.progress_state + ".tmp." + std::to_string(static_cast<long long>(::getpid()));
      {
        std::ofstream out(temporary);
        if (out) out << json.str() << '\n';
      }
      std::error_code error;
      std::filesystem::rename(temporary, params_.progress_state, error);
      if (error) {
        std::filesystem::remove(params_.progress_state, error);
        error.clear();
        std::filesystem::rename(temporary, params_.progress_state, error);
      }
    }
  }

  const Params& params_;
  const std::vector<ftlm_checkpoint::ManifestEntry>& manifest_;
  ftlm_checkpoint::CheckpointData& checkpoint_;
  std::mutex& checkpoint_mutex_;
  std::chrono::steady_clock::time_point start_;
  ftlm_checkpoint::CompletionStatus initial_status_;
  mutable std::mutex state_mutex_;
  std::mutex output_mutex_;
  bool has_current_block_ = false;
  BlockKey current_block_;
  std::uint64_t current_dimension_ = 0;
  std::map<int, int> active_steps_;
  std::map<int, std::chrono::steady_clock::time_point> active_started_;
  std::uint64_t new_samples_ = 0;
  bool has_last_durable_sample_ = false;
  BlockKey last_durable_block_;
  int last_durable_sample_id_ = -1;
  double last_sample_duration_seconds_ = -1.0;
  std::atomic<bool> stopped_{false};
  std::mutex monitor_wait_mutex_;
  std::condition_variable monitor_cv_;
  std::thread monitor_;
};

int run_v2_checkpoint_job(
    const Params& params,
    const Lattice& lattice,
    const std::vector<double>& mu_values) {
  if (params.checkpoint.empty()) die("v2 checkpoint execution requires --checkpoint.");
  g_stop_requested.store(false);
  g_signal_requested.store(false);
#ifdef SIGUSR1
  std::signal(SIGUSR1, stop_signal_handler);
#endif
  std::signal(SIGTERM, stop_signal_handler);

  ftlm_checkpoint::CheckpointData checkpoint =
      ftlm_checkpoint::read_checkpoint(params.checkpoint, false);
  if (checkpoint.format != ftlm_checkpoint::Format::kV2) {
    die("Internal error: v2 runner received a non-v2 checkpoint.");
  }
  if (checkpoint.trailing_partial_record) {
    std::cout << "FTLM_PROGRESS event=REPAIR_TRAILING_RECORD valid_bytes="
              << checkpoint.valid_bytes << std::endl;
    ftlm_checkpoint::truncate_to_valid_bytes(params.checkpoint, checkpoint.valid_bytes);
    checkpoint = ftlm_checkpoint::read_checkpoint(params.checkpoint, false);
  }

  const auto manifest = ftlm_checkpoint::build_manifest(
      params.lx, params.ly, params.exact_block_threshold, params.lanczos_steps);
  std::map<BlockKey, const ftlm_checkpoint::ManifestEntry*> manifest_by_key;
  for (const auto& entry : manifest) manifest_by_key[entry.key] = &entry;

  std::mutex checkpoint_mutex;
  ftlm_checkpoint::Writer writer(params.checkpoint);
  ProgressReporter progress(params, manifest, checkpoint, checkpoint_mutex);
  const auto is_block_complete = [&](const ftlm_checkpoint::ManifestEntry& entry) {
    std::lock_guard<std::mutex> lock(checkpoint_mutex);
    return ftlm_checkpoint::block_complete(checkpoint, entry, params.samples);
  };
  {
    std::lock_guard<std::mutex> lock(checkpoint_mutex);
    std::cout << "FTLM_PROGRESS event=RESUME "
              << ftlm_checkpoint::status_text(
                     ftlm_checkpoint::completion_status(checkpoint, manifest, params.samples))
              << std::endl;
  }

  bool stopped_early = false;
  bool selected_block_found = false;
  for (int n_up = 0; n_up <= lattice.sites && !stopped_early; ++n_up) {
    for (int n_down = 0; n_down <= lattice.sites && !stopped_early; ++n_down) {
      if (params.only_block &&
          (n_up != params.only_block_nup || n_down != params.only_block_ndown)) {
        continue;
      }
      const SectorKey sector_key{n_up, n_down};
      if (!params.only_block) {
        int completed_samples = 0;
        {
          std::lock_guard<std::mutex> lock(checkpoint_mutex);
          const auto marker = checkpoint.sector_complete_samples.find(sector_key);
          if (marker != checkpoint.sector_complete_samples.end()) {
            completed_samples = marker->second;
          }
        }
        if (completed_samples >= params.samples) {
          std::cout << "FTLM_PROGRESS event=SECTOR_REUSED sector="
                    << n_up << ',' << n_down << " completed_R=" << completed_samples
                    << std::endl;
          continue;
        }
      }

      std::vector<const ftlm_checkpoint::ManifestEntry*> sector_entries;
      std::uint64_t sector_dimension = 0;
      bool sector_needs_work = false;
      for (const auto& entry : manifest) {
        if (entry.key.n_up != n_up || entry.key.n_down != n_down) continue;
        sector_dimension += entry.basis_dim;
        if (params.only_block &&
            (entry.key.mx != params.only_block_mx || entry.key.my != params.only_block_my)) {
          continue;
        }
        sector_entries.push_back(&entry);
        if (!is_block_complete(entry)) {
          sector_needs_work = true;
        }
      }
      if (sector_entries.empty()) continue;
      if (sector_dimension > params.max_sector_dim) {
        std::ostringstream message;
        message << "Sector (Nup=" << n_up << ", Ndown=" << n_down
                << ") has full dimension " << sector_dimension
                << ", which exceeds --max-sector-dim=" << params.max_sector_dim;
        die(message.str());
      }
      if (!sector_needs_work) {
        if (!params.only_block) {
          writer.append_sector_complete(sector_key, params.samples);
          std::lock_guard<std::mutex> lock(checkpoint_mutex);
          checkpoint.sector_complete_samples[sector_key] = params.samples;
        }
        continue;
      }
      if (g_stop_requested.load() || progress.runtime_expired()) {
        stopped_early = true;
        break;
      }

      std::cout << "FTLM_PROGRESS event=SECTOR_START sector=" << n_up << ',' << n_down
                << " full_dim=" << sector_dimension << std::endl;
      const SectorBasis sector = build_sector_basis(lattice, n_up, n_down, params);
      if (sector.full_dim != sector_dimension) {
        die("Analytic and constructed sector dimensions disagree.");
      }

      for (int my = 0; my < lattice.ly && !stopped_early; ++my) {
        for (int mx = 0; mx < lattice.lx && !stopped_early; ++mx) {
          if (params.only_block &&
              (mx != params.only_block_mx || my != params.only_block_my)) {
            continue;
          }
          const BlockKey key{n_up, n_down, mx, my};
          const auto manifest_it = manifest_by_key.find(key);
          if (manifest_it == manifest_by_key.end()) continue;
          const ftlm_checkpoint::ManifestEntry& entry = *manifest_it->second;
          selected_block_found = selected_block_found || params.only_block;
          if (is_block_complete(entry)) {
            continue;
          }
          if (g_stop_requested.load() || progress.runtime_expired()) {
            stopped_early = true;
            break;
          }

          const MomentumBlock block = build_momentum_block(sector, lattice, mx, my);
          if (block.basis_dim != entry.basis_dim) {
            die("Analytic and constructed momentum-block dimensions disagree for " +
                ftlm_checkpoint::block_key_text(key));
          }
          progress.set_block(key, block.basis_dim);

          if (entry.exact) {
            ftlm_checkpoint::ExactRecord record;
            record.key = key;
            record.particles = block.particles;
            record.basis_dim = block.basis_dim;
            record.eigenvalues = diagonalize_block_exact(block);
            writer.append_exact(record);
            {
              std::lock_guard<std::mutex> lock(checkpoint_mutex);
              record.eigenvalues.clear();
              checkpoint.exact_blocks[key] = record;
            }
            progress.event("CHECKPOINTED_EXACT");
          } else {
            std::vector<int> missing_samples;
            {
              std::lock_guard<std::mutex> lock(checkpoint_mutex);
              const auto existing = checkpoint.samples.find(key);
              for (int sample_id = 0; sample_id < params.samples; ++sample_id) {
                if (existing == checkpoint.samples.end() ||
                    existing->second.find(sample_id) == existing->second.end()) {
                  missing_samples.push_back(sample_id);
                }
              }
            }
            std::exception_ptr thread_error;
            std::mutex thread_error_mutex;
            const std::uint64_t block_seed =
                make_block_seed(params, n_up, n_down, mx, my);
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic, 1) num_threads(params.threads > 0 ? params.threads : omp_get_max_threads())
#endif
            for (int index = 0; index < static_cast<int>(missing_samples.size()); ++index) {
              if (g_stop_requested.load()) continue;
              const int sample_id = missing_samples[static_cast<std::size_t>(index)];
              progress.start_sample(sample_id);
              try {
                LanczosSample sample = run_ftlm_sample(
                    block,
                    params,
                    block_seed,
                    sample_id,
                    [&](int step) { progress.update_sample_step(sample_id, step); });
                ftlm_checkpoint::SampleRecord record;
                record.key = key;
                record.particles = block.particles;
                record.basis_dim = block.basis_dim;
                record.sample_id = sample_id;
                record.eigenvalues = std::move(sample.eigenvalues);
                record.overlaps = std::move(sample.overlaps);
                writer.append_sample(record);
                {
                  std::lock_guard<std::mutex> lock(checkpoint_mutex);
                  record.eigenvalues.clear();
                  record.overlaps.clear();
                  checkpoint.samples[key][sample_id] = record;
                }
                progress.finish_sample(sample_id);
              } catch (...) {
                {
                  std::lock_guard<std::mutex> lock(thread_error_mutex);
                  if (!thread_error) thread_error = std::current_exception();
                }
                g_stop_requested.store(true);
              }
            }
            if (thread_error) std::rethrow_exception(thread_error);
          }

          if (is_block_complete(entry)) {
            writer.append_block_complete(key, params.samples);
            {
              std::lock_guard<std::mutex> lock(checkpoint_mutex);
              checkpoint.block_complete_samples[key] = params.samples;
            }
          }
          progress.clear_block();
          if (g_stop_requested.load() || progress.runtime_expired()) {
            stopped_early = true;
          }
        }
      }

      if (!params.only_block && !stopped_early) {
        bool sector_complete = true;
        for (const auto* entry : sector_entries) {
          if (!is_block_complete(*entry)) {
            sector_complete = false;
            break;
          }
        }
        if (sector_complete) {
          writer.append_sector_complete(sector_key, params.samples);
          {
            std::lock_guard<std::mutex> lock(checkpoint_mutex);
            checkpoint.sector_complete_samples[sector_key] = params.samples;
          }
          progress.event("SECTOR_COMPLETE");
        }
      }
    }
  }

  if (params.only_block && !selected_block_found) {
    die("--only-block target has zero dimension or was not found.");
  }

  ftlm_checkpoint::CompletionStatus status;
  {
    std::lock_guard<std::mutex> lock(checkpoint_mutex);
    status = ftlm_checkpoint::completion_status(checkpoint, manifest, params.samples);
  }
  if (params.only_block) {
    const BlockKey selected{
        params.only_block_nup, params.only_block_ndown,
        params.only_block_mx, params.only_block_my};
    const auto entry = manifest_by_key.find(selected);
    const bool complete = entry != manifest_by_key.end() &&
        is_block_complete(*entry->second);
    progress.event(complete ? "ONLY_BLOCK_COMPLETE" : "ONLY_BLOCK_INCOMPLETE");
    progress.stop();
    std::cout << "FTLM_STATUS " << ftlm_checkpoint::status_text(status) << std::endl;
    return complete ? 0 : 75;
  }

  if (!status.complete) {
    progress.event("RUN_INCOMPLETE");
    progress.stop();
    std::cout << "FTLM_STATUS " << ftlm_checkpoint::status_text(status) << std::endl;
    return 75;
  }

  bool append_run_marker = false;
  {
    std::lock_guard<std::mutex> lock(checkpoint_mutex);
    append_run_marker = checkpoint.run_complete_samples < params.samples;
  }
  if (append_run_marker) {
    writer.append_run_complete(params.samples);
    std::lock_guard<std::mutex> lock(checkpoint_mutex);
    checkpoint.run_complete_samples = params.samples;
  }
  progress.event("RUN_COMPLETE");
  progress.stop();
  const ftlm_checkpoint::CheckpointData reduction_checkpoint =
      ftlm_checkpoint::read_checkpoint(params.checkpoint, true);
  const auto grid = ftlm_checkpoint::reduce_checkpoint(
      reduction_checkpoint,
      manifest,
      params.samples,
      lattice.sites,
      params.beta_values,
      mu_values);
  ftlm_checkpoint::write_thermo_csv(
      params.output, params.beta_values, mu_values, grid);
  std::cout << "Wrote " << params.output << std::endl;
  return 0;
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
    if ((params.only_block || params.max_runtime_minutes > 0.0 ||
         !params.progress_jsonl.empty() || !params.progress_state.empty()) &&
        params.checkpoint.empty()) {
      die("Runtime/progress/only-block operation requires --checkpoint.");
    }
    const ftlm_checkpoint::Format checkpoint_format =
        validate_or_write_checkpoint_metadata(params);
    if (checkpoint_format == ftlm_checkpoint::Format::kV2) {
      return run_v2_checkpoint_job(params, lattice, mu_values);
    }
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
      std::vector<std::vector<double>> all_log_partitions(
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
          all_log_partitions[ibeta][imu] =
              std::log(z_scaled[imu]) - beta * emin_by_mu[imu];
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
          all_partitions,
          all_log_partitions);
    } else {
      std::vector<double> densities(mu_values.size(), 0.0);
      std::vector<double> charge_correlations(mu_values.size(), 0.0);
      std::vector<double> compressibilities(mu_values.size(), 0.0);
      std::vector<double> partitions(mu_values.size(), 0.0);
      std::vector<double> log_partitions(mu_values.size(), 0.0);
      for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
        partitions[imu] = single_z_scaled[imu];
        log_partitions[imu] =
            std::log(single_z_scaled[imu]) - params.beta * single_emin_by_mu[imu];
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
          partitions,
          log_partitions);
    }
    std::cout << "Wrote " << params.output << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }
}
