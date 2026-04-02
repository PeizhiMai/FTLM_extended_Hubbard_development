#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
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
  double mu_min = -4.0;
  double mu_max = 10.0;
  int mu_count = 61;
  int samples = 5;
  int lanczos_steps = 80;
  std::uint64_t seed = 12345;
  std::size_t max_sector_dim = 2000000;
  std::string output = "n_vs_mu.csv";
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
  double min_shifted_energy = std::numeric_limits<double>::infinity();
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
      << "  --mu-min X            minimum chemical potential\n"
      << "  --mu-max X            maximum chemical potential\n"
      << "  --mu-count N          number of mu points\n"
      << "  --samples N           FTLM random vectors per k block\n"
      << "  --lanczos-steps N     Lanczos steps per random vector\n"
      << "  --seed N              random seed\n"
      << "  --max-sector-dim N    abort if a full sector exceeds this size\n"
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
    } else if (arg == "--seed") {
      p.seed = parse_value<std::uint64_t>(next(arg));
    } else if (arg == "--max-sector-dim") {
      p.max_sector_dim = parse_value<std::size_t>(next(arg));
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
  if (p.beta <= 0.0) {
    die("--beta must be positive.");
  }
  if (p.samples <= 0 || p.lanczos_steps <= 0) {
    die("--samples and --lanczos-steps must be positive.");
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
  const std::uint64_t occupied = state.up | state.down;
  const int doublons = popcount64(state.up & state.down);
  int nn_pairs = 0;
  for (const auto& bond : lattice.unique_bonds) {
    const int ni = ((occupied >> bond.first) & 1ULL) ? 1 : 0;
    const int nj = ((occupied >> bond.second) & 1ULL) ? 1 : 0;
    nn_pairs += ni * nj;
  }
  return u * static_cast<double>(doublons) + v * static_cast<double>(nn_pairs);
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
    std::mt19937_64& rng) {
  LanczosSpectrum spectrum;
  spectrum.basis_dim = block.basis_dim;
  spectrum.particles = block.particles;

  const int krylov_dim =
      std::min<int>(params.lanczos_steps, static_cast<int>(block.basis_dim));
  std::normal_distribution<double> normal(0.0, 1.0);

  std::vector<Complex> q_prev(block.basis_dim, 0.0);
  std::vector<Complex> q_cur(block.basis_dim, 0.0);
  std::vector<Complex> q_next(block.basis_dim, 0.0);
  std::vector<Complex> hq(block.basis_dim, 0.0);

  for (int sample = 0; sample < params.samples; ++sample) {
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

    for (int sorted = 0; sorted < actual_steps; ++sorted) {
      const int original = order[static_cast<std::size_t>(sorted)];
      const double v1 = eigenvectors[0][static_cast<std::size_t>(original)];
      spectrum.eigenvalues.push_back(eigenvalues[static_cast<std::size_t>(original)]);
      spectrum.sample_weights.push_back((v1 * v1) / static_cast<double>(params.samples));
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
    thermo.z += weight;
    thermo.n_total += static_cast<double>(spectrum.particles) * weight;
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
    thermo.z += weight;
    thermo.n_total += static_cast<double>(particles) * weight;
  }
  return thermo;
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
  out << "nup,ndown,mx,my,basis_dim,mu,z_ftlm,z_ed,z_ratio,n_ftlm_total,n_ed_total,n_diff,min_ftlm,min_ed,sum_sample_weights\n";
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
        << ftlm.min_shifted_energy << "," << exact.min_shifted_energy << ","
        << sum_sample_weights << "\n";
  }
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
    const std::vector<double>& partitions) {
  std::ofstream out(path);
  if (!out) {
    die("Failed to open output file: " + path);
  }
  out << "mu,n,partition_like\n";
  out << std::setprecision(15);
  for (std::size_t i = 0; i < mu_values.size(); ++i) {
    out << mu_values[i] << "," << densities[i] << "," << partitions[i] << "\n";
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Params params = parse_args(argc, argv);
    const Lattice lattice = build_lattice(params.lx, params.ly);
    const std::vector<double> mu_values =
        linspace(params.mu_min, params.mu_max, params.mu_count);
    const bool trace_blocks = !params.trace_partition_csv.empty();
    const bool debug_selected_block =
        params.debug_block_nup >= 0 && params.debug_block_ndown >= 0 &&
        params.debug_block_mx >= 0 && params.debug_block_my >= 0;

    std::cout << "Reference used: "
              << "reference_ftlm_hub_cond/fltm_hub_cond/hubTri2Dcond_omp.f"
              << " and cond_spect_omp.f\n";
    std::cout << "Running FTLM with translation-symmetry k blocks on a "
              << params.lx << "x" << params.ly << " rectangular lattice ("
              << lattice.sites << " sites)\n";
    std::cout << "beta=" << params.beta << " U=" << params.u << " V=" << params.v
              << " tx=" << params.tx << " ty=" << params.ty << " tp=" << params.tp
              << " phix=" << params.phix << " phiy=" << params.phiy << "\n";

    std::mt19937_64 rng(params.seed);
    std::vector<double> densities(mu_values.size(), 0.0);
    std::vector<double> partitions(mu_values.size(), 0.0);
    std::vector<double> emin_by_mu(
        mu_values.size(), std::numeric_limits<double>::infinity());
    std::vector<double> z_scaled(mu_values.size(), 0.0);
    std::vector<double> n_scaled(mu_values.size(), 0.0);
    std::ofstream trace_out;
    if (trace_blocks) {
      trace_out.open(params.trace_partition_csv);
      if (!trace_out) {
        die("Failed to open trace CSV: " + params.trace_partition_csv);
      }
      trace_out << "nup,ndown,mx,my,particles,basis_dim,mu,"
                << "z_ftlm,z_ed,z_ratio,n_ftlm_total,n_ed_total,n_total_diff,"
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

            const LanczosSpectrum spectrum = run_ftlm_block(block, params, rng);
            const bool need_exact_block =
                trace_blocks || (debug_selected_block &&
                                 block_matches_debug(params, n_up, n_down, mx, my));
            std::vector<double> exact_eigenvalues;
            if (need_exact_block) {
              exact_eigenvalues = diagonalize_block_exact(block);
            }

            if (trace_blocks) {
              const double sum_sample_weights =
                  std::accumulate(spectrum.sample_weights.begin(), spectrum.sample_weights.end(), 0.0);
              for (double mu : mu_values) {
                const BlockThermo ftlm = raw_block_thermo_from_ftlm(spectrum, params.beta, mu);
                const BlockThermo exact =
                    raw_block_thermo_from_exact(exact_eigenvalues, spectrum.particles, params.beta, mu);
                const double z_ratio = (exact.z > 0.0) ? (ftlm.z / exact.z) : 0.0;
                trace_out << n_up << "," << n_down << "," << mx << "," << my << ","
                          << spectrum.particles << "," << spectrum.basis_dim << ","
                          << mu << "," << ftlm.z << "," << exact.z << "," << z_ratio
                          << "," << ftlm.n_total << "," << exact.n_total << ","
                          << (ftlm.n_total - exact.n_total) << ","
                          << ftlm.min_shifted_energy << "," << exact.min_shifted_energy << ","
                          << sum_sample_weights << "\n";
              }
            }

            if (debug_selected_block && block_matches_debug(params, n_up, n_down, mx, my)) {
              const double sum_sample_weights =
                  std::accumulate(spectrum.sample_weights.begin(), spectrum.sample_weights.end(), 0.0);
              std::cout << "Selected block debug Nup=" << n_up
                        << " Ndown=" << n_down
                        << " mx=" << mx
                        << " my=" << my
                        << " basis_dim=" << spectrum.basis_dim
                        << " sum_sample_weights=" << sum_sample_weights << "\n";
              for (double mu : mu_values) {
                const BlockThermo ftlm = raw_block_thermo_from_ftlm(spectrum, params.beta, mu);
                const BlockThermo exact =
                    raw_block_thermo_from_exact(exact_eigenvalues, spectrum.particles, params.beta, mu);
                const double z_ratio = (exact.z > 0.0) ? (ftlm.z / exact.z) : 0.0;
                std::cout << "  debug mu=" << std::setw(12) << mu
                          << " z_ftlm=" << std::setw(14) << ftlm.z
                          << " z_ed=" << std::setw(14) << exact.z
                          << " ratio=" << std::setw(14) << z_ratio
                          << " n_diff=" << std::setw(14) << (ftlm.n_total - exact.n_total)
                          << "\n";
              }
              if (!params.debug_block_csv.empty()) {
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
              }
            }

            for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
              const double mu = mu_values[imu];
              double block_min = std::numeric_limits<double>::infinity();
              for (double eigenvalue : spectrum.eigenvalues) {
                block_min = std::min(
                    block_min,
                    eigenvalue - mu * static_cast<double>(spectrum.particles));
              }
              if (!std::isfinite(block_min)) {
                continue;
              }

              if (!std::isfinite(emin_by_mu[imu])) {
                emin_by_mu[imu] = block_min;
              } else if (block_min < emin_by_mu[imu]) {
                const double scale =
                    std::exp(-params.beta * (emin_by_mu[imu] - block_min));
                z_scaled[imu] *= scale;
                n_scaled[imu] *= scale;
                emin_by_mu[imu] = block_min;
              }

              for (std::size_t i = 0; i < spectrum.eigenvalues.size(); ++i) {
                const double shifted =
                    spectrum.eigenvalues[i] - mu * static_cast<double>(spectrum.particles) -
                    emin_by_mu[imu];
                const double weight =
                    static_cast<double>(spectrum.basis_dim) * spectrum.sample_weights[i] *
                    std::exp(-params.beta * shifted);
                z_scaled[imu] += weight;
                n_scaled[imu] += static_cast<double>(spectrum.particles) * weight;
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

    if (debug_selected_block && !params.debug_block_csv.empty() && !wrote_selected_block_debug) {
      die("Selected block debug target was not found among active FTLM blocks.");
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
