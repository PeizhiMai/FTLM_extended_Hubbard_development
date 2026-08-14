#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ftlm_checkpoint {

constexpr int kCheckpointVersion = 2;
constexpr int kDepthCheckpointVersion = 3;
constexpr int kRngVersion = 1;
constexpr const char* kRngAlgorithm = "mt19937_64_box_muller53_v1";
constexpr const char* kSeedDerivation = "splitmix64_block_sample_v1";

struct BlockKey {
  int n_up = 0;
  int n_down = 0;
  int mx = 0;
  int my = 0;

  bool operator==(const BlockKey& other) const {
    return n_up == other.n_up && n_down == other.n_down &&
           mx == other.mx && my == other.my;
  }

  bool operator<(const BlockKey& other) const {
    if (n_up != other.n_up) return n_up < other.n_up;
    if (n_down != other.n_down) return n_down < other.n_down;
    if (my != other.my) return my < other.my;
    return mx < other.mx;
  }
};

struct BlockKeyHash {
  std::size_t operator()(const BlockKey& key) const;
};

struct SectorKey {
  int n_up = 0;
  int n_down = 0;

  bool operator==(const SectorKey& other) const {
    return n_up == other.n_up && n_down == other.n_down;
  }
};

struct SectorKeyHash {
  std::size_t operator()(const SectorKey& key) const;
};

enum class Format {
  kMissing,
  kEmpty,
  kV1,
  kV2,
  kV3,
  kUnknown,
};

struct SampleRecord {
  BlockKey key;
  int particles = 0;
  std::uint64_t basis_dim = 0;
  int sample_id = -1;
  // Zero for legacy/v2 records.  V3 records tag every compact Ritz spectrum
  // by the Lanczos prefix length from which it was constructed.
  int lanczos_steps = 0;
  std::vector<double> eigenvalues;
  std::vector<double> overlaps;
  std::uint64_t checksum = 0;
};

struct ExactRecord {
  BlockKey key;
  int particles = 0;
  std::uint64_t basis_dim = 0;
  std::vector<double> eigenvalues;
  std::uint64_t checksum = 0;
};

// Legacy v1 checkpoints stored all samples for a block in one already-normalized
// record. They remain reducible only at the sample count in their metadata.
struct LegacyBlockRecord {
  BlockKey key;
  int particles = 0;
  std::uint64_t basis_dim = 0;
  double trace_prefactor = 1.0;
  bool exact = false;
  std::vector<double> eigenvalues;
  std::vector<double> weights;
};

struct CheckpointData {
  Format format = Format::kMissing;
  std::unordered_map<BlockKey, std::map<int, SampleRecord>, BlockKeyHash> samples;
  // V3 layout: block -> permanent sample ID -> Lanczos prefix -> spectrum.
  std::unordered_map<
      BlockKey,
      std::map<int, std::map<int, SampleRecord>>,
      BlockKeyHash> depth_samples;
  std::unordered_map<BlockKey, ExactRecord, BlockKeyHash> exact_blocks;
  std::unordered_map<BlockKey, int, BlockKeyHash> block_complete_samples;
  std::unordered_map<SectorKey, int, SectorKeyHash> sector_complete_samples;
  int run_complete_samples = 0;
  std::unordered_map<BlockKey, std::map<int, int>, BlockKeyHash>
      block_complete_depth_samples;
  std::unordered_map<SectorKey, std::map<int, int>, SectorKeyHash>
      sector_complete_depth_samples;
  std::map<int, int> run_complete_depth_samples;
  std::unordered_map<BlockKey, LegacyBlockRecord, BlockKeyHash> legacy_blocks;
  std::uint64_t valid_bytes = 0;
  bool trailing_partial_record = false;
};

struct Metadata {
  std::map<std::string, std::string> values;

  bool has(const std::string& key) const;
  std::string get(const std::string& key) const;
  int get_int(const std::string& key) const;
  std::uint64_t get_u64(const std::string& key) const;
  double get_double(const std::string& key) const;
};

struct ManifestEntry {
  BlockKey key;
  std::uint64_t basis_dim = 0;
  bool exact = false;
  int particles = 0;
  int krylov_steps = 0;
};

struct CompletionStatus {
  int target_samples = 0;
  int target_lanczos_steps = 0;
  int minimum_complete_samples = 0;
  std::uint64_t expected_sample_records = 0;
  std::uint64_t durable_sample_records = 0;
  std::uint64_t expected_blocks = 0;
  std::uint64_t completed_blocks = 0;
  long double total_weight = 0.0L;
  long double durable_weight = 0.0L;
  double weighted_fraction = 0.0;
  bool complete = false;
  bool has_next_missing = false;
  BlockKey next_block;
  int next_sample_id = -1;
  int next_lanczos_steps = 0;
};

struct ThermoPoint {
  double density = 0.0;
  double charge_correlation = 0.0;
  double compressibility = 0.0;
  double partition_like = 0.0;
  double log_partition = 0.0;
};

using ThermoGrid = std::vector<std::vector<ThermoPoint>>;

Format detect_format(const std::string& path);
Metadata read_metadata(const std::string& path);
void write_metadata_if_missing(const std::string& path, const std::string& text);
// Index-only mode retains record keys, IDs, checksums, and dimensions while
// dropping compact spectra.  Production jobs use it so checkpoint growth does
// not increase the resident spectrum payload; reducers request full payloads.
CheckpointData read_checkpoint(
    const std::string& path,
    bool load_spectra = true,
    int spectrum_lanczos_steps = 0);
void initialize_v2_checkpoint(const std::string& path);
void initialize_v3_checkpoint(const std::string& path);
void truncate_to_valid_bytes(const std::string& path, std::uint64_t valid_bytes);

class Writer {
 public:
  explicit Writer(std::string path);
  void append_sample(const SampleRecord& record);
  // Append all depth records for one random vector under one open/fsync.  Each
  // prefix remains an independently checksummed record on disk.
  void append_sample_bundle(const std::vector<SampleRecord>& records);
  void append_exact(const ExactRecord& record);
  void append_block_complete(const BlockKey& key, int samples, int lanczos_steps = 0);
  void append_sector_complete(const SectorKey& key, int samples, int lanczos_steps = 0);
  void append_run_complete(int samples, int lanczos_steps = 0);

 private:
  std::string path_;
  Format format_ = Format::kUnknown;
  std::mutex mutex_;
};

std::vector<ManifestEntry> build_manifest(
    int lx,
    int ly,
    std::size_t exact_threshold,
    int lanczos_steps);

int contiguous_samples(
    const CheckpointData& data,
    const BlockKey& key,
    int lanczos_steps = 0);

bool has_sample(
    const CheckpointData& data,
    const BlockKey& key,
    int sample_id,
    int lanczos_steps = 0);

bool block_complete(
    const CheckpointData& data,
    const ManifestEntry& entry,
    int target_samples,
    int lanczos_steps = 0);

CompletionStatus completion_status(
    const CheckpointData& data,
    const std::vector<ManifestEntry>& manifest,
    int target_samples,
    int lanczos_steps = 0);

ThermoGrid reduce_checkpoint(
    const CheckpointData& data,
    const std::vector<ManifestEntry>& manifest,
    int target_samples,
    int sites,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    int legacy_samples = 0,
    int lanczos_steps = 0);

void write_thermo_csv(
    const std::string& path,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    const ThermoGrid& grid);

std::string block_key_text(const BlockKey& key);
std::string status_text(const CompletionStatus& status);

}  // namespace ftlm_checkpoint
