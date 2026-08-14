#include "thermo_checkpoint.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <complex>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <type_traits>

#include <fcntl.h>
#include <unistd.h>

namespace ftlm_checkpoint {
namespace {

constexpr std::array<char, 8> kV1Magic = {'F', 'T', 'L', 'M', 'C', 'P', '1', '\n'};
constexpr std::array<char, 8> kV2Magic = {'F', 'T', 'L', 'M', 'C', 'P', '2', '\n'};
constexpr std::array<char, 8> kV3Magic = {'F', 'T', 'L', 'M', 'C', 'P', '3', '\n'};
constexpr std::array<char, 8> kV1Record = {'B', 'L', 'K', 'R', 'E', 'C', '1', '\n'};
constexpr std::array<char, 8> kV2Record = {'F', '2', 'R', 'E', 'C', '1', '\n', '\0'};
constexpr std::array<char, 8> kV3Record = {'F', '3', 'R', 'E', 'C', '1', '\n', '\0'};

enum class RecordKind : std::uint32_t {
  kSample = 1,
  kExact = 2,
  kBlockComplete = 3,
  kSectorComplete = 4,
  kRunComplete = 5,
};

[[noreturn]] void fail(const std::string& message) {
  throw std::runtime_error(message);
}

std::uint64_t fnv1a64(const std::vector<std::uint8_t>& bytes) {
  std::uint64_t hash = 1469598103934665603ULL;
  for (std::uint8_t byte : bytes) {
    hash ^= static_cast<std::uint64_t>(byte);
    hash *= 1099511628211ULL;
  }
  return hash;
}

template <typename T>
void append_pod(std::vector<std::uint8_t>& out, const T& value) {
  static_assert(std::is_trivially_copyable<T>::value, "POD required");
  const auto* ptr = reinterpret_cast<const std::uint8_t*>(&value);
  out.insert(out.end(), ptr, ptr + sizeof(T));
}

template <typename T>
T read_pod(const std::vector<std::uint8_t>& bytes, std::size_t& cursor) {
  static_assert(std::is_trivially_copyable<T>::value, "POD required");
  if (cursor + sizeof(T) > bytes.size()) {
    fail("Checkpoint payload is truncated.");
  }
  T value{};
  std::memcpy(&value, bytes.data() + cursor, sizeof(T));
  cursor += sizeof(T);
  return value;
}

void append_key(std::vector<std::uint8_t>& out, const BlockKey& key) {
  append_pod(out, static_cast<std::int32_t>(key.n_up));
  append_pod(out, static_cast<std::int32_t>(key.n_down));
  append_pod(out, static_cast<std::int32_t>(key.mx));
  append_pod(out, static_cast<std::int32_t>(key.my));
}

BlockKey read_key(const std::vector<std::uint8_t>& bytes, std::size_t& cursor) {
  BlockKey key;
  key.n_up = read_pod<std::int32_t>(bytes, cursor);
  key.n_down = read_pod<std::int32_t>(bytes, cursor);
  key.mx = read_pod<std::int32_t>(bytes, cursor);
  key.my = read_pod<std::int32_t>(bytes, cursor);
  return key;
}

void append_double_vector(std::vector<std::uint8_t>& out, const std::vector<double>& values) {
  append_pod(out, static_cast<std::uint64_t>(values.size()));
  if (!values.empty()) {
    const auto* ptr = reinterpret_cast<const std::uint8_t*>(values.data());
    out.insert(out.end(), ptr, ptr + values.size() * sizeof(double));
  }
}

std::vector<double> read_double_vector(
    const std::vector<std::uint8_t>& bytes,
    std::size_t& cursor) {
  const std::uint64_t count = read_pod<std::uint64_t>(bytes, cursor);
  if (count > 100000000ULL || count > (bytes.size() - cursor) / sizeof(double)) {
    fail("Checkpoint vector length is implausible.");
  }
  std::vector<double> values(static_cast<std::size_t>(count));
  const std::size_t byte_count = values.size() * sizeof(double);
  if (byte_count > 0) {
    std::memcpy(values.data(), bytes.data() + cursor, byte_count);
    cursor += byte_count;
  }
  return values;
}

void write_all(int fd, const void* data, std::size_t size, const std::string& path) {
  const auto* ptr = static_cast<const std::uint8_t*>(data);
  std::size_t written = 0;
  while (written < size) {
    const ssize_t result = ::write(fd, ptr + written, size - written);
    if (result < 0) {
      if (errno == EINTR) continue;
      fail("Failed writing checkpoint " + path + ": " + std::strerror(errno));
    }
    written += static_cast<std::size_t>(result);
  }
}

void write_record_frame(
    int fd,
    const std::string& path,
    Format format,
    RecordKind kind,
    const std::vector<std::uint8_t>& payload) {
  const std::uint32_t kind_value = static_cast<std::uint32_t>(kind);
  const std::uint64_t payload_size = static_cast<std::uint64_t>(payload.size());
  const std::uint64_t checksum = fnv1a64(payload);
  const auto& record_magic = format == Format::kV3 ? kV3Record : kV2Record;
  write_all(fd, record_magic.data(), record_magic.size(), path);
  write_all(fd, &kind_value, sizeof(kind_value), path);
  write_all(fd, &payload_size, sizeof(payload_size), path);
  write_all(fd, &checksum, sizeof(checksum), path);
  if (!payload.empty()) write_all(fd, payload.data(), payload.size(), path);
}

void append_records(
    const std::string& path,
    Format format,
    const std::vector<std::pair<RecordKind, std::vector<std::uint8_t>>>& records) {
  if (records.empty()) return;
  const int fd = ::open(path.c_str(), O_WRONLY | O_APPEND);
  if (fd < 0) {
    fail("Failed to append checkpoint " + path + ": " + std::strerror(errno));
  }
  try {
    for (const auto& record : records) {
      write_record_frame(fd, path, format, record.first, record.second);
    }
    if (::fsync(fd) != 0) {
      fail("Failed to fsync checkpoint " + path + ": " + std::strerror(errno));
    }
  } catch (...) {
    ::close(fd);
    throw;
  }
  if (::close(fd) != 0) {
    fail("Failed to close checkpoint " + path + ": " + std::strerror(errno));
  }
}

void append_record(
    const std::string& path,
    Format format,
    RecordKind kind,
    const std::vector<std::uint8_t>& payload) {
  append_records(path, format, {{kind, payload}});
}

template <typename T>
bool read_binary(std::istream& in, T& value) {
  return static_cast<bool>(in.read(reinterpret_cast<char*>(&value), sizeof(T)));
}

double log_add(double lhs, double rhs) {
  if (!std::isfinite(lhs)) return rhs;
  if (!std::isfinite(rhs)) return lhs;
  const double hi = std::max(lhs, rhs);
  const double lo = std::min(lhs, rhs);
  return hi + std::log1p(std::exp(lo - hi));
}

std::vector<int> permutation_cycle_lengths(int lx, int ly, int dx, int dy) {
  const int sites = lx * ly;
  std::vector<bool> seen(static_cast<std::size_t>(sites), false);
  std::vector<int> lengths;
  for (int site = 0; site < sites; ++site) {
    if (seen[static_cast<std::size_t>(site)]) continue;
    int current = site;
    int length = 0;
    while (!seen[static_cast<std::size_t>(current)]) {
      seen[static_cast<std::size_t>(current)] = true;
      ++length;
      const int x = current % lx;
      const int y = current / lx;
      current = ((y + dy) % ly) * lx + ((x + dx) % lx);
    }
    lengths.push_back(length);
  }
  return lengths;
}

long long fermion_translation_trace(const std::vector<int>& cycles, int particles) {
  std::vector<long long> polynomial(1, 1);
  for (int length : cycles) {
    std::vector<long long> next(polynomial.size() + static_cast<std::size_t>(length), 0);
    const long long coefficient = (length % 2 == 0) ? -1LL : 1LL;
    for (std::size_t i = 0; i < polynomial.size(); ++i) {
      next[i] += polynomial[i];
      next[i + static_cast<std::size_t>(length)] += coefficient * polynomial[i];
    }
    polynomial.swap(next);
  }
  if (particles < 0 || static_cast<std::size_t>(particles) >= polynomial.size()) return 0;
  return polynomial[static_cast<std::size_t>(particles)];
}

std::string trim(const std::string& text) {
  const std::size_t first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return {};
  const std::size_t last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

}  // namespace

std::size_t BlockKeyHash::operator()(const BlockKey& key) const {
  std::size_t value = std::hash<int>{}(key.n_up);
  value ^= std::hash<int>{}(key.n_down) + 0x9e3779b9U + (value << 6U) + (value >> 2U);
  value ^= std::hash<int>{}(key.mx) + 0x9e3779b9U + (value << 6U) + (value >> 2U);
  value ^= std::hash<int>{}(key.my) + 0x9e3779b9U + (value << 6U) + (value >> 2U);
  return value;
}

std::size_t SectorKeyHash::operator()(const SectorKey& key) const {
  std::size_t value = std::hash<int>{}(key.n_up);
  value ^= std::hash<int>{}(key.n_down) + 0x9e3779b9U + (value << 6U) + (value >> 2U);
  return value;
}

bool Metadata::has(const std::string& key) const {
  return values.find(key) != values.end();
}

std::string Metadata::get(const std::string& key) const {
  const auto it = values.find(key);
  if (it == values.end()) fail("Checkpoint metadata is missing key: " + key);
  return it->second;
}

int Metadata::get_int(const std::string& key) const {
  return std::stoi(get(key));
}

std::uint64_t Metadata::get_u64(const std::string& key) const {
  return static_cast<std::uint64_t>(std::stoull(get(key)));
}

double Metadata::get_double(const std::string& key) const {
  return std::stod(get(key));
}

Format detect_format(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) return Format::kMissing;
  std::array<char, 8> magic{};
  if (!in.read(magic.data(), static_cast<std::streamsize>(magic.size()))) {
    return Format::kEmpty;
  }
  if (magic == kV1Magic) return Format::kV1;
  if (magic == kV2Magic) return Format::kV2;
  if (magic == kV3Magic) return Format::kV3;
  return Format::kUnknown;
}

Metadata read_metadata(const std::string& path) {
  std::ifstream in(path);
  if (!in) fail("Failed to open checkpoint metadata: " + path);
  Metadata metadata;
  std::string line;
  while (std::getline(in, line)) {
    if (trim(line).empty()) continue;
    const std::size_t equals = line.find('=');
    if (equals == std::string::npos) fail("Malformed checkpoint metadata line: " + line);
    metadata.values[trim(line.substr(0, equals))] = trim(line.substr(equals + 1));
  }
  return metadata;
}

void write_metadata_if_missing(const std::string& path, const std::string& text) {
  if (std::filesystem::exists(path)) return;
  std::ofstream out(path);
  if (!out) fail("Failed to write checkpoint metadata: " + path);
  out << text;
}

void initialize_v2_checkpoint(const std::string& path) {
  const Format format = detect_format(path);
  if (format == Format::kV2) return;
  if (format != Format::kMissing && format != Format::kEmpty) {
    fail("Cannot initialize v2 checkpoint over an existing non-v2 file: " + path);
  }
  const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0664);
  if (fd < 0) fail("Failed to initialize checkpoint " + path + ": " + std::strerror(errno));
  try {
    write_all(fd, kV2Magic.data(), kV2Magic.size(), path);
    if (::fsync(fd) != 0) fail("Failed to fsync new checkpoint: " + path);
  } catch (...) {
    ::close(fd);
    throw;
  }
  ::close(fd);
}

void initialize_v3_checkpoint(const std::string& path) {
  const Format format = detect_format(path);
  if (format == Format::kV3) return;
  if (format != Format::kMissing && format != Format::kEmpty) {
    fail("Cannot initialize v3 checkpoint over an existing non-v3 file: " + path);
  }
  const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0664);
  if (fd < 0) fail("Failed to initialize checkpoint " + path + ": " + std::strerror(errno));
  try {
    write_all(fd, kV3Magic.data(), kV3Magic.size(), path);
    if (::fsync(fd) != 0) fail("Failed to fsync new checkpoint: " + path);
  } catch (...) {
    ::close(fd);
    throw;
  }
  ::close(fd);
}

CheckpointData read_checkpoint(
    const std::string& path,
    bool load_spectra,
    int spectrum_lanczos_steps) {
  CheckpointData data;
  data.format = detect_format(path);
  if (data.format == Format::kMissing || data.format == Format::kEmpty) return data;
  if (data.format == Format::kUnknown) fail("Unrecognized checkpoint format: " + path);

  std::ifstream in(path, std::ios::binary);
  if (!in) fail("Failed to read checkpoint: " + path);
  std::array<char, 8> magic{};
  in.read(magic.data(), static_cast<std::streamsize>(magic.size()));
  data.valid_bytes = magic.size();

  if (data.format == Format::kV1) {
    while (true) {
      std::array<char, 8> record_magic{};
      if (!in.read(record_magic.data(), static_cast<std::streamsize>(record_magic.size()))) break;
      if (record_magic != kV1Record) break;
      LegacyBlockRecord record;
      std::uint8_t exact_flag = 0;
      std::uint64_t eigen_count = 0;
      std::uint64_t weight_count = 0;
      if (!read_binary(in, record.key.n_up) || !read_binary(in, record.key.n_down) ||
          !read_binary(in, record.key.mx) || !read_binary(in, record.key.my) ||
          !read_binary(in, record.particles) || !read_binary(in, record.basis_dim) ||
          !read_binary(in, record.trace_prefactor) || !read_binary(in, exact_flag) ||
          !read_binary(in, eigen_count) || !read_binary(in, weight_count)) {
        break;
      }
      if (eigen_count > 100000000ULL || weight_count > 100000000ULL) {
        fail("Implausibly large v1 checkpoint record: " + path);
      }
      record.exact = exact_flag != 0;
      record.eigenvalues.resize(static_cast<std::size_t>(eigen_count));
      record.weights.resize(static_cast<std::size_t>(weight_count));
      if (!in.read(reinterpret_cast<char*>(record.eigenvalues.data()),
                   static_cast<std::streamsize>(record.eigenvalues.size() * sizeof(double))) ||
          !in.read(reinterpret_cast<char*>(record.weights.data()),
                   static_cast<std::streamsize>(record.weights.size() * sizeof(double)))) {
        break;
      }
      data.legacy_blocks[record.key] = std::move(record);
      data.valid_bytes = static_cast<std::uint64_t>(in.tellg());
    }
    return data;
  }

  while (true) {
    const std::streampos record_start = in.tellg();
    std::array<char, 8> record_magic{};
    if (!in.read(record_magic.data(), static_cast<std::streamsize>(record_magic.size()))) {
      if (in.gcount() != 0 || !in.eof()) data.trailing_partial_record = true;
      break;
    }
    const auto& expected_record_magic =
        data.format == Format::kV3 ? kV3Record : kV2Record;
    if (record_magic != expected_record_magic) {
      fail("Corrupted checkpoint record magic at byte " +
           std::to_string(static_cast<long long>(record_start)) + ": " + path);
    }
    std::uint32_t kind_value = 0;
    std::uint64_t payload_size = 0;
    std::uint64_t expected_checksum = 0;
    if (!read_binary(in, kind_value) || !read_binary(in, payload_size) ||
        !read_binary(in, expected_checksum)) {
      data.trailing_partial_record = true;
      break;
    }
    if (payload_size > (1ULL << 32U)) fail("Implausibly large checkpoint payload.");
    std::vector<std::uint8_t> payload(static_cast<std::size_t>(payload_size));
    if (!payload.empty() &&
        !in.read(reinterpret_cast<char*>(payload.data()),
                 static_cast<std::streamsize>(payload.size()))) {
      data.trailing_partial_record = true;
      break;
    }
    const std::uint64_t checksum = fnv1a64(payload);
    if (checksum != expected_checksum) {
      fail("Checksum mismatch in complete checkpoint record: " + path);
    }

    std::size_t cursor = 0;
    const RecordKind kind = static_cast<RecordKind>(kind_value);
    if (kind == RecordKind::kSample) {
      SampleRecord record;
      record.key = read_key(payload, cursor);
      record.particles = read_pod<std::int32_t>(payload, cursor);
      record.basis_dim = read_pod<std::uint64_t>(payload, cursor);
      record.sample_id = read_pod<std::int32_t>(payload, cursor);
      if (data.format == Format::kV3) {
        record.lanczos_steps = read_pod<std::int32_t>(payload, cursor);
      }
      record.eigenvalues = read_double_vector(payload, cursor);
      record.overlaps = read_double_vector(payload, cursor);
      record.checksum = checksum;
      if (record.sample_id < 0 ||
          (data.format == Format::kV3 && record.lanczos_steps <= 0) ||
          record.eigenvalues.size() != record.overlaps.size()) {
        fail("Invalid checkpoint sample record.");
      }
      if (data.format == Format::kV3) {
        auto& by_depth = data.depth_samples[record.key][record.sample_id];
        const auto existing = by_depth.find(record.lanczos_steps);
        if (existing != by_depth.end() && existing->second.checksum != checksum) {
          fail("Conflicting duplicate sample record for " + block_key_text(record.key) +
               " sample=" + std::to_string(record.sample_id) +
               " m=" + std::to_string(record.lanczos_steps));
        }
        if (!load_spectra ||
            (spectrum_lanczos_steps > 0 &&
             record.lanczos_steps != spectrum_lanczos_steps)) {
          record.eigenvalues.clear();
          record.overlaps.clear();
        }
        by_depth[record.lanczos_steps] = std::move(record);
      } else {
        auto& by_id = data.samples[record.key];
        const auto existing = by_id.find(record.sample_id);
        if (existing != by_id.end() && existing->second.checksum != checksum) {
          fail("Conflicting duplicate sample record for " + block_key_text(record.key) +
               " sample=" + std::to_string(record.sample_id));
        }
        if (!load_spectra) {
          record.eigenvalues.clear();
          record.overlaps.clear();
        }
        by_id[record.sample_id] = std::move(record);
      }
    } else if (kind == RecordKind::kExact) {
      ExactRecord record;
      record.key = read_key(payload, cursor);
      record.particles = read_pod<std::int32_t>(payload, cursor);
      record.basis_dim = read_pod<std::uint64_t>(payload, cursor);
      record.eigenvalues = read_double_vector(payload, cursor);
      record.checksum = checksum;
      const auto existing = data.exact_blocks.find(record.key);
      if (existing != data.exact_blocks.end() && existing->second.checksum != checksum) {
        fail("Conflicting duplicate exact record for " + block_key_text(record.key));
      }
      if (!load_spectra) record.eigenvalues.clear();
      data.exact_blocks[record.key] = std::move(record);
    } else if (kind == RecordKind::kBlockComplete) {
      const BlockKey key = read_key(payload, cursor);
      const int lanczos_steps = data.format == Format::kV3
          ? read_pod<std::int32_t>(payload, cursor) : 0;
      const int samples = read_pod<std::int32_t>(payload, cursor);
      if (data.format == Format::kV3) {
        auto& achieved = data.block_complete_depth_samples[key][lanczos_steps];
        achieved = std::max(achieved, samples);
      } else {
        data.block_complete_samples[key] =
            std::max(data.block_complete_samples[key], samples);
      }
    } else if (kind == RecordKind::kSectorComplete) {
      SectorKey key;
      key.n_up = read_pod<std::int32_t>(payload, cursor);
      key.n_down = read_pod<std::int32_t>(payload, cursor);
      const int lanczos_steps = data.format == Format::kV3
          ? read_pod<std::int32_t>(payload, cursor) : 0;
      const int samples = read_pod<std::int32_t>(payload, cursor);
      if (data.format == Format::kV3) {
        auto& achieved = data.sector_complete_depth_samples[key][lanczos_steps];
        achieved = std::max(achieved, samples);
      } else {
        data.sector_complete_samples[key] =
            std::max(data.sector_complete_samples[key], samples);
      }
    } else if (kind == RecordKind::kRunComplete) {
      if (data.format == Format::kV3) {
        const int lanczos_steps = read_pod<std::int32_t>(payload, cursor);
        const int samples = read_pod<std::int32_t>(payload, cursor);
        data.run_complete_depth_samples[lanczos_steps] =
            std::max(data.run_complete_depth_samples[lanczos_steps], samples);
      } else {
        data.run_complete_samples =
            std::max(data.run_complete_samples, read_pod<std::int32_t>(payload, cursor));
      }
    } else {
      fail("Unknown checkpoint record kind: " + std::to_string(kind_value));
    }
    if (cursor != payload.size()) fail("Unexpected trailing bytes in checkpoint record.");
    data.valid_bytes = static_cast<std::uint64_t>(in.tellg());
  }
  return data;
}

void truncate_to_valid_bytes(const std::string& path, std::uint64_t valid_bytes) {
  if (::truncate(path.c_str(), static_cast<off_t>(valid_bytes)) != 0) {
    fail("Failed to truncate incomplete checkpoint tail " + path + ": " +
         std::strerror(errno));
  }
}

Writer::Writer(std::string path) : path_(std::move(path)) {
  format_ = detect_format(path_);
  if (format_ == Format::kMissing || format_ == Format::kEmpty) {
    initialize_v2_checkpoint(path_);
    format_ = Format::kV2;
  }
  if (format_ != Format::kV2 && format_ != Format::kV3) {
    fail("Writer requires a v2 or v3 checkpoint: " + path_);
  }
}

void Writer::append_sample(const SampleRecord& record) {
  append_sample_bundle({record});
}

void Writer::append_sample_bundle(const std::vector<SampleRecord>& records) {
  std::vector<std::pair<RecordKind, std::vector<std::uint8_t>>> encoded;
  encoded.reserve(records.size());
  for (const SampleRecord& record : records) {
    if (record.sample_id < 0 || record.eigenvalues.size() != record.overlaps.size()) {
      fail("Invalid sample record passed to checkpoint writer.");
    }
    std::vector<std::uint8_t> payload;
    append_key(payload, record.key);
    append_pod(payload, static_cast<std::int32_t>(record.particles));
    append_pod(payload, record.basis_dim);
    append_pod(payload, static_cast<std::int32_t>(record.sample_id));
    if (format_ == Format::kV3) {
      if (record.lanczos_steps <= 0) {
        fail("V3 sample record requires positive Lanczos steps.");
      }
      append_pod(payload, static_cast<std::int32_t>(record.lanczos_steps));
    }
    append_double_vector(payload, record.eigenvalues);
    append_double_vector(payload, record.overlaps);
    encoded.emplace_back(RecordKind::kSample, std::move(payload));
  }
  std::lock_guard<std::mutex> lock(mutex_);
  append_records(path_, format_, encoded);
}

void Writer::append_exact(const ExactRecord& record) {
  std::vector<std::uint8_t> payload;
  append_key(payload, record.key);
  append_pod(payload, static_cast<std::int32_t>(record.particles));
  append_pod(payload, record.basis_dim);
  append_double_vector(payload, record.eigenvalues);
  std::lock_guard<std::mutex> lock(mutex_);
  append_record(path_, format_, RecordKind::kExact, payload);
}

void Writer::append_block_complete(const BlockKey& key, int samples, int lanczos_steps) {
  std::vector<std::uint8_t> payload;
  append_key(payload, key);
  if (format_ == Format::kV3) {
    if (lanczos_steps <= 0) fail("V3 block marker requires positive Lanczos steps.");
    append_pod(payload, static_cast<std::int32_t>(lanczos_steps));
  }
  append_pod(payload, static_cast<std::int32_t>(samples));
  std::lock_guard<std::mutex> lock(mutex_);
  append_record(path_, format_, RecordKind::kBlockComplete, payload);
}

void Writer::append_sector_complete(const SectorKey& key, int samples, int lanczos_steps) {
  std::vector<std::uint8_t> payload;
  append_pod(payload, static_cast<std::int32_t>(key.n_up));
  append_pod(payload, static_cast<std::int32_t>(key.n_down));
  if (format_ == Format::kV3) {
    if (lanczos_steps <= 0) fail("V3 sector marker requires positive Lanczos steps.");
    append_pod(payload, static_cast<std::int32_t>(lanczos_steps));
  }
  append_pod(payload, static_cast<std::int32_t>(samples));
  std::lock_guard<std::mutex> lock(mutex_);
  append_record(path_, format_, RecordKind::kSectorComplete, payload);
}

void Writer::append_run_complete(int samples, int lanczos_steps) {
  std::vector<std::uint8_t> payload;
  if (format_ == Format::kV3) {
    if (lanczos_steps <= 0) fail("V3 run marker requires positive Lanczos steps.");
    append_pod(payload, static_cast<std::int32_t>(lanczos_steps));
  }
  append_pod(payload, static_cast<std::int32_t>(samples));
  std::lock_guard<std::mutex> lock(mutex_);
  append_record(path_, format_, RecordKind::kRunComplete, payload);
}

std::vector<ManifestEntry> build_manifest(
    int lx,
    int ly,
    std::size_t exact_threshold,
    int lanczos_steps) {
  if (lx <= 0 || ly <= 0 || lanczos_steps <= 0) fail("Invalid manifest parameters.");
  const int sites = lx * ly;
  const double pi = std::acos(-1.0);
  std::vector<std::vector<long long>> trace_by_shift(
      static_cast<std::size_t>(lx * ly),
      std::vector<long long>(static_cast<std::size_t>(sites + 1), 0));
  for (int dy = 0; dy < ly; ++dy) {
    for (int dx = 0; dx < lx; ++dx) {
      const auto cycles = permutation_cycle_lengths(lx, ly, dx, dy);
      auto& trace = trace_by_shift[static_cast<std::size_t>(dy * lx + dx)];
      for (int particles = 0; particles <= sites; ++particles) {
        trace[static_cast<std::size_t>(particles)] =
            fermion_translation_trace(cycles, particles);
      }
    }
  }

  std::vector<ManifestEntry> manifest;
  for (int n_up = 0; n_up <= sites; ++n_up) {
    for (int n_down = 0; n_down <= sites; ++n_down) {
      for (int my = 0; my < ly; ++my) {
        for (int mx = 0; mx < lx; ++mx) {
          std::complex<long double> projection = 0.0L;
          for (int dy = 0; dy < ly; ++dy) {
            for (int dx = 0; dx < lx; ++dx) {
              const double angle = -2.0 * pi *
                  (static_cast<double>(mx * dx) / static_cast<double>(lx) +
                   static_cast<double>(my * dy) / static_cast<double>(ly));
              const std::complex<long double> phase(std::cos(angle), std::sin(angle));
              const auto& trace = trace_by_shift[static_cast<std::size_t>(dy * lx + dx)];
              projection += phase * static_cast<long double>(trace[static_cast<std::size_t>(n_up)]) *
                            static_cast<long double>(trace[static_cast<std::size_t>(n_down)]);
            }
          }
          projection /= static_cast<long double>(sites);
          const long long dimension = std::llround(projection.real());
          if (dimension < 0 || std::abs(projection.imag()) > 1e-5L ||
              std::abs(projection.real() - static_cast<long double>(dimension)) > 1e-5L) {
            fail("Failed to construct an integer momentum-block manifest.");
          }
          if (dimension == 0) continue;
          ManifestEntry entry;
          entry.key = {n_up, n_down, mx, my};
          entry.basis_dim = static_cast<std::uint64_t>(dimension);
          entry.exact = entry.basis_dim <= exact_threshold;
          entry.particles = n_up + n_down;
          entry.krylov_steps = std::min<int>(lanczos_steps, static_cast<int>(entry.basis_dim));
          manifest.push_back(entry);
        }
      }
    }
  }
  return manifest;
}

bool has_sample(
    const CheckpointData& data,
    const BlockKey& key,
    int sample_id,
    int lanczos_steps) {
  if (data.format == Format::kV3) {
    if (lanczos_steps <= 0) fail("V3 sample lookup requires positive Lanczos steps.");
    const auto block = data.depth_samples.find(key);
    if (block == data.depth_samples.end()) return false;
    const auto sample = block->second.find(sample_id);
    return sample != block->second.end() &&
        sample->second.find(lanczos_steps) != sample->second.end();
  }
  const auto block = data.samples.find(key);
  return block != data.samples.end() && block->second.find(sample_id) != block->second.end();
}

int contiguous_samples(
    const CheckpointData& data,
    const BlockKey& key,
    int lanczos_steps) {
  int count = 0;
  while (has_sample(data, key, count, lanczos_steps)) ++count;
  return count;
}

bool block_complete(
    const CheckpointData& data,
    const ManifestEntry& entry,
    int target_samples,
    int lanczos_steps) {
  if (data.format == Format::kV1) {
    return data.legacy_blocks.find(entry.key) != data.legacy_blocks.end();
  }
  if (entry.exact) return data.exact_blocks.find(entry.key) != data.exact_blocks.end();
  return contiguous_samples(data, entry.key, lanczos_steps) >= target_samples;
}

CompletionStatus completion_status(
    const CheckpointData& data,
    const std::vector<ManifestEntry>& manifest,
    int target_samples,
    int lanczos_steps) {
  if (target_samples <= 0) fail("Target sample count must be positive.");
  if (data.format == Format::kV3 && lanczos_steps <= 0) {
    fail("V3 completion status requires positive Lanczos steps.");
  }
  CompletionStatus status;
  status.target_samples = target_samples;
  status.target_lanczos_steps = lanczos_steps;
  status.expected_blocks = manifest.size();
  status.minimum_complete_samples = std::numeric_limits<int>::max();
  for (const ManifestEntry& entry : manifest) {
    if (entry.exact) {
      if (block_complete(data, entry, target_samples, lanczos_steps)) {
        ++status.completed_blocks;
      } else if (!status.has_next_missing) {
        status.has_next_missing = true;
        status.next_block = entry.key;
        status.next_sample_id = -1;
      }
      continue;
    }
    status.expected_sample_records += static_cast<std::uint64_t>(target_samples);
    const long double sample_weight =
        static_cast<long double>(entry.basis_dim) * entry.krylov_steps;
    status.total_weight += sample_weight * target_samples;
    int contiguous = 0;
    if (data.format == Format::kV1) {
      contiguous = data.legacy_blocks.find(entry.key) != data.legacy_blocks.end()
          ? target_samples : 0;
    } else {
      contiguous = contiguous_samples(data, entry.key, lanczos_steps);
    }
    status.minimum_complete_samples = std::min(status.minimum_complete_samples, contiguous);
    int present_for_target = 0;
    int first_missing = -1;
    if (data.format == Format::kV1) {
      present_for_target = contiguous;
    } else {
      for (int sample_id = 0; sample_id < target_samples; ++sample_id) {
        const bool present = has_sample(data, entry.key, sample_id, lanczos_steps);
        if (present) {
          ++present_for_target;
        } else if (first_missing < 0) {
          first_missing = sample_id;
        }
      }
    }
    status.durable_sample_records += static_cast<std::uint64_t>(present_for_target);
    status.durable_weight += sample_weight * present_for_target;
    if (present_for_target == target_samples) {
      ++status.completed_blocks;
    } else if (!status.has_next_missing) {
      status.has_next_missing = true;
      status.next_block = entry.key;
      status.next_sample_id = first_missing >= 0 ? first_missing : contiguous;
      status.next_lanczos_steps = lanczos_steps;
    }
  }
  if (status.minimum_complete_samples == std::numeric_limits<int>::max()) {
    status.minimum_complete_samples = target_samples;
  }
  status.weighted_fraction = status.total_weight > 0.0L
      ? static_cast<double>(status.durable_weight / status.total_weight)
      : (status.completed_blocks == status.expected_blocks ? 1.0 : 0.0);
  status.complete = status.completed_blocks == status.expected_blocks;
  return status;
}

ThermoGrid reduce_checkpoint(
    const CheckpointData& data,
    const std::vector<ManifestEntry>& manifest,
    int target_samples,
    int sites,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    int legacy_samples,
    int lanczos_steps) {
  if (sites <= 0 || beta_values.empty() || mu_values.empty()) fail("Invalid reducer grid.");
  if (data.format == Format::kV1 && target_samples != legacy_samples) {
    fail("A v1 checkpoint can only be reduced at its original sample count " +
         std::to_string(legacy_samples) + ".");
  }
  const CompletionStatus status =
      completion_status(data, manifest, target_samples, lanczos_steps);
  if (!status.complete) {
    std::string next = "; next missing " + block_key_text(status.next_block);
    next += status.next_sample_id >= 0
        ? " sample=" + std::to_string(status.next_sample_id)
        : " exact record";
    if (data.format == Format::kV3) {
      next += " m=" + std::to_string(lanczos_steps);
    }
    fail("Checkpoint is incomplete for R=" + std::to_string(target_samples) + next);
  }

  ThermoGrid grid(
      beta_values.size(),
      std::vector<ThermoPoint>(mu_values.size()));
  const int max_particles = 2 * sites;
  for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
    const double beta = beta_values[ibeta];
    std::vector<double> log_q(static_cast<std::size_t>(max_particles + 1),
                              -std::numeric_limits<double>::infinity());
    std::vector<double> min_energy(static_cast<std::size_t>(max_particles + 1),
                                   std::numeric_limits<double>::infinity());

    if (data.format == Format::kV1) {
      for (const auto& item : data.legacy_blocks) {
        const LegacyBlockRecord& record = item.second;
        const int particles = record.particles;
        for (std::size_t i = 0; i < record.eigenvalues.size(); ++i) {
          double coefficient = record.trace_prefactor;
          if (!record.exact) coefficient *= record.weights.at(i);
          if (!(coefficient > 0.0)) continue;
          const double term = std::log(coefficient) - beta * record.eigenvalues[i];
          log_q[static_cast<std::size_t>(particles)] =
              log_add(log_q[static_cast<std::size_t>(particles)], term);
          min_energy[static_cast<std::size_t>(particles)] =
              std::min(min_energy[static_cast<std::size_t>(particles)], record.eigenvalues[i]);
        }
      }
    } else {
      // Traverse exact blocks in manifest order.  Iterating the unordered map
      // directly makes the last few floating-point bits depend on its bucket
      // layout, which would violate exact lower-R reproducibility after a
      // checkpoint is extended and reloaded.
      for (const ManifestEntry& entry : manifest) {
        if (!entry.exact) continue;
        const auto exact = data.exact_blocks.find(entry.key);
        if (exact == data.exact_blocks.end()) {
          fail("Missing exact block during reduction.");
        }
        const ExactRecord& record = exact->second;
        for (double energy : record.eigenvalues) {
          const int particles = record.particles;
          log_q[static_cast<std::size_t>(particles)] =
              log_add(log_q[static_cast<std::size_t>(particles)], -beta * energy);
          min_energy[static_cast<std::size_t>(particles)] =
              std::min(min_energy[static_cast<std::size_t>(particles)], energy);
        }
      }
      for (const ManifestEntry& entry : manifest) {
        if (entry.exact) continue;
        const double log_prefactor =
            std::log(static_cast<double>(entry.basis_dim)) - std::log(target_samples);
        for (int sample_id = 0; sample_id < target_samples; ++sample_id) {
          const SampleRecord* record_ptr = nullptr;
          if (data.format == Format::kV3) {
            const auto block = data.depth_samples.find(entry.key);
            if (block != data.depth_samples.end()) {
              const auto sample = block->second.find(sample_id);
              if (sample != block->second.end()) {
                const auto depth = sample->second.find(lanczos_steps);
                if (depth != sample->second.end()) record_ptr = &depth->second;
              }
            }
          } else {
            const auto block = data.samples.find(entry.key);
            if (block != data.samples.end()) {
              const auto sample = block->second.find(sample_id);
              if (sample != block->second.end()) record_ptr = &sample->second;
            }
          }
          if (record_ptr == nullptr) fail("Missing sample record during reduction.");
          const SampleRecord& record = *record_ptr;
          for (std::size_t i = 0; i < record.eigenvalues.size(); ++i) {
            const double overlap = record.overlaps[i];
            if (!(overlap > 0.0)) continue;
            const double term = log_prefactor + std::log(overlap) -
                                beta * record.eigenvalues[i];
            log_q[static_cast<std::size_t>(entry.particles)] =
                log_add(log_q[static_cast<std::size_t>(entry.particles)], term);
            min_energy[static_cast<std::size_t>(entry.particles)] =
                std::min(min_energy[static_cast<std::size_t>(entry.particles)],
                         record.eigenvalues[i]);
          }
        }
      }
    }

    for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
      const double mu = mu_values[imu];
      double max_log_weight = -std::numeric_limits<double>::infinity();
      double min_grand_energy = std::numeric_limits<double>::infinity();
      for (int particles = 0; particles <= max_particles; ++particles) {
        if (std::isfinite(log_q[static_cast<std::size_t>(particles)])) {
          max_log_weight = std::max(
              max_log_weight,
              log_q[static_cast<std::size_t>(particles)] + beta * mu * particles);
        }
        if (std::isfinite(min_energy[static_cast<std::size_t>(particles)])) {
          min_grand_energy = std::min(
              min_grand_energy,
              min_energy[static_cast<std::size_t>(particles)] - mu * particles);
        }
      }
      if (!std::isfinite(max_log_weight)) fail("Reducer found an empty partition function.");
      long double z_scaled = 0.0L;
      long double n_scaled = 0.0L;
      long double n2_scaled = 0.0L;
      for (int particles = 0; particles <= max_particles; ++particles) {
        const double log_weight = log_q[static_cast<std::size_t>(particles)] +
                                  beta * mu * particles;
        if (!std::isfinite(log_weight)) continue;
        const long double weight = std::exp(log_weight - max_log_weight);
        z_scaled += weight;
        n_scaled += weight * particles;
        n2_scaled += weight * particles * particles;
      }
      const double log_z = max_log_weight + std::log(static_cast<double>(z_scaled));
      const double mean_n = static_cast<double>(n_scaled / z_scaled);
      const double mean_n2 = static_cast<double>(n2_scaled / z_scaled);
      const double charge = std::max(
          0.0, (mean_n2 - mean_n * mean_n) / static_cast<double>(sites));
      ThermoPoint& point = grid[ibeta][imu];
      point.density = mean_n / static_cast<double>(sites);
      point.charge_correlation = charge;
      point.compressibility = beta * charge;
      point.log_partition = log_z;
      const double shifted_log = log_z + beta * min_grand_energy;
      point.partition_like = shifted_log < 700.0
          ? std::exp(shifted_log) : std::numeric_limits<double>::infinity();
    }
  }
  return grid;
}

void write_thermo_csv(
    const std::string& path,
    const std::vector<double>& beta_values,
    const std::vector<double>& mu_values,
    const ThermoGrid& grid) {
  std::ofstream out(path);
  if (!out) fail("Failed to write thermodynamic CSV: " + path);
  out << "beta,mu,n,charge_correlation,compressibility,partition_like,log_partition\n";
  out << std::setprecision(15);
  for (std::size_t ibeta = 0; ibeta < beta_values.size(); ++ibeta) {
    for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
      const ThermoPoint& point = grid.at(ibeta).at(imu);
      out << beta_values[ibeta] << ',' << mu_values[imu] << ','
          << point.density << ',' << point.charge_correlation << ','
          << point.compressibility << ',' << point.partition_like << ','
          << point.log_partition << '\n';
    }
  }
}

std::string block_key_text(const BlockKey& key) {
  std::ostringstream out;
  out << "(Nup=" << key.n_up << ",Ndown=" << key.n_down
      << ",mx=" << key.mx << ",my=" << key.my << ')';
  return out.str();
}

std::string status_text(const CompletionStatus& status) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(2)
      << "target_R=" << status.target_samples;
  if (status.target_lanczos_steps > 0) {
    out << " target_m=" << status.target_lanczos_steps;
  }
  out
      << " min_complete_R=" << status.minimum_complete_samples
      << " samples=" << status.durable_sample_records << '/'
      << status.expected_sample_records
      << " blocks=" << status.completed_blocks << '/' << status.expected_blocks
      << " weighted_progress=" << 100.0 * status.weighted_fraction << '%'
      << " state=" << (status.complete ? "complete" : "incomplete");
  if (status.has_next_missing) {
    out << " next=" << block_key_text(status.next_block);
    if (status.next_sample_id >= 0) {
      out << ",sample=" << status.next_sample_id;
      if (status.next_lanczos_steps > 0) {
        out << ",m=" << status.next_lanczos_steps;
      }
    } else {
      out << ",exact-record";
    }
  }
  return out.str();
}

}  // namespace ftlm_checkpoint
