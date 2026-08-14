#include "thermo_checkpoint.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct Params {
  std::string checkpoint;
  int samples = 0;
  std::vector<double> beta_values;
  double mu_min = -3.0;
  double mu_max = 4.0;
  int mu_count = 281;
  std::string output;
  bool status = false;
  std::string status_json;
};

template <typename T>
T parse_value(const std::string& text) {
  std::istringstream in(text);
  T value{};
  in >> value;
  if (!in || !in.eof()) throw std::runtime_error("Failed to parse value: " + text);
  return value;
}

std::vector<double> parse_double_list(const std::string& text) {
  std::vector<double> values;
  std::istringstream in(text);
  std::string item;
  while (std::getline(in, item, ',')) {
    if (item.empty()) throw std::runtime_error("Empty beta-list entry.");
    values.push_back(parse_value<double>(item));
  }
  return values;
}

std::vector<double> linspace(double first, double last, int count) {
  std::vector<double> values(static_cast<std::size_t>(count));
  const double step = (last - first) / static_cast<double>(count - 1);
  for (int i = 0; i < count; ++i) {
    values[static_cast<std::size_t>(i)] = first + step * static_cast<double>(i);
  }
  return values;
}

void print_help(const char* argv0) {
  std::cout
      << "Usage: " << argv0 << " --checkpoint PATH --samples R [options]\n"
      << "  --beta B              one inverse temperature\n"
      << "  --beta-list LIST      comma-separated inverse temperatures\n"
      << "  --mu-min X            minimum chemical potential\n"
      << "  --mu-max X            maximum chemical potential\n"
      << "  --mu-count N          number of chemical-potential points\n"
      << "  --output PATH         thermodynamic CSV output\n"
      << "  --status              print checkpoint completion status\n"
      << "  --status-json PATH    write machine-readable status JSON\n";
}

Params parse_args(int argc, char** argv) {
  Params params;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto next = [&]() {
      if (++i >= argc) throw std::runtime_error("Missing value for " + arg);
      return std::string(argv[i]);
    };
    if (arg == "--help") {
      print_help(argv[0]);
      std::exit(0);
    } else if (arg == "--checkpoint") {
      params.checkpoint = next();
    } else if (arg == "--samples") {
      params.samples = parse_value<int>(next());
    } else if (arg == "--beta") {
      params.beta_values = {parse_value<double>(next())};
    } else if (arg == "--beta-list") {
      params.beta_values = parse_double_list(next());
    } else if (arg == "--mu-min") {
      params.mu_min = parse_value<double>(next());
    } else if (arg == "--mu-max") {
      params.mu_max = parse_value<double>(next());
    } else if (arg == "--mu-count") {
      params.mu_count = parse_value<int>(next());
    } else if (arg == "--output") {
      params.output = next();
    } else if (arg == "--status") {
      params.status = true;
    } else if (arg == "--status-json") {
      params.status_json = next();
    } else {
      throw std::runtime_error("Unknown argument: " + arg);
    }
  }
  if (params.checkpoint.empty()) throw std::runtime_error("--checkpoint is required.");
  if (params.samples <= 0) throw std::runtime_error("--samples must be positive.");
  if (params.mu_count < 2) throw std::runtime_error("--mu-count must be at least 2.");
  if (!params.status && params.status_json.empty() && params.output.empty()) {
    throw std::runtime_error("--output is required unless --status or --status-json is used.");
  }
  if (!params.output.empty() && params.beta_values.empty()) {
    throw std::runtime_error("--beta or --beta-list is required for reduction.");
  }
  return params;
}

struct MissingWork {
  ftlm_checkpoint::BlockKey key;
  bool exact = false;
  std::vector<int> sample_ids;
};

std::vector<MissingWork> missing_work(
    const ftlm_checkpoint::CheckpointData& data,
    const std::vector<ftlm_checkpoint::ManifestEntry>& manifest,
    int samples) {
  std::vector<MissingWork> missing;
  for (const auto& entry : manifest) {
    if (data.format == ftlm_checkpoint::Format::kV1) {
      if (data.legacy_blocks.find(entry.key) == data.legacy_blocks.end()) {
        missing.push_back({entry.key, true, {}});
      }
      continue;
    }
    if (entry.exact) {
      if (data.exact_blocks.find(entry.key) == data.exact_blocks.end()) {
        missing.push_back({entry.key, true, {}});
      }
      continue;
    }
    MissingWork work;
    work.key = entry.key;
    const auto block = data.samples.find(entry.key);
    for (int sample_id = 0; sample_id < samples; ++sample_id) {
      if (block == data.samples.end() ||
          block->second.find(sample_id) == block->second.end()) {
        work.sample_ids.push_back(sample_id);
      }
    }
    if (!work.sample_ids.empty()) missing.push_back(std::move(work));
  }
  return missing;
}

std::string compact_ids(const std::vector<int>& ids) {
  std::ostringstream out;
  for (std::size_t begin = 0; begin < ids.size();) {
    std::size_t end = begin;
    while (end + 1 < ids.size() && ids[end + 1] == ids[end] + 1) ++end;
    if (begin > 0) out << ',';
    out << ids[begin];
    if (end > begin) out << '-' << ids[end];
    begin = end + 1;
  }
  return out.str();
}

void print_missing_work(const std::vector<MissingWork>& missing) {
  for (const MissingWork& work : missing) {
    std::cout << "FTLM_MISSING " << ftlm_checkpoint::block_key_text(work.key) << ' ';
    if (work.exact) {
      std::cout << "exact_record\n";
    } else {
      std::cout << "sample_ids=" << compact_ids(work.sample_ids) << "\n";
    }
  }
}

void write_status_json(
    const std::string& path,
    const ftlm_checkpoint::CompletionStatus& status,
    const std::vector<MissingWork>& missing) {
  if (path.empty()) return;
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Failed to write status JSON: " + path);
  out << std::setprecision(15)
      << "{\n"
      << "  \"target_R\": " << status.target_samples << ",\n"
      << "  \"minimum_complete_R\": " << status.minimum_complete_samples << ",\n"
      << "  \"durable_sample_records\": " << status.durable_sample_records << ",\n"
      << "  \"expected_sample_records\": " << status.expected_sample_records << ",\n"
      << "  \"completed_blocks\": " << status.completed_blocks << ",\n"
      << "  \"expected_blocks\": " << status.expected_blocks << ",\n"
      << "  \"weighted_progress\": " << status.weighted_fraction << ",\n"
      << "  \"complete\": " << (status.complete ? "true" : "false") << ",\n"
      << "  \"next_missing\": ";
  if (status.has_next_missing) {
    out << "{\"n_up\":" << status.next_block.n_up
        << ",\"n_down\":" << status.next_block.n_down
        << ",\"mx\":" << status.next_block.mx
        << ",\"my\":" << status.next_block.my
        << ",\"sample_id\":";
    if (status.next_sample_id >= 0) out << status.next_sample_id;
    else out << "null";
    out << ",\"exact\":" << (status.next_sample_id < 0 ? "true" : "false") << "},\n";
  } else {
    out << "null,\n";
  }
  out << "  \"missing_blocks\": [";
  for (std::size_t i = 0; i < missing.size(); ++i) {
    const MissingWork& work = missing[i];
    if (i > 0) out << ',';
    out << "\n    {\"n_up\":" << work.key.n_up
        << ",\"n_down\":" << work.key.n_down
        << ",\"mx\":" << work.key.mx
        << ",\"my\":" << work.key.my
        << ",\"exact\":" << (work.exact ? "true" : "false")
        << ",\"missing_sample_ids\":[";
    for (std::size_t j = 0; j < work.sample_ids.size(); ++j) {
      if (j > 0) out << ',';
      out << work.sample_ids[j];
    }
    out << "]}";
  }
  if (!missing.empty()) out << '\n';
  out << "  ]\n";
  out << "}\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Params params = parse_args(argc, argv);
    const auto data = ftlm_checkpoint::read_checkpoint(
        params.checkpoint, !params.output.empty());
    if (data.format == ftlm_checkpoint::Format::kMissing ||
        data.format == ftlm_checkpoint::Format::kEmpty) {
      throw std::runtime_error("Checkpoint is missing or empty: " + params.checkpoint);
    }
    const auto metadata = ftlm_checkpoint::read_metadata(params.checkpoint + ".meta");
    if (data.format == ftlm_checkpoint::Format::kV2) {
      if (metadata.get_int("format") != ftlm_checkpoint::kCheckpointVersion ||
          metadata.get_int("rng_version") != ftlm_checkpoint::kRngVersion ||
          metadata.get("rng_algorithm") != ftlm_checkpoint::kRngAlgorithm ||
          metadata.get("seed_derivation") != ftlm_checkpoint::kSeedDerivation) {
        throw std::runtime_error(
            "Unsupported v2 checkpoint RNG/format metadata; use the matching executable.");
      }
    }
    const int lx = metadata.get_int("lx");
    const int ly = metadata.get_int("ly");
    const int lanczos_steps = metadata.get_int("lanczos_steps");
    const std::size_t threshold = static_cast<std::size_t>(
        metadata.get_u64("exact_block_threshold"));
    const int legacy_samples = data.format == ftlm_checkpoint::Format::kV1
        ? metadata.get_int("samples") : 0;
    const auto manifest =
        ftlm_checkpoint::build_manifest(lx, ly, threshold, lanczos_steps);
    const auto status =
        ftlm_checkpoint::completion_status(data, manifest, params.samples);
    const auto missing = missing_work(data, manifest, params.samples);
    std::cout << "FTLM_STATUS " << ftlm_checkpoint::status_text(status) << "\n";
    if (params.status) print_missing_work(missing);
    write_status_json(params.status_json, status, missing);
    if ((params.status || !params.status_json.empty()) && params.output.empty()) return 0;

    const auto mu_values = linspace(params.mu_min, params.mu_max, params.mu_count);
    const auto grid = ftlm_checkpoint::reduce_checkpoint(
        data,
        manifest,
        params.samples,
        lx * ly,
        params.beta_values,
        mu_values,
        legacy_samples);
    ftlm_checkpoint::write_thermo_csv(
        params.output, params.beta_values, mu_values, grid);
    std::cout << "Wrote " << params.output << "\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << "\n";
    return 1;
  }
}
