#include "thermo_checkpoint.hpp"

#include <fstream>
#include <algorithm>
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
  int lanczos_steps = 0;
  std::vector<int> sample_values;
  std::vector<int> lanczos_step_values;
  std::vector<double> beta_values;
  double mu_min = -3.0;
  double mu_max = 4.0;
  int mu_count = 281;
  std::string output;
  std::string output_template;
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

std::vector<int> parse_int_list(const std::string& text) {
  std::vector<int> values;
  std::istringstream in(text);
  std::string item;
  while (std::getline(in, item, ',')) {
    if (item.empty()) throw std::runtime_error("Empty integer-list entry.");
    values.push_back(parse_value<int>(item));
  }
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
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
      << "  --samples-list LIST   batch-reduce several R values in one checkpoint read\n"
      << "  --lanczos-steps M     select a saved Lanczos prefix (required for v3)\n"
      << "  --lanczos-steps-list LIST\n"
      << "                        batch-reduce several saved prefixes\n"
      << "  --beta B              one inverse temperature\n"
      << "  --beta-list LIST      comma-separated inverse temperatures\n"
      << "  --mu-min X            minimum chemical potential\n"
      << "  --mu-max X            maximum chemical potential\n"
      << "  --mu-count N          number of chemical-potential points\n"
      << "  --output PATH         thermodynamic CSV output\n"
      << "  --output-template P   batch path containing {m} and {R} placeholders\n"
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
    } else if (arg == "--samples-list") {
      params.sample_values = parse_int_list(next());
    } else if (arg == "--lanczos-steps") {
      params.lanczos_steps = parse_value<int>(next());
    } else if (arg == "--lanczos-steps-list") {
      params.lanczos_step_values = parse_int_list(next());
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
    } else if (arg == "--output-template") {
      params.output_template = next();
    } else if (arg == "--status") {
      params.status = true;
    } else if (arg == "--status-json") {
      params.status_json = next();
    } else {
      throw std::runtime_error("Unknown argument: " + arg);
    }
  }
  if (params.checkpoint.empty()) throw std::runtime_error("--checkpoint is required.");
  const bool batch = !params.sample_values.empty() ||
      !params.lanczos_step_values.empty() || !params.output_template.empty();
  if (!batch && params.samples <= 0) throw std::runtime_error("--samples must be positive.");
  if (batch) {
    if (params.sample_values.empty() && params.samples > 0) {
      params.sample_values = {params.samples};
    }
    if (params.lanczos_step_values.empty() && params.lanczos_steps > 0) {
      params.lanczos_step_values = {params.lanczos_steps};
    }
    if (params.sample_values.empty() || params.lanczos_step_values.empty()) {
      throw std::runtime_error(
          "Batch reduction requires samples and Lanczos-step values.");
    }
    if (params.output_template.empty() ||
        params.output_template.find("{m}") == std::string::npos ||
        params.output_template.find("{R}") == std::string::npos) {
      throw std::runtime_error(
          "--output-template must contain both {m} and {R} in batch mode.");
    }
    if (params.status || !params.status_json.empty() || !params.output.empty()) {
      throw std::runtime_error(
          "Batch reduction cannot be combined with status or single --output mode.");
    }
    for (int samples : params.sample_values) {
      if (samples <= 0) throw std::runtime_error("Batch sample values must be positive.");
    }
    for (int steps : params.lanczos_step_values) {
      if (steps <= 0) throw std::runtime_error("Batch Lanczos steps must be positive.");
    }
  }
  if (params.lanczos_steps < 0) {
    throw std::runtime_error("--lanczos-steps must be positive when provided.");
  }
  if (params.mu_count < 2) throw std::runtime_error("--mu-count must be at least 2.");
  if (!batch && !params.status && params.status_json.empty() && params.output.empty()) {
    throw std::runtime_error("--output is required unless --status or --status-json is used.");
  }
  if ((!params.output.empty() || batch) && params.beta_values.empty()) {
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
    int samples,
    int lanczos_steps) {
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
    for (int sample_id = 0; sample_id < samples; ++sample_id) {
      if (!ftlm_checkpoint::has_sample(data, entry.key, sample_id, lanczos_steps)) {
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
      << "  \"target_m\": ";
  if (status.target_lanczos_steps > 0) out << status.target_lanczos_steps;
  else out << "null";
  out << ",\n"
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
    out << ",\"lanczos_steps\":";
    if (status.next_lanczos_steps > 0) out << status.next_lanczos_steps;
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

std::string zero_padded(int value) {
  std::ostringstream out;
  out << std::setw(3) << std::setfill('0') << value;
  return out.str();
}

void replace_all(std::string& text, const std::string& pattern, const std::string& value) {
  std::size_t position = 0;
  while ((position = text.find(pattern, position)) != std::string::npos) {
    text.replace(position, pattern.size(), value);
    position += value.size();
  }
}

std::string batch_output_path(const std::string& pattern, int samples, int steps) {
  std::string path = pattern;
  replace_all(path, "{m}", zero_padded(steps));
  replace_all(path, "{R}", zero_padded(samples));
  return path;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Params params = parse_args(argc, argv);
    const bool batch = !params.output_template.empty();
    const auto data = ftlm_checkpoint::read_checkpoint(
        params.checkpoint,
        !params.output.empty() || batch,
        batch ? 0 : params.lanczos_steps);
    if (data.format == ftlm_checkpoint::Format::kMissing ||
        data.format == ftlm_checkpoint::Format::kEmpty) {
      throw std::runtime_error("Checkpoint is missing or empty: " + params.checkpoint);
    }
    const auto metadata = ftlm_checkpoint::read_metadata(params.checkpoint + ".meta");
    if (data.format == ftlm_checkpoint::Format::kV2 ||
        data.format == ftlm_checkpoint::Format::kV3) {
      const int expected_format = data.format == ftlm_checkpoint::Format::kV3
          ? ftlm_checkpoint::kDepthCheckpointVersion
          : ftlm_checkpoint::kCheckpointVersion;
      if (metadata.get_int("format") != expected_format ||
          metadata.get_int("rng_version") != ftlm_checkpoint::kRngVersion ||
          metadata.get("rng_algorithm") != ftlm_checkpoint::kRngAlgorithm ||
          metadata.get("seed_derivation") != ftlm_checkpoint::kSeedDerivation) {
        throw std::runtime_error(
            "Unsupported checkpoint RNG/format metadata; use the matching executable.");
      }
    }
    const int lx = metadata.get_int("lx");
    const int ly = metadata.get_int("ly");
    int lanczos_steps = params.lanczos_steps;
    if (data.format == ftlm_checkpoint::Format::kV3 && !batch) {
      if (lanczos_steps <= 0) {
        throw std::runtime_error("--lanczos-steps M is required for a v3 checkpoint.");
      }
    } else if (data.format != ftlm_checkpoint::Format::kV3) {
      const int checkpoint_steps = metadata.get_int("lanczos_steps");
      if (lanczos_steps > 0 && lanczos_steps != checkpoint_steps) {
        throw std::runtime_error(
            "Requested --lanczos-steps does not match this legacy/v2 checkpoint.");
      }
      lanczos_steps = checkpoint_steps;
    }
    const std::size_t threshold = static_cast<std::size_t>(
        metadata.get_u64("exact_block_threshold"));
    const int legacy_samples = data.format == ftlm_checkpoint::Format::kV1
        ? metadata.get_int("samples") : 0;

    if (batch) {
      if (data.format != ftlm_checkpoint::Format::kV3) {
        const int checkpoint_steps = metadata.get_int("lanczos_steps");
        for (int steps : params.lanczos_step_values) {
          if (steps != checkpoint_steps) {
            throw std::runtime_error(
                "A legacy/v2 checkpoint can only batch-reduce its immutable Lanczos depth.");
          }
        }
      }
      const auto mu_values = linspace(params.mu_min, params.mu_max, params.mu_count);
      for (int steps : params.lanczos_step_values) {
        const auto manifest =
            ftlm_checkpoint::build_manifest(lx, ly, threshold, steps);
        for (int samples : params.sample_values) {
          const auto status = ftlm_checkpoint::completion_status(
              data, manifest, samples, steps);
          std::cout << "FTLM_STATUS " << ftlm_checkpoint::status_text(status) << "\n";
          const auto grid = ftlm_checkpoint::reduce_checkpoint(
              data,
              manifest,
              samples,
              lx * ly,
              params.beta_values,
              mu_values,
              legacy_samples,
              steps);
          const std::string output =
              batch_output_path(params.output_template, samples, steps);
          ftlm_checkpoint::write_thermo_csv(
              output, params.beta_values, mu_values, grid);
          std::cout << "Wrote " << output << "\n";
        }
      }
      return 0;
    }

    const auto manifest =
        ftlm_checkpoint::build_manifest(lx, ly, threshold, lanczos_steps);
    const auto status =
        ftlm_checkpoint::completion_status(data, manifest, params.samples, lanczos_steps);
    const auto missing = missing_work(data, manifest, params.samples, lanczos_steps);
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
        legacy_samples,
        lanczos_steps);
    ftlm_checkpoint::write_thermo_csv(
        params.output, params.beta_values, mu_values, grid);
    std::cout << "Wrote " << params.output << "\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << "\n";
    return 1;
  }
}
