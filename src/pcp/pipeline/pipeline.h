#pragma once

#include <filesystem>
#include <functional>
#include <optional>
#include <string>

#include <nlohmann/json.hpp>

#include "pcp/blocks/spline_fit.h"
#include "pcp/polycube/l1_solver.h"
#include "pcp/structure/final_map.h"
#include "pcp/structure/height_simplify.h"
#include "pcp/tet/tetrahedralize.h"

namespace pcp {

enum class Stage { Validate = 0, Tet = 1, L1 = 2, Structure = 3, Simplify = 4, Blocks = 5, Fit = 6 };
Stage parse_stage(const std::string& s);  // validate | tet | l1 | structure | simplify | blocks | fit
const char* stage_name(Stage s);

// All tunable parameters of the pipeline.
struct PipelineParams {
  TetParams tet;
  L1Params l1;
  FinalMapParams final_map;
  SimplifyParams simplify;
  FitParams fit;
  std::string feature_edges;  // optional path of a feature edge file

  nlohmann::json to_json() const;
  // Overwrites the fields present in j ({"tet": {...}, "l1": {...}, "final_map": {...},
  // "simplify": {...}, "fit": {...}, "feature_edges": "..."}). Throws pcp::Error on unknown keys.
  void update_from_json(const nlohmann::json& j);
};

enum class Verbosity { Quiet, Normal, Verbose };

struct PipelineOptions {
  std::filesystem::path input;
  std::filesystem::path output_dir;
  Stage until = Stage::Fit;
  std::optional<std::filesystem::path> resume_from;  // saved tet mesh (.tmsh)
  Verbosity verbosity = Verbosity::Normal;
};

// Runs validation, tetrahedralization, l1 deformation, structure extraction (labels, patch
// graph, final polycube map), height simplification (Chen 2019, unless disabled), block
// decomposition and B-spline fitting, writing all results and stats.json into output_dir.
// Throws pcp::Error on failure; the statistics written so far are still saved.
nlohmann::json run_pipeline(const PipelineOptions& opt, const PipelineParams& params,
                            const std::function<void(const std::string&)>& console = {});

}  // namespace pcp
