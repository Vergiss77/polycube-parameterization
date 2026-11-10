// End-to-end acceptance tests on the three provided models (ctest label "models").
#include <doctest/doctest.h>

#include <cmath>
#include <fstream>

#include <Eigen/Geometry>

#include "helpers.h"
#include "pcp/core/constants.h"
#include "pcp/io/off.h"
#include "pcp/io/tmsh.h"
#include "pcp/pipeline/pipeline.h"
#include "pcp/structure/labeling.h"

using namespace pcp;
namespace fs = std::filesystem;

namespace {

nlohmann::json run(const fs::path& input, const std::string& name, PipelineParams params = {}) {
  PipelineOptions opt;
  opt.input = input;
  opt.output_dir = testutil::temp_dir() / name;
  opt.verbosity = Verbosity::Quiet;
  nlohmann::json stats = run_pipeline(opt, params);
  MESSAGE(name << ": error " << stats.value("polycube_error", -1.0) << " patches "
               << stats.value("patches", -1) << " corners " << stats.value("corners", -1)
               << " arap " << stats.value("arap_mean", -1.0) << " inverted "
               << stats.value("inverted_tets", -1) << " time "
               << stats["timing_seconds"].value("total", -1.0));
  return stats;
}

double rotation_angle_deg(const nlohmann::json& stats) {
  return stats["rotation"]["angle_deg"].get<double>();
}

bool is_identity(const nlohmann::json& m) {
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      if (std::abs(m[r][c].get<double>() - (r == c ? 1.0 : 0.0)) > 1e-12) return false;
  return true;
}

// Acceptance thresholds. TAILSITTER uses relaxed values (a known limitation recorded in
// design.md): error < 5e-3 and at most 0.01% inverted tets.
struct Acceptance {
  double max_error = 1e-3;
  double max_inverted_fraction = 0.0;
};

void check_common(const nlohmann::json& s, Acceptance acc = {}) {
  CHECK(s["polycube_error"].get<double>() < acc.max_error);
  CHECK(s["l1"]["area_rel_diff"].get<double>() < 1e-6);
  CHECK(s["l1"]["min_normal_slack"].get<double>() > 1e-6);
  CHECK(s["topology"]["ok"].get<bool>());
  const long long inverted = s["inverted_tets"].get<long long>();
  const long long tets = s["tets"].get<long long>();
  CHECK(static_cast<double>(inverted) <= acc.max_inverted_fraction * static_cast<double>(tets));
}

}  // namespace

TEST_CASE("model_ispar" * doctest::test_suite("models")) {
  const fs::path input = testutil::models_dir() / "I-SPAR-1B59CD84.off";
  nlohmann::json s = run(input, "model_ispar");
  check_common(s);
  CHECK(s["l1"]["arap_mean"].get<double>() < 1e-3);
  CHECK(s["patches"].get<int>() == 14);
  CHECK(s["corners"].get<int>() == 24);

  // Labels of the deformed boundary agree with the original normals.
  SurfaceMesh surface = read_off(input);
  TetMesh l1 = read_tmsh(testutil::temp_dir() / "model_ispar" / "02_l1.tmsh");
  CHECK(initial_labels(l1.X, l1.B) == initial_labels(surface.V, surface.F));

  // Repeat run is bit-identical.
  run(input, "model_ispar_repeat");
  TetMesh a = read_tmsh(testutil::temp_dir() / "model_ispar" / "04_polycube.tmsh");
  TetMesh b = read_tmsh(testutil::temp_dir() / "model_ispar_repeat" / "04_polycube.tmsh");
  CHECK(a.X == b.X);
}

TEST_CASE("model_clevis" * doctest::test_suite("models")) {
  const fs::path input = testutil::models_dir() / "CLEVIS-D64EEBFA.off";
  nlohmann::json s = run(input, "model_clevis");
  check_common(s);
  const int patches = s["patches"].get<int>();

  // Orientation independence: the same model rotated by 30 degrees about an arbitrary axis.
  SurfaceMesh m = read_off(input);
  Eigen::Matrix3d R =
      Eigen::AngleAxisd(30.0 * kDegToRad, Eigen::Vector3d(1, 2, 3).normalized()).toRotationMatrix();
  m.V = m.V * R.transpose();
  const fs::path rotated = testutil::temp_dir() / "clevis_rotated.off";
  write_off(rotated, m);
  nlohmann::json sr = run(rotated, "model_clevis_rotated");
  CHECK(sr["polycube_error"].get<double>() < 1e-3);
  CHECK(sr["patches"].get<int>() == patches);

  // Rotation disabled: identity is reported (a short schedule is enough to check this).
  PipelineParams p;
  p.l1.global_rotation = false;
  p.l1.max_stages = 1;
  PipelineOptions opt;
  opt.input = input;
  opt.output_dir = testutil::temp_dir() / "model_clevis_norot";
  opt.until = Stage::L1;
  opt.verbosity = Verbosity::Quiet;
  nlohmann::json sn = run_pipeline(opt, p);
  CHECK(is_identity(sn["rotation"]["matrix"]));
}

TEST_CASE("model_tailsitter" * doctest::test_suite("models")) {
  const fs::path input = testutil::models_dir() / "TAILSITTER.off";
  nlohmann::json s = run(input, "model_tailsitter");
  check_common(s, Acceptance{5e-3, 1e-4});
  // The arms run along the 45-degree diagonals: the global rotation is 45 degrees about Z.
  const double angle = rotation_angle_deg(s);
  const double off = std::abs(std::fmod(angle, 90.0) - 45.0);
  CHECK(off < 5.0);
  auto axis = s["rotation"]["axis"];
  CHECK(std::abs(std::abs(axis[2].get<double>()) - 1.0) < 1e-2);
  // Genus preserved by the patch structure.
  CHECK(s["topology"]["euler_complex"].get<long long>() == -8);
}

TEST_CASE("model_tailsitter_beta" * doctest::test_suite("models")) {
  const fs::path input = testutil::models_dir() / "TAILSITTER.off";
  PipelineParams lo, hi;
  lo.l1.beta = 0.2;
  hi.l1.beta = 4.0;
  nlohmann::json a = run(input, "model_tailsitter_beta_lo", lo);
  nlohmann::json b = run(input, "model_tailsitter_beta_hi", hi);
  CHECK(b["corners"].get<int>() <= a["corners"].get<int>());
}
