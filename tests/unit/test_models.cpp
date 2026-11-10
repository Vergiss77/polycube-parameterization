// End-to-end acceptance tests on the three provided models (ctest label "models").
#include <doctest/doctest.h>

#include <cmath>

#include <Eigen/Geometry>

#include "helpers.h"
#include "pcp/core/constants.h"
#include "pcp/io/gismo_xml.h"
#include "pcp/io/off.h"
#include "pcp/io/text.h"
#include "pcp/io/tmsh.h"
#include "pcp/pipeline/pipeline.h"
#include "pcp/structure/labeling.h"

using namespace pcp;
namespace fs = std::filesystem;

namespace {

nlohmann::json run(const fs::path& input, const std::string& name, PipelineParams params = {},
                   Stage until = Stage::Fit) {
  PipelineOptions opt;
  opt.input = input;
  opt.output_dir = testutil::temp_dir() / name;
  opt.until = until;
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

// Height simplification checks (spec polycube-simplification).
void check_simplify(const nlohmann::json& s, bool strictly_fewer_cuts, Acceptance acc = {}) {
  REQUIRE(s.contains("simplify"));
  const nlohmann::json& sp = s["simplify"];
  REQUIRE(sp["enabled"].get<bool>());
  for (auto it = sp.begin(); it != sp.end(); ++it)
    if (!it.value().is_object()) CHECK_MESSAGE(!it.value().is_null(), "simplify." << it.key() << " is NaN");
  const int before = sp["cuts_before"]["total"].get<int>(), after = sp["cuts_after"]["total"].get<int>();
  if (strictly_fewer_cuts) CHECK(after < before);
  else CHECK(after <= before);
  const double smax = sp["max_stretch"].get<double>();
  CHECK(sp["min_ratio"].get<double>() >= 1.0 / smax - 1e-9);
  CHECK(sp["max_ratio"].get<double>() <= smax + 1e-9);
  CHECK(sp["max_planarity_error_rel"].get<double>() <= 1e-9);
  CHECK(sp["topology"]["ok"].get<bool>());
  const long long inverted = sp["inverted_tets"].get<long long>();
  CHECK(static_cast<double>(inverted) <= acc.max_inverted_fraction * s["tets"].get<double>());
  CHECK(sp["merges"].get<int>() <= sp["merges_proposed"].get<int>());
  CHECK(sp["inverted_kept"].get<long long>() <= sp["inverted_reference"].get<long long>());
  MESSAGE("simplify: cuts " << before << " -> " << after << ", merges " << sp["merges"] << " of "
                            << sp["merges_proposed"] << ", edges " << sp["chain_edges"] << " + "
                            << sp["separation_edges"] << " + " << sp["width_edges"]
                            << ", ratio [" << sp["min_ratio"] << ", " << sp["max_ratio"]
                            << "], arap " << sp["arap_mean"] << " / " << sp["arap_max"]
                            << ", inverted " << inverted << " (before untangle "
                            << sp["inverted_before_untangle"] << ")");
}

// Block decomposition and B-spline fitting checks (specs polycube-block-decomposition and
// bspline-volume-fitting). max_error < 0 skips the max error check; same for surface_rms.
struct FitAcceptance {
  double max_error = -1.0;
  double surface_rms = -1.0;
  bool positive_jacobian = false;
};

void check_blocks_and_fit(const nlohmann::json& s, const fs::path& dir, FitAcceptance acc) {
  REQUIRE(s.contains("blocks"));
  REQUIRE(s.contains("fit"));
  const nlohmann::json& bl = s["blocks"];
  const nlohmann::json& fit = s["fit"];
  for (const nlohmann::json* j : {&bl, &fit})
    for (auto it = j->begin(); it != j->end(); ++it)
      CHECK_MESSAGE(!it.value().is_null(), it.key() << " is NaN");
  CHECK(bl["volume_rel_diff"].get<double>() <= 1e-9);
  // Cut coordinates per axis = simplified cut count per axis.
  if (s["simplify"]["enabled"].get<bool>()) {
    const auto& ca = s["simplify"]["cuts_after"];
    CHECK(bl["cuts"][0].get<int>() == ca["x"].get<int>());
    CHECK(bl["cuts"][1].get<int>() == ca["y"].get<int>());
    CHECK(bl["cuts"][2].get<int>() == ca["z"].get<int>());
  }
  // Every face classified, interior faces symmetric.
  const auto bj = nlohmann::json::parse(text::read_file(dir / "06_blocks.json"));
  const auto& blocks = bj["blocks"];
  CHECK(blocks.size() == bl["blocks"].get<size_t>());
  static const char* side[6] = {"x_lo", "x_hi", "y_lo", "y_hi", "z_lo", "z_hi"};
  for (size_t i = 0; i < blocks.size(); ++i)
    for (int f = 0; f < 6; ++f) {
      const auto& face = blocks[i]["faces"][side[f]];
      CHECK((face.contains("neighbor") != face.contains("patch")));
      if (face.contains("neighbor"))
        CHECK(blocks[face["neighbor"].get<size_t>()]["faces"][side[f ^ 1]]["neighbor"].get<size_t>() == i);
    }
  CHECK(fit["c0_max_gap"].get<double>() <= 1e-12);
  // G+Smo XML read back: patches, interfaces and boundary sides match the block layout.
  const GismoMultiPatch mp = read_gismo_xml(dir / "07_bspline.xml");
  CHECK(mp.volumes.size() == blocks.size());
  CHECK(int(mp.interfaces.size()) == bl["interior_faces"].get<int>());
  CHECK(int(mp.boundary.size()) == bl["boundary_faces"].get<int>());
  if (acc.max_error >= 0) CHECK(fit["max_error"].get<double>() <= acc.max_error);
  if (acc.surface_rms >= 0) CHECK(fit["surface_rms_error"].get<double>() <= acc.surface_rms);
  if (acc.positive_jacobian) CHECK(fit["nonpositive_jacobian"].get<long long>() == 0);
  MESSAGE("blocks " << bl["blocks"] << " (grid cells " << bl["grid_cells"] << ", thinnest "
                    << bl["thinnest_block_rel"] << "); fit: control points " << fit["control_points"]
                    << ", error rms " << fit["rms_error"] << " max " << fit["max_error"] << ", surface rms "
                    << fit["surface_rms_error"] << ", distance max " << fit["surface_distance_max"]
                    << ", scaled Jacobian min " << fit["min_scaled_jacobian"] << " mean "
                    << fit["mean_scaled_jacobian"] << ", non-positive " << fit["nonpositive_jacobian"] << " of "
                    << fit["eval_points"] << " (" << fit["blocks_with_nonpositive"] << " blocks), unlocated "
                    << fit["unlocated_samples"] << ", " << fit["seconds"] << " s");
}

}  // namespace

TEST_CASE("model_ispar" * doctest::test_suite("models")) {
  const fs::path input = testutil::models_dir() / "I-SPAR-1B59CD84.off";
  nlohmann::json s = run(input, "model_ispar");
  check_common(s);
  CHECK(s["l1"]["arap_mean"].get<double>() < 1e-3);
  CHECK(s["patches"].get<int>() == 14);
  CHECK(s["corners"].get<int>() == 24);
  check_simplify(s, false);
  check_blocks_and_fit(s, testutil::temp_dir() / "model_ispar", FitAcceptance{1e-3, -1.0, true});

  // Labels of the deformed boundary agree with the original normals.
  SurfaceMesh surface = read_off(input);
  TetMesh l1 = read_tmsh(testutil::temp_dir() / "model_ispar" / "02_l1.tmsh");
  CHECK(initial_labels(l1.X, l1.B) == initial_labels(surface.V, surface.F));

  // Repeat run is bit-identical.
  run(input, "model_ispar_repeat");
  TetMesh a = read_tmsh(testutil::temp_dir() / "model_ispar" / "04_polycube.tmsh");
  TetMesh b = read_tmsh(testutil::temp_dir() / "model_ispar_repeat" / "04_polycube.tmsh");
  CHECK(a.X == b.X);
  // Simplification is deterministic too.
  auto sa = nlohmann::json::parse(text::read_file(testutil::temp_dir() / "model_ispar" / "05_simplify.json"));
  auto sb = nlohmann::json::parse(text::read_file(testutil::temp_dir() / "model_ispar_repeat" / "05_simplify.json"));
  CHECK(sa["merges"] == sb["merges"]);
  CHECK(sa["patches"] == sb["patches"]);
  // So is the fitted multipatch.
  CHECK(text::read_file(testutil::temp_dir() / "model_ispar" / "07_bspline.xml") ==
        text::read_file(testutil::temp_dir() / "model_ispar_repeat" / "07_bspline.xml"));
}

TEST_CASE("model_clevis" * doctest::test_suite("models")) {
  const fs::path input = testutil::models_dir() / "CLEVIS-D64EEBFA.off";
  nlohmann::json s = run(input, "model_clevis");
  check_common(s);
  check_simplify(s, true);
  check_blocks_and_fit(s, testutil::temp_dir() / "model_clevis", FitAcceptance{-1.0, 5e-3, false});
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
  check_simplify(s, true, Acceptance{5e-3, 1e-4});
  check_blocks_and_fit(s, testutil::temp_dir() / "model_tailsitter", FitAcceptance{-1.0, 5e-3, false});
  // A larger stretch bound merges at least as much (re-run on the saved graph, no re-mapping).
  auto sj = nlohmann::json::parse(text::read_file(testutil::temp_dir() / "model_tailsitter" / "05_simplify.json"));
  const HeightGraph hg = HeightGraph::from_json(sj["graph_normalized"]);
  const SimplifyResult lo = simplify_heights(hg, 1.5), hi = simplify_heights(hg, 3.0);
  MESSAGE("tailsitter cuts with s = 1.5: " << lo.cuts_total_after() << ", s = 3: " << hi.cuts_total_after());
  CHECK(hi.cuts_total_after() <= lo.cuts_total_after());
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
  // This test compares the l1 structure only: stop after the structure stage.
  nlohmann::json a = run(input, "model_tailsitter_beta_lo", lo, Stage::Structure);
  nlohmann::json b = run(input, "model_tailsitter_beta_hi", hi, Stage::Structure);
  CHECK(b["corners"].get<int>() <= a["corners"].get<int>());
}
