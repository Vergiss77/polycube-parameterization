#include "pcp/pipeline/pipeline.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <algorithm>
#include <string>

#include <Eigen/Geometry>

#include "pcp/blocks/block_layout.h"
#include "pcp/core/constants.h"
#include "pcp/io/gismo_xml.h"
#include "pcp/geom/normalize.h"
#include "pcp/io/off.h"
#include "pcp/io/text.h"
#include "pcp/io/tmsh.h"
#include "pcp/io/validate.h"
#include "pcp/io/vtu.h"
#include "pcp/polycube/feature_edges.h"
#include "pcp/polycube/rotation.h"
#include "pcp/structure/labeling.h"
#include "pcp/structure/patch_graph.h"

namespace pcp {

namespace fs = std::filesystem;

Stage parse_stage(const std::string& s) {
  if (s == "validate") return Stage::Validate;
  if (s == "tet") return Stage::Tet;
  if (s == "l1") return Stage::L1;
  if (s == "structure") return Stage::Structure;
  if (s == "simplify") return Stage::Simplify;
  if (s == "blocks") return Stage::Blocks;
  if (s == "fit") return Stage::Fit;
  throw Error("unknown stage '" + s + "' (expected validate, tet, l1, structure, simplify, blocks or fit)");
}

const char* stage_name(Stage s) {
  switch (s) {
    case Stage::Validate: return "validate";
    case Stage::Tet: return "tet";
    case Stage::L1: return "l1";
    case Stage::Structure: return "structure";
    case Stage::Simplify: return "simplify";
    case Stage::Blocks: return "blocks";
    default: return "fit";
  }
}

nlohmann::json PipelineParams::to_json() const {
  return {{"tet", tet.to_json()},
          {"l1", l1.to_json()},
          {"final_map", final_map.to_json()},
          {"simplify", simplify.to_json()},
          {"fit", fit.to_json()},
          {"feature_edges", feature_edges}};
}

void PipelineParams::update_from_json(const nlohmann::json& j) {
  for (auto it = j.begin(); it != j.end(); ++it) {
    const std::string& k = it.key();
    if (k == "tet")
      tet.update_from_json(it.value());
    else if (k == "l1")
      l1.update_from_json(it.value());
    else if (k == "final_map")
      final_map.update_from_json(it.value());
    else if (k == "simplify")
      simplify.update_from_json(it.value());
    else if (k == "fit")
      fit.update_from_json(it.value());
    else if (k == "feature_edges")
      feature_edges = it.value().get<std::string>();
    else
      throw Error("unknown parameter section '" + k + "'");
  }
}

namespace {

class Logger {
 public:
  Logger(const fs::path& file, Verbosity v, std::function<void(const std::string&)> console)
      : v_(v), console_(std::move(console)) {
    f_ = std::fopen(file.string().c_str(), "wb");
    if (!f_) throw Error("cannot write log file: " + file.string());
  }
  ~Logger() {
    if (f_) std::fclose(f_);
  }
  void info(const std::string& s) { emit(s, v_ != Verbosity::Quiet); }
  void detail(const std::string& s) { emit(s, v_ == Verbosity::Verbose); }

 private:
  void emit(const std::string& s, bool to_console) {
    std::fputs(s.c_str(), f_);
    std::fputc('\n', f_);
    std::fflush(f_);
    if (to_console && console_) console_(s);
  }
  FILE* f_ = nullptr;
  Verbosity v_;
  std::function<void(const std::string&)> console_;
};

void write_json(const fs::path& p, const nlohmann::json& j) { text::write_file(p, j.dump(2) + "\n"); }

double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

nlohmann::json matrix_json(const Eigen::Matrix3d& R) {
  return {{R(0, 0), R(0, 1), R(0, 2)}, {R(1, 0), R(1, 1), R(1, 2)}, {R(2, 0), R(2, 1), R(2, 2)}};
}

}  // namespace

nlohmann::json run_pipeline(const PipelineOptions& opt, const PipelineParams& params,
                            const std::function<void(const std::string&)>& console) {
  if (!fs::exists(opt.input)) throw Error("input file does not exist: " + opt.input.string());
  fs::create_directories(opt.output_dir);
  const fs::path out = opt.output_dir;
  Logger log(out / "log.txt", opt.verbosity, console);
  write_json(out / "params.json", params.to_json());

  nlohmann::json stats;
  nlohmann::json timing;
  stats["input"] = opt.input.string();
  auto save_stats = [&] {
    stats["timing_seconds"] = timing;
    write_json(out / "stats.json", stats);
  };
  const auto t_all = std::chrono::steady_clock::now();

  try {
    // -- validate
    auto t0 = std::chrono::steady_clock::now();
    SurfaceMesh surface = read_off(opt.input);
    SurfaceReport rep = validate_surface(surface);
    write_json(out / "00_validate.json", rep.to_json());
    stats["input_vertices"] = rep.n_vertices;
    stats["input_faces"] = rep.n_faces;
    stats["genus"] = rep.genus;
    timing["validate"] = seconds_since(t0);
    log.info("validate: " + std::to_string(rep.n_vertices) + " vertices, " +
             std::to_string(rep.n_faces) + " faces, genus " + std::to_string(rep.genus));
    if (!rep.ok) {
      std::string msg = "input surface is invalid:";
      for (const auto& e : rep.errors) msg += " " + e + ";";
      throw Error(msg);
    }
    save_stats();
    if (opt.until == Stage::Validate) return stats;

    // -- tetrahedralize (or resume)
    t0 = std::chrono::steady_clock::now();
    TetMesh tets;
    if (opt.resume_from) {
      tets = read_tmsh(*opt.resume_from);
      if (tets.B.rows() != surface.F.rows())
        throw Error("resume mesh boundary does not match the input surface");
      log.info("tet: resumed from " + opt.resume_from->string());
    } else {
      TetReport tr;
      tets = tetrahedralize(surface, params.tet, &tr);
      stats["tetrahedralization"] = tr.to_json();
      log.detail("tet: TetGen switches " + tr.switches + ", dihedral range [" +
                 std::to_string(tr.min_dihedral_deg) + ", " + std::to_string(tr.max_dihedral_deg) +
                 "] deg, " + std::to_string(tr.n_small_dihedral) + " tets below 5 deg");
      write_tmsh(out / "01_tet.tmsh", tets);
      write_vtu_tets(out / "01_tet.vtu", tets.X, tets.T);
    }
    stats["tets"] = tets.T.rows();
    stats["nodes"] = tets.X.rows();
    timing["tet"] = seconds_since(t0);
    log.info("tet: " + std::to_string(tets.X.rows()) + " nodes, " + std::to_string(tets.T.rows()) +
             " tets");
    save_stats();
    if (opt.until == Stage::Tet) return stats;

    // -- l1 polycube deformation
    t0 = std::chrono::steady_clock::now();
    Normalization norm = normalization_for(tets.X, tets.B);
    const MatX3d rest = norm.apply(tets.X);
    DeformMesh mesh = DeformMesh::build(rest, tets.T, tets.B);
    std::vector<FeatureEdgeEnergy::Edge> edges;
    if (!params.feature_edges.empty()) edges = read_feature_edges(params.feature_edges, mesh.topo);
    L1Result l1 = solve_l1_polycube(rest, mesh, params.l1, edges,
                                    [&](const std::string& s) { log.info("l1: " + s); });
    TetMesh l1_mesh = tets;
    l1_mesh.X = norm.inverse(l1.X);
    write_tmsh(out / "02_l1.tmsh", l1_mesh);
    write_vtu_tets(out / "02_l1.vtu", l1_mesh.X, tets.T, {CellField::scalar("arap", l1.arap_per_tet)});
    Eigen::Vector3d rot_axis = Eigen::AngleAxisd(l1.rotation).axis();
    stats["rotation"] = {{"matrix", matrix_json(l1.rotation)},
                         {"angle_deg", rotation_angle(l1.rotation) * kRadToDeg},
                         {"axis", {rot_axis.x(), rot_axis.y(), rot_axis.z()}}};
    nlohmann::json stages = nlohmann::json::array();
    for (const StageLog& s : l1.stage_logs)
      stages.push_back({{"stage", s.stage}, {"alpha", s.alpha}, {"eps", s.eps}, {"error", s.error},
                        {"arap", s.arap}, {"iterations", s.iterations}, {"seconds", s.seconds},
                        {"cg_iterations", s.cg_iterations}, {"cg_max", s.cg_max},
                        {"end", s.end_reason}});
    nlohmann::json stalls = nlohmann::json::array();
    for (const StallEvent& e : l1.stalls)
      stalls.push_back({{"stage", e.stage}, {"iteration", e.iteration},
                        {"slim_triangles", e.slim_triangles}, {"gradient_norm", e.gradient_norm}});
    stats["l1"] = {{"polycube_error", l1.error},
                   {"reached_threshold", l1.reached_threshold},
                   {"arap_mean", l1.arap_mean},
                   {"arap_max", l1.arap_max},
                   {"area_rel_diff", l1.area_rel_diff},
                   {"min_normal_slack", l1.min_normal_slack},
                   {"inverted_tets", l1.inverted_tets},
                   {"stages", l1.stages},
                   {"iterations", l1.iterations},
                   {"stage_log", stages},
                   {"stalls", stalls},
                   {"solver_timing", l1.timing.to_json()}};
    stats["polycube_error"] = l1.error;
    stats["stages"] = l1.stages;
    stats["iterations"] = l1.iterations;
    timing["l1"] = seconds_since(t0);
    log.detail("l1: solver timing " + l1.timing.to_json().dump() + ", ARAP mean " +
               std::to_string(l1.arap_mean) + ", area change " + std::to_string(l1.area_rel_diff));
    log.info("l1: polycube error " + std::to_string(l1.error) + " after " +
             std::to_string(l1.iterations) + " iterations");
    save_stats();
    if (opt.until == Stage::L1) return stats;

    // -- structure: labels, cleanup, patch graph, final map
    t0 = std::chrono::steady_clock::now();
    Eigen::VectorXi labels = initial_labels(l1.X, tets.B);
    CleanupReport cleanup = cleanup_labels(l1.X, tets.B, mesh.topo, labels);
    PatchGraph graph = build_patch_graph(l1.X, tets.B, mesh.topo, labels);
    TopologyCheck topo_check = check_polycube_topology(tets.B, graph);
    stats["cleanup"] = cleanup.to_json();
    stats["topology"] = topo_check.to_json();
    log.detail("structure: cleanup " + cleanup.to_json().dump() + ", Euler characteristic " +
               std::to_string(topo_check.euler_complex) + " (surface " +
               std::to_string(topo_check.euler_surface) + "), non-disk patches " +
               std::to_string(topo_check.non_disk_patches));
    write_vtu_surface(out / "03_labels_input.vtu", surface.V, surface.F,
                      {CellField::integer("label", labels), CellField::integer("patch", graph.patches.face_patch)});
    log.info("structure: " + std::to_string(graph.patches.label.size()) + " patches, " +
             std::to_string(graph.corners.size()) + " corners, " + std::to_string(graph.chains.size()) +
             " chains");
    if (!topo_check.ok) {
      save_stats();
      std::string msg = "polycube topology check failed:";
      for (const auto& v : topo_check.violations) msg += " " + v + ";";
      throw Error(msg);
    }

    FinalMapResult fm = compute_final_map(rest, l1.X, tets.T, tets.B, graph, graph.height,
                                          params.final_map);
    // Patch graph with the exact heights of the final polycube (in output units).
    const MatX3d polycube = norm.inverse(fm.X);
    PatchGraph final_graph = build_patch_graph(polycube, tets.B, mesh.topo, labels);
    write_json(out / "03_patch_graph.json", final_graph.to_json(polycube));
    write_vtu_surface(out / "03_labels_polycube.vtu", polycube, tets.B,
                      {CellField::integer("label", labels), CellField::integer("patch", graph.patches.face_patch)});
    TetMesh pc_mesh = tets;
    pc_mesh.X = polycube;
    write_tmsh(out / "04_polycube.tmsh", pc_mesh);
    Eigen::VectorXi inverted(tets.T.rows());
    for (Index t = 0; t < tets.T.rows(); ++t) inverted(t) = tet_volume(polycube, tets.T, t) > 0.0 ? 0 : 1;
    write_vtu_tets(out / "04_polycube.vtu", polycube, tets.T,
                   {CellField::scalar("arap", fm.arap_per_tet), CellField::integer("inverted", inverted)});

    const double diag = (polycube.colwise().maxCoeff() - polycube.colwise().minCoeff()).norm();
    stats["final_map"] = {{"arap_mean", fm.arap_mean},
                          {"arap_max", fm.arap_max},
                          {"inverted_tets", fm.inverted_tets},
                          {"inverted_before_untangle", fm.inverted_before_untangle},
                          {"iterations", fm.iterations},
                          {"max_planarity_error_rel", fm.max_planarity_error / norm.scale / diag},
                          {"polycube_error", polycube_error(fm.X, tets.B)}};
    stats["patches"] = graph.patches.label.size();
    stats["corners"] = graph.corners.size();
    stats["chains"] = graph.chains.size();
    stats["arap_mean"] = fm.arap_mean;
    stats["arap_max"] = fm.arap_max;
    stats["inverted_tets"] = fm.inverted_tets;
    timing["structure"] = seconds_since(t0);
    log.detail("final map: " + std::to_string(fm.iterations) + " ARAP iterations in " +
               std::to_string(fm.seconds) + " s");
    log.info("final map: arap mean " + std::to_string(fm.arap_mean) + ", max " +
             std::to_string(fm.arap_max) + ", inverted tets " + std::to_string(fm.inverted_tets));
    save_stats();
    if (opt.until == Stage::Structure) {
      timing["total"] = seconds_since(t_all);
      save_stats();
      return stats;
    }

    // -- simplify: Chen 2019 height simplification and re-mapping with the new heights
    // Axis-aligned (normalized) polycube and patch heights used by the later stages.
    MatX3d poly_X = fm.X;
    std::vector<double> poly_heights = fm.heights;
    if (!params.simplify.enabled) {
      stats["simplify"] = {{"enabled", false}};
      log.info("simplify: disabled");
    } else {
      t0 = std::chrono::steady_clock::now();
      // Patch graph of the exact (normalized) polycube: chain ends follow its heights.
      PatchGraph g_exact = build_patch_graph(fm.X, tets.B, mesh.topo, labels);
      SimplifiedPolycube sp = simplify_and_remap(rest, fm.X, tets.T, tets.B, g_exact, fm.heights,
                                                 params.simplify, params.final_map);
      const SimplifyResult& sr = sp.simplify;
      const MatX3d simplified = norm.inverse(sp.map.X);
      auto out_height = [&](int axis, double h) { return h / norm.scale + norm.center(axis); };

      nlohmann::json patches_j = nlohmann::json::array();
      for (size_t p = 0; p < g_exact.patches.label.size(); ++p) {
        const int axis = g_exact.patches.label[p] / 2;
        patches_j.push_back({{"id", p},
                             {"axis", axis},
                             {"label", g_exact.patches.label[p]},
                             {"group", sr.group_of_patch[p]},
                             {"height_before", out_height(axis, fm.heights[p])},
                             {"height_after", out_height(axis, sr.patch_height[p])}});
      }
      nlohmann::json merges_j = nlohmann::json::array();
      for (const MergeRecord& m : sr.merges) {
        const int axis = sp.graph.groups[static_cast<size_t>(m.a)].axis;
        merges_j.push_back({{"a", m.a}, {"b", m.b}, {"axis", axis},
                            {"height", out_height(axis, m.height)}, {"cost", m.cost}});
      }
      nlohmann::json sj = sr.to_json();
      sj["merges"] = merges_j;
      sj["patches"] = patches_j;
      sj["graph_normalized"] = sp.graph.to_json();
      sj["max_stretch"] = params.simplify.max_stretch;
      write_json(out / "05_simplify.json", sj);

      TetMesh s_mesh = tets;
      s_mesh.X = simplified;
      write_tmsh(out / "05_polycube_simplified.tmsh", s_mesh);
      Eigen::VectorXi inv2(tets.T.rows());
      for (Index t = 0; t < tets.T.rows(); ++t) inv2(t) = tet_volume(simplified, tets.T, t) > 0.0 ? 0 : 1;
      write_vtu_tets(out / "05_polycube_simplified.vtu", simplified, tets.T,
                     {CellField::scalar("arap", sp.map.arap_per_tet), CellField::integer("inverted", inv2)});

      const double sdiag = (simplified.colwise().maxCoeff() - simplified.colwise().minCoeff()).norm();
      nlohmann::json st = {
          {"enabled", true},
          {"max_stretch", params.simplify.max_stretch},
          {"cuts_before", sj["cuts_before"]},
          {"cuts_after", sj["cuts_after"]},
          {"merges", sr.merges.size()},
          {"merges_proposed", sr.merges_proposed},
          {"chain_edges", sr.chain_edges},
          {"separation_edges", sr.separation_edges},
          {"width_edges", sr.width_edges},
          {"inverted_reference", sp.inverted_reference},
          {"inverted_kept", sp.inverted_kept},
          {"probe_maps", sp.probe_maps},
          {"chains_without_ends", sp.graph.chains_without_ends},
          {"min_ratio", sr.min_ratio},
          {"max_ratio", sr.max_ratio},
          {"arap_mean", sp.map.arap_mean},
          {"arap_max", sp.map.arap_max},
          {"inverted_tets", sp.map.inverted_tets},
          {"inverted_before_untangle", sp.map.inverted_before_untangle},
          {"max_planarity_error_rel", sp.map.max_planarity_error / norm.scale / sdiag},
          {"topology", sp.topology.to_json()},
          {"seconds_simplify", sp.seconds_simplify},
          {"seconds_map", sp.map.seconds}};
      stats["simplify"] = st;
      timing["simplify"] = seconds_since(t0);
      poly_X = sp.map.X;
      poly_heights = sr.patch_height;
      log.info("simplify: cuts " + std::to_string(sr.cuts_total_before()) + " -> " +
               std::to_string(sr.cuts_total_after()) + " (" + std::to_string(sr.merges.size()) +
               " of " + std::to_string(sr.merges_proposed) + " merges kept, " +
               std::to_string(sr.chain_edges) + " chain + " + std::to_string(sr.separation_edges) +
               " separation + " + std::to_string(sr.width_edges) + " width edges, ratio [" +
               std::to_string(sr.min_ratio) + ", " + std::to_string(sr.max_ratio) + "])");
      log.detail("simplify: rollback inverted tets " + std::to_string(sp.inverted_kept) + " (reference " +
                 std::to_string(sp.inverted_reference) + ", " + std::to_string(sp.probe_maps) +
                 " probe maps)");
      log.info("simplify: re-mapped polycube arap mean " + std::to_string(sp.map.arap_mean) + ", max " +
               std::to_string(sp.map.arap_max) + ", inverted tets " +
               std::to_string(sp.map.inverted_tets));
      if (!sp.topology.ok) {
        std::string msg = "simplified polycube topology check failed:";
        for (const auto& v : sp.topology.violations) msg += " " + v + ";";
        throw Error(msg);
      }
    }
    save_stats();
    if (opt.until == Stage::Simplify) {
      timing["total"] = seconds_since(t_all);
      save_stats();
      return stats;
    }

    // -- blocks: cut the polycube along the cut planes into cuboid blocks
    t0 = std::chrono::steady_clock::now();
    const auto cuts = cut_coordinates(poly_X, graph.patches.label, poly_heights);
    const BlockLayout layout =
        build_block_layout(poly_X, tets.B, graph.patches.face_patch, graph.patches.label, cuts);
    const Eigen::VectorXi tet_block = assign_tets_to_blocks(layout, poly_X, tets.T);
    write_json(out / "06_blocks.json", layout.to_json(norm.scale, norm.center));
    write_vtu_tets(out / "06_blocks_polycube.vtu", norm.inverse(poly_X), tets.T,
                   {CellField::integer("block", tet_block)});
    write_vtu_tets(out / "06_blocks_input.vtu", tets.X, tets.T, {CellField::integer("block", tet_block)});
    double tet_vol = 0.0;
    for (Index t = 0; t < tets.T.rows(); ++t) tet_vol += tet_volume(poly_X, tets.T, t);
    const double pdiag = (poly_X.colwise().maxCoeff() - poly_X.colwise().minCoeff()).norm();
    double thinnest = std::numeric_limits<double>::infinity();
    for (const Block& b : layout.blocks) thinnest = std::min(thinnest, (b.hi - b.lo).minCoeff());
    stats["blocks"] = {{"blocks", layout.blocks.size()},
                       {"cuts", {cuts[0].size(), cuts[1].size(), cuts[2].size()}},
                       {"grid_cells", layout.grid_cells()},
                       {"interior_faces", layout.interior_faces()},
                       {"boundary_faces", layout.boundary_faces()},
                       {"volume_rel_diff", std::abs(layout.blocks_volume() - tet_vol) / std::abs(tet_vol)},
                       {"thinnest_block_rel", thinnest / pdiag},
                       {"seconds", seconds_since(t0)}};
    timing["blocks"] = seconds_since(t0);
    log.info("blocks: " + std::to_string(layout.blocks.size()) + " blocks in a " +
             std::to_string(layout.intervals(0)) + " x " + std::to_string(layout.intervals(1)) + " x " +
             std::to_string(layout.intervals(2)) + " grid, " + std::to_string(layout.interior_faces()) +
             " interior and " + std::to_string(layout.boundary_faces()) + " boundary faces");
    save_stats();
    if (opt.until == Stage::Blocks) {
      timing["total"] = seconds_since(t_all);
      save_stats();
      return stats;
    }

    // -- fit: one trivariate B-spline volume per block, shared control points (C0)
    t0 = std::chrono::steady_clock::now();
    const FitParams& fp = params.fit;
    const SplineSpace space = build_spline_space(layout, fp.degree, fp.element_size * pdiag);
    const FitSamples samples = sample_blocks(layout, space, fp.samples_per_element, poly_X, tets.T, tets.X);
    const double t_sample = seconds_since(t0);
    const FitSolution sol = fit_volumes(layout, space, samples, fp.smoothing);
    const double t_solve = seconds_since(t0) - t_sample;
    const double in_diag = (tets.X.colwise().maxCoeff() - tets.X.colwise().minCoeff()).norm();
    const FitQuality quality = evaluate_fit(layout, space, sol, samples, tets.X, tets.B, in_diag);
    write_gismo_xml(out / "07_bspline.xml", gismo_topology(layout, sol.volumes));
    const SampleMesh smesh = sample_mesh(layout, space, sol);
    write_vtu_hexes(out / "07_bspline_samples.vtu", smesh.X, smesh.H, {CellField::integer("block", smesh.block)},
                    {PointField::scalar("scaled_jacobian", smesh.scaled_jacobian)});
    long long ncoefs = 0;
    for (const TensorBSpline3& v : sol.volumes) ncoefs += static_cast<long long>(v.coefs.size());
    nlohmann::json fj = quality.to_json();
    fj["degree"] = fp.degree;
    fj["control_points"] = space.n_unknowns;
    fj["control_points_per_patch_total"] = ncoefs;
    fj["samples"] = samples.block.size();
    fj["unlocated_samples"] = samples.unlocated;
    fj["seconds_sampling"] = t_sample;
    fj["seconds_solve"] = t_solve;
    fj["seconds"] = seconds_since(t0);
    stats["fit"] = fj;
    timing["fit"] = seconds_since(t0);
    log.info("fit: " + std::to_string(space.n_unknowns) + " control points, " +
             std::to_string(samples.block.size()) + " samples (" + std::to_string(samples.unlocated) +
             " unlocated), error rms " + std::to_string(quality.rms_error) + " max " +
             std::to_string(quality.max_error) + ", surface rms " + std::to_string(quality.surface_rms_error));
    log.info("fit: scaled Jacobian min " + std::to_string(quality.min_scaled_jacobian) + " mean " +
             std::to_string(quality.mean_scaled_jacobian) + ", non-positive " +
             std::to_string(quality.nonpositive_jacobian) + " of " + std::to_string(quality.eval_points) +
             " points in " + std::to_string(quality.blocks_with_nonpositive) + " blocks, C0 gap " +
             std::to_string(quality.c0_max_gap));
  } catch (const Error& e) {
    timing["total"] = seconds_since(t_all);
    stats["error"] = e.what();
    save_stats();
    log.info(std::string("error: ") + e.what());
    throw;
  }
  timing["total"] = seconds_since(t_all);
  save_stats();
  return stats;
}

}  // namespace pcp
