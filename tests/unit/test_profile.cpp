#include <doctest/doctest.h>
#include <cstdlib>
#include <sstream>

#include "helpers.h"
#include "pcp/geom/normalize.h"
#include "pcp/io/off.h"
#include "pcp/polycube/l1_solver.h"
#include "pcp/tet/tetrahedralize.h"

using namespace pcp;

// Manual profiling helper: a few iterations on CLEVIS. Run with -tc=profile_l1.
TEST_CASE("profile_l1" * doctest::test_suite("profile") * doctest::skip()) {
  SurfaceMesh s = read_off(testutil::models_dir() / "CLEVIS-D64EEBFA.off");
  TetMesh t = tetrahedralize(s, TetParams{});
  Normalization n = normalization_for(t.X, t.B);
  MatX3d rest = n.apply(t.X);
  DeformMesh mesh = DeformMesh::build(rest, t.T, t.B);
  L1Params p;
  p.max_stages = 1;
  p.max_iters_per_stage = 3;
  p.global_rotation = false;
  if (const char* b = std::getenv("PCP_BETA")) p.beta = std::atof(b);
  L1Result r = solve_l1_polycube(rest, mesh, p, {}, [](const std::string& s) { MESSAGE(s); });
  MESSAGE("iters " << r.iterations << " total " << r.seconds << " " << r.timing.to_json().dump());
}

#include <chrono>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include "pcp/polycube/energies.h"

TEST_CASE("profile_linear" * doctest::test_suite("profile") * doctest::skip()) {
  SurfaceMesh s = read_off(testutil::models_dir() / "CLEVIS-D64EEBFA.off");
  TetMesh t = tetrahedralize(s, TetParams{});
  Normalization n = normalization_for(t.X, t.B);
  MatX3d rest = n.apply(t.X);
  DeformMesh mesh = DeformMesh::build(rest, t.T, t.B);
  MatX3d X = rest;
  L1NormalEnergy e1(mesh);
  ArapEnergy arap(mesh);
  ComplexityEnergy eta(mesh);
  for (double alpha : {0.1, 100.0}) {
    e1.weight = alpha;
    e1.eps = alpha > 1 ? 0.01 : 1.0;
    e1.begin_iteration(X);
    eta.begin_iteration(X);
    BlockSparse H = mesh.make_hessian();
    e1.add_hessian(X, H);
    arap.add_hessian(X, H);
    eta.add_hessian(X, H);
    Eigen::SparseMatrix<double> A = H.matrix();
    double dm = 0;
    for (Index i = 0; i < A.rows(); ++i) dm += A.coeff(i, i);
    dm /= A.rows();
    for (Index i = 0; i < A.rows(); ++i) A.coeffRef(i, i) += 1e-8 * dm;
    Eigen::VectorXd b = Eigen::VectorXd::Random(A.rows());
    auto now = [] { return std::chrono::steady_clock::now(); };
    auto sec = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };
    MESSAGE("alpha " << alpha << " dofs " << A.rows() << " nnz " << A.nonZeros());
    {
      Eigen::SparseMatrix<double, Eigen::RowMajor> Ar = A;
      Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
                               Eigen::IncompleteCholesky<double>> cg;
      auto t0 = now();
      cg.compute(A);
      cg.setTolerance(1e-8);
      cg.setMaxIterations(5000);
      auto t1 = now();
      Eigen::VectorXd x = cg.solve(b);
      auto t2 = now();
      MESSAGE("  CG+IC: setup " << sec(t0, t1) << " solve " << sec(t1, t2) << " iters " << cg.iterations()
                                << " err " << cg.error() << " info " << int(cg.info()));
    }
    {
      Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
      cg.setTolerance(1e-8);
      cg.setMaxIterations(20000);
      auto t0 = now();
      cg.compute(A);
      Eigen::VectorXd x = cg.solve(b);
      auto t1 = now();
      MESSAGE("  CG+diag: " << sec(t0, t1) << " iters " << cg.iterations() << " err " << cg.error());
    }
    if (alpha < 1) {
      Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
      auto t0 = now();
      llt.analyzePattern(A);
      auto t1 = now();
      llt.factorize(A);
      auto t2 = now();
      MESSAGE("  LLT analyze " << sec(t0, t1) << " factorize " << sec(t1, t2) << " nnz(L) "
                               << llt.matrixL().nestedExpression().nonZeros());
    }
  }
}

#include "pcp/io/off.h"
#include "pcp/polycube/rotation.h"

// Full run on PCP_MODEL with optional PCP_BETA; writes the deformed boundary to PCP_OUT.off.
TEST_CASE("profile_full" * doctest::test_suite("profile") * doctest::skip()) {
  const char* model = std::getenv("PCP_MODEL");
  REQUIRE(model);
  SurfaceMesh s = read_off(testutil::models_dir() / model);
  TetMesh t = tetrahedralize(s, TetParams{});
  Normalization n = normalization_for(t.X, t.B);
  MatX3d rest = n.apply(t.X);
  DeformMesh mesh = DeformMesh::build(rest, t.T, t.B);
  L1Params p;
  if (const char* b = std::getenv("PCP_BETA")) p.beta = std::atof(b);
  L1Result r = solve_l1_polycube(rest, mesh, p, {}, [](const std::string& s) { MESSAGE(s); });
  MESSAGE("error " << r.error << " arap " << r.arap_mean << " max " << r.arap_max << " iters "
                   << r.iterations << " time " << r.seconds << " rot "
                   << rotation_angle(r.rotation) * 57.29578 << " stalls " << r.stalls.size() << " "
                   << r.timing.to_json().dump());
  if (const char* out = std::getenv("PCP_OUT")) {
    SurfaceMesh o;
    o.V = n.inverse(r.X);
    o.F = t.B;
    write_off(std::string(out) + ".off", o);
  }
}

// Sphere with the iterative solver forced, to watch CG iteration counts across the schedule.
TEST_CASE("profile_sphere_cg" * doctest::test_suite("profile") * doctest::skip()) {
  SurfaceMesh s = testutil::make_icosphere(4);
  TetMesh t = tetrahedralize(s, TetParams{});
  Normalization n = normalization_for(t.X, t.B);
  MatX3d rest = n.apply(t.X);
  DeformMesh mesh = DeformMesh::build(rest, t.T, t.B);
  L1Params p;
  p.linear_solver = std::getenv("PCP_SOLVER") ? std::getenv("PCP_SOLVER") : "cg";
  if (const char* m = std::getenv("PCP_CGMAX")) p.cg_max_iterations = std::atoi(m);
  p.log_iterations = std::getenv("PCP_ITERLOG") != nullptr;
  L1Result r = solve_l1_polycube(rest, mesh, p, {}, [](const std::string& s) { MESSAGE(s); });
  MESSAGE("nodes " << t.X.rows() << " error " << r.error << " iters " << r.iterations << " time "
                   << r.seconds << " " << r.timing.to_json().dump());
}

#include "pcp/io/text.h"
#include "pcp/io/tmsh.h"
#include "pcp/structure/final_map.h"
#include "pcp/structure/labeling.h"
#include "pcp/structure/patch_graph.h"

// Final-map experiments on a saved pipeline output directory PCP_DIR.
TEST_CASE("profile_final_map" * doctest::test_suite("profile") * doctest::skip()) {
  const char* dir = std::getenv("PCP_DIR");
  REQUIRE(dir);
  TetMesh t = read_tmsh(std::string(dir) + "/01_tet.tmsh");
  TetMesh l1 = read_tmsh(std::string(dir) + "/02_l1.tmsh");
  Normalization n = normalization_for(t.X, t.B);
  MatX3d rest = n.apply(t.X);
  MatX3d X = (l1.X.rowwise() - n.center) * n.scale;  // same frame as the solver output
  BoundaryTopology topo = build_boundary_topology(t.B);
  Eigen::VectorXi labels = initial_labels(X, t.B);
  cleanup_labels(X, t.B, topo, labels);
  PatchGraph g = build_patch_graph(X, t.B, topo, labels);
  for (bool lap : {true})
    for (double tol : {1e-7}) {
      FinalMapParams p;
      p.laplacian_init = lap;
      p.rel_tolerance = tol;
      p.max_iterations = 1000;
      FinalMapResult r = compute_final_map(rest, X, t.T, t.B, g, g.height, p);
      MESSAGE("laplacian_init " << lap << " tol " << tol << ": iters " << r.iterations << " inverted "
                                << r.inverted_tets << " arap mean " << r.arap_mean << " max "
                                << r.arap_max << " before untangle " << r.inverted_before_untangle
                                << " time " << r.seconds);
    }
}

#include <map>
#include <set>
#include "pcp/tet/tetrahedralize.h"
TEST_CASE("profile_inverted" * doctest::test_suite("profile") * doctest::skip()) {
  const char* dir = std::getenv("PCP_DIR");
  REQUIRE(dir);
  TetMesh t = read_tmsh(std::string(dir) + "/01_tet.tmsh");
  TetMesh pc = read_tmsh(std::string(dir) + "/04_polycube.tmsh");
  nlohmann::json g = nlohmann::json::parse(text::read_file(std::string(dir) + "/03_patch_graph.json"));
  std::vector<std::set<int>> vpatch(t.X.rows());
  for (const auto& p : g["patches"])
    for (int f : p["faces"])
      for (int k = 0; k < 3; ++k) vpatch[t.B(f, k)].insert(p["id"].get<int>());
  int all_same = 0, n = 0;
  for (Index i = 0; i < t.T.rows(); ++i) {
    // Patches shared by all 4 nodes.
    std::set<int> common = vpatch[t.T(i, 0)];
    for (int k = 1; k < 4; ++k) {
      std::set<int> c2;
      for (int p : common) if (vpatch[t.T(i, k)].count(p)) c2.insert(p);
      common = c2;
    }
    bool inv = tet_volume(pc.X, pc.T, i) <= 0.0;
    if (!common.empty()) ++all_same;
    if (inv) {
      ++n;
      MESSAGE("inverted tet " << i << " volume " << tet_volume(pc.X, pc.T, i) << " boundary nodes "
              << (int)!vpatch[t.T(i,0)].empty() + !vpatch[t.T(i,1)].empty() + !vpatch[t.T(i,2)].empty() + !vpatch[t.T(i,3)].empty()
              << " common patches " << common.size());
    }
  }
  MESSAGE("tets with all 4 nodes on one patch: " << all_same << ", inverted: " << n);
}

// Lists the same-axis adjacent patch pairs of a saved l1 result (PCP_DIR).
TEST_CASE("profile_structure" * doctest::test_suite("profile") * doctest::skip()) {
  const char* dir = std::getenv("PCP_DIR");
  REQUIRE(dir);
  TetMesh l1 = read_tmsh(std::string(dir) + "/02_l1.tmsh");
  BoundaryTopology topo = build_boundary_topology(l1.B);
  Eigen::VectorXi labels = initial_labels(l1.X, l1.B);
  CleanupReport rep = cleanup_labels(l1.X, l1.B, topo, labels);
  PatchGraph g = build_patch_graph(l1.X, l1.B, topo, labels);
  static const char* nm[6] = {"+X", "-X", "+Y", "-Y", "+Z", "-Z"};
  auto centroid = [&](int p) {
    Eigen::RowVector3d c = Eigen::RowVector3d::Zero();
    for (int f : g.patches.faces[p]) c += (l1.X.row(l1.B(f, 0)) + l1.X.row(l1.B(f, 1)) + l1.X.row(l1.B(f, 2))) / 3.0;
    return (c / static_cast<double>(g.patches.faces[p].size())).eval();
  };
  for (size_t p = 0; p < g.patches.label.size(); ++p)
    for (int q : g.patches.neighbors[p])
      if (static_cast<size_t>(q) > p && label_axis(g.patches.label[p]) == label_axis(g.patches.label[q]))
        MESSAGE(p << " " << nm[g.patches.label[p]] << " (" << g.patches.faces[p].size() << " faces, area "
                  << g.area[p] << ") at " << centroid(static_cast<int>(p)) << "  |  " << q << " "
                  << nm[g.patches.label[q]] << " (" << g.patches.faces[q].size() << " faces) at "
                  << centroid(q));
}
TEST_CASE("profile_structure_check" * doctest::test_suite("profile") * doctest::skip()) {
  TetMesh l1 = read_tmsh(std::string(std::getenv("PCP_DIR")) + "/02_l1.tmsh");
  BoundaryTopology topo = build_boundary_topology(l1.B);
  Eigen::VectorXi labels = initial_labels(l1.X, l1.B);
  CleanupReport rep = cleanup_labels(l1.X, l1.B, topo, labels);
  PatchGraph g = build_patch_graph(l1.X, l1.B, topo, labels);
  TopologyCheck tc = check_polycube_topology(l1.B, g);
  MESSAGE(rep.to_json().dump() << " patches " << g.patches.label.size() << " corners " << g.corners.size()
                               << " " << tc.to_json().dump());
}
TEST_CASE("profile_heights" * doctest::test_suite("profile") * doctest::skip()) {
  const std::string dir = std::getenv("PCP_DIR");
  TetMesh t = read_tmsh(dir + "/01_tet.tmsh");
  TetMesh l1 = read_tmsh(dir + "/02_l1.tmsh");
  Normalization n = normalization_for(t.X, t.B);
  MatX3d rest = n.apply(t.X);
  MatX3d X = (l1.X.rowwise() - n.center) * n.scale;
  BoundaryTopology topo = build_boundary_topology(t.B);
  Eigen::VectorXi labels = initial_labels(X, t.B);
  cleanup_labels(X, t.B, topo, labels);
  PatchGraph g = build_patch_graph(X, t.B, topo, labels);
  FinalMapParams p;
  p.untangle = false;
  FinalMapResult r = compute_final_map(rest, X, t.T, t.B, g, g.height, p);
  // Chains whose lower/upper patch order flipped between the l1 result and the final map.
  int flipped = 0, collapsed = 0;
  for (const Chain& c : g.chains) {
    if (c.lower_patch < 0 || c.upper_patch < 0 || c.lower_patch == c.upper_patch) continue;
    double before = g.height[c.upper_patch] - g.height[c.lower_patch];
    double after = r.heights[c.upper_patch] - r.heights[c.lower_patch];
    if (after <= 0) ++flipped;
    else if (after < 0.2 * before) ++collapsed;
    if (after < 0.5 * before)
      MESSAGE("chain axis " << c.axis << " length before " << before / n.scale << " after " << after / n.scale);
  }
  std::map<std::string, int> where;
  for (Index i = 0; i < t.T.rows(); ++i) {
    if (tet_volume(r.X, t.T, i) > 0.0) continue;
    Eigen::RowVector3d c = (t.X.row(t.T(i, 0)) + t.X.row(t.T(i, 1)) + t.X.row(t.T(i, 2)) + t.X.row(t.T(i, 3))) / 4.0;
    double rad = std::hypot(c(0), c(1));
    int nb = 0;
    for (int k = 0; k < 4; ++k) nb += t.T(i, k) < 28326 ? 1 : 0;  // input surface vertices come first
    char key[96];
    std::snprintf(key, sizeof(key), "r %2.0f-%2.0f  z %3.1f  boundary nodes %d", std::floor(rad / 5) * 5,
                  std::floor(rad / 5) * 5 + 5, std::round(c(2) * 2) / 2, nb);
    ++where[key];
  }
  for (auto& [k, v] : where) MESSAGE(v << "x  " << k);
  MESSAGE("inverted " << r.inverted_tets << " flipped chains " << flipped << " collapsed (<20%) " << collapsed);
}

#include "pcp/geom/boundary_topology.h"
#include "pcp/structure/height_simplify.h"

// Height simplification on a finished run: PCP_DIR=<output dir with 01_tet.tmsh,
// 03_patch_graph.json and 04_polycube.tmsh>, optional PCP_STRETCH=<s list, comma separated> and
// PCP_REMAP=1 to also re-map the volume for each s (as the pipeline's simplify stage does).
TEST_CASE("profile_simplify" * doctest::test_suite("profile") * doctest::skip()) {
  const char* dir_c = std::getenv("PCP_DIR");
  REQUIRE(dir_c != nullptr);
  const std::filesystem::path dir(dir_c);
  TetMesh tet = read_tmsh(dir / "01_tet.tmsh");
  TetMesh pc = read_tmsh(dir / "04_polycube.tmsh");
  const auto pj = nlohmann::json::parse(text::read_file(dir / "03_patch_graph.json"));
  Eigen::VectorXi labels = Eigen::VectorXi::Constant(pc.B.rows(), -1);
  static const std::map<std::string, int> names = {{"+X", 0}, {"-X", 1}, {"+Y", 2}, {"-Y", 3}, {"+Z", 4}, {"-Z", 5}};
  for (const auto& p : pj["patches"])
    for (int f : p["faces"].get<std::vector<int>>()) labels(f) = names.at(p["label"].get<std::string>());
  REQUIRE(labels.minCoeff() >= 0);
  const Normalization norm = normalization_for(tet.X, tet.B);
  const MatX3d rest = norm.apply(tet.X), X = norm.apply(pc.X);
  BoundaryTopology topo = build_boundary_topology(pc.B);
  PatchGraph g = build_patch_graph(X, pc.B, topo, labels);
  HeightGraph chains_only = build_height_graph(g, g.height);
  HeightGraph hg = chains_only;
  add_separation_edges(hg, g, g.height, X, pc.B);
  add_width_edges(hg, g, g.height, X, pc.B);
  MESSAGE("width edges " << hg.width_edges);
  int per_axis[3] = {0, 0, 0};
  for (const HeightGroup& grp : hg.groups) ++per_axis[grp.axis];
  MESSAGE("patches " << g.patches.label.size() << " groups " << hg.groups.size() << " (x "
                     << per_axis[0] << " y " << per_axis[1] << " z " << per_axis[2] << ") chain edges "
                     << chains_only.edges.size() << " separation edges " << hg.separation_edges);
  const bool remap = std::getenv("PCP_REMAP") != nullptr;
  std::string list = std::getenv("PCP_STRETCH") ? std::getenv("PCP_STRETCH") : "1.5,2,3";
  std::stringstream ss(list);
  for (std::string item; std::getline(ss, item, ',');) {
    SimplifyParams sp;
    sp.max_stretch = std::atof(item.c_str());
    SimplifyResult a = simplify_heights(chains_only, sp.max_stretch), r = simplify_heights(hg, sp.max_stretch);
    std::string extra;
    if (remap) {
      FinalMapParams fp;
      fp.untangle = std::getenv("PCP_NOUNTANGLE") == nullptr;
      SimplifiedPolycube out = simplify_and_remap(rest, X, tet.T, pc.B, g, g.height, sp, fp);
      extra = " | arap " + std::to_string(out.map.arap_mean) + " / " + std::to_string(out.map.arap_max) +
              ", inverted " + std::to_string(out.map.inverted_tets) + " (before untangle " +
              std::to_string(out.map.inverted_before_untangle) + "), topology " +
              (out.topology.ok ? "ok" : "FAILED") + ", map " + std::to_string(out.map.seconds) + " s";
    }
    if (std::getenv("PCP_MERGES"))
      for (const MergeRecord& m : r.merges) {
        const HeightGroup &ga = hg.groups[size_t(m.a)], &gb = hg.groups[size_t(m.b)];
        auto oh = [&](double h) { return h / norm.scale + norm.center(ga.axis); };
        MESSAGE("  merge axis " << ga.axis << " groups " << m.a << " (h " << oh(ga.height) << ", "
                                << ga.patches.size() << " patches, label " << g.patches.label[size_t(ga.patches[0])]
                                << ") + " << m.b << " (h " << oh(gb.height) << ", label "
                                << g.patches.label[size_t(gb.patches[0])] << ") -> " << oh(m.height)
                                << " cost " << m.cost);
      }
    if (const char* ks = std::getenv("PCP_PREFIX")) {
      std::stringstream kk(ks);
      for (std::string it2; std::getline(kk, it2, ',');) {
        const size_t k = std::min(r.merges.size(), size_t(std::atoi(it2.c_str())));
        HeightState st = HeightState::initial(hg);
        for (size_t m = 0; m < k; ++m) {
          for (auto& x : st.rep) if (x == r.merges[m].b) x = r.merges[m].a;
          st.height[size_t(r.merges[m].a)] = r.merges[m].height;
        }
        std::vector<double> ph(g.patches.label.size());
        for (size_t gi = 0; gi < hg.groups.size(); ++gi)
          for (int p : hg.groups[gi].patches) ph[size_t(p)] = st.height[size_t(st.rep[gi])];
        FinalMapParams fp;
        fp.heights_fixed = true;
        fp.untangle = false;
        FinalMapResult fm = compute_final_map(rest, X, tet.T, pc.B, g, ph, fp);
        int cuts = 0;
        for (size_t gi = 0; gi < hg.groups.size(); ++gi) cuts += st.rep[gi] == int(gi);
        MESSAGE("  prefix " << k << " (last cost " << (k ? r.merges[k - 1].cost : 0.0) << "): cuts " << cuts
                            << ", arap " << fm.arap_mean << " / " << fm.arap_max << ", inverted " << fm.inverted_tets);
      }
    }
    MESSAGE("s " << sp.max_stretch << ": cuts " << r.cuts_total_before() << " -> " << r.cuts_total_after()
                 << " (x " << r.cuts_after[0] << " y " << r.cuts_after[1] << " z " << r.cuts_after[2]
                 << "), merges " << r.merges.size() << ", ratio [" << r.min_ratio << ", " << r.max_ratio
                 << "]; chains only: " << a.cuts_total_after() << extra);
  }
}
