#include <doctest/doctest.h>
#include <cstdlib>

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
