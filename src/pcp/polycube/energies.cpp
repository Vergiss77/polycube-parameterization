#include "pcp/polycube/energies.h"

#include <cmath>
#include <limits>

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include "pcp/geom/triangle.h"

namespace pcp {

namespace {

using Vec9 = Eigen::Matrix<double, 9, 1>;
using Vec12 = Eigen::Matrix<double, 12, 1>;
using Mat9 = Eigen::Matrix<double, 9, 9>;
using Mat12 = Eigen::Matrix<double, 12, 12>;

inline void add_face_gradient(const MatX3i& B, Index f, const Vec9& gf, Eigen::VectorXd& g) {
  for (int a = 0; a < 3; ++a)
    for (int k = 0; k < 3; ++k) g(3 * B(f, a) + k) += gf(3 * a + k);
}

inline void add_pair_gradient(const FacePair& p, const Vec12& gp, Eigen::VectorXd& g) {
  for (int a = 0; a < 4; ++a)
    for (int k = 0; k < 3; ++k) g(3 * p.nodes[static_cast<size_t>(a)] + k) += gp(3 * a + k);
}

inline std::array<int, 3> face_nodes(const MatX3i& B, Index f) {
  return {B(f, 0), B(f, 1), B(f, 2)};
}

Mat9 psd_part(const Mat9& H) {
  Eigen::SelfAdjointEigenSolver<Mat9> es(H);
  Vec9 ev = es.eigenvalues().cwiseMax(0.0);
  return es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();
}

}  // namespace

// ---------------------------------------------------------------------------------------
// Shared helpers

double polycube_error(const MatX3d& X, const MatX3i& B) {
  double l1 = 0.0, area = 0.0;
  Eigen::Vector3d x0, x1, x2;
  for (Index f = 0; f < B.rows(); ++f) {
    face_vertices(X, B, f, x0, x1, x2);
    Eigen::Vector3d a = tri::area_vector(x0, x1, x2);
    l1 += a.lpNorm<1>();
    area += a.norm();
  }
  return (l1 - area) / area;
}

double area_constraint(const MatX3d& X, const MatX3i& B, double target) {
  double area = 0.0;
  Eigen::Vector3d x0, x1, x2;
  for (Index f = 0; f < B.rows(); ++f) {
    face_vertices(X, B, f, x0, x1, x2);
    area += tri::area_vector(x0, x1, x2).norm();
  }
  return area - target;
}

void area_constraint_gradient(const MatX3d& X, const MatX3i& B, Eigen::VectorXd& g) {
  g.setZero(3 * X.rows());
  Eigen::Vector3d x0, x1, x2;
  for (Index f = 0; f < B.rows(); ++f) {
    face_vertices(X, B, f, x0, x1, x2);
    Eigen::Vector3d a = tri::area_vector(x0, x1, x2);
    double len = a.norm();
    if (len <= 0.0) continue;
    Vec9 gf = tri::area_vector_jacobian(x0, x1, x2).transpose() * (a / len);
    add_face_gradient(B, f, gf, g);
  }
}

bool pair_normals(const MatX3d& X, const FacePair& p, Eigen::Vector3d& n0, Eigen::Vector3d& n1,
                  Eigen::Matrix<double, 3, 12>& J0, Eigen::Matrix<double, 3, 12>& J1) {
  const Eigen::Vector3d xa = X.row(p.nodes[0]), xb = X.row(p.nodes[1]), xc = X.row(p.nodes[2]),
                        xd = X.row(p.nodes[3]);
  // Face 0 is (v0, v1, o0); face 1 is (v1, v0, o1) in its own orientation.
  Eigen::Vector3d a0 = tri::area_vector(xa, xb, xc), a1 = tri::area_vector(xb, xa, xd);
  double l0 = a0.norm(), l1 = a1.norm();
  if (!(l0 > 0.0) || !(l1 > 0.0)) return false;
  n0 = a0 / l0;
  n1 = a1 / l1;
  tri::Mat39 j0 = tri::unit_normal_jacobian(xa, xb, xc);
  tri::Mat39 j1 = tri::unit_normal_jacobian(xb, xa, xd);
  J0.setZero();
  J0.block<3, 9>(0, 0) = j0;
  J1.setZero();
  J1.block<3, 3>(0, 3) = j1.block<3, 3>(0, 0);
  J1.block<3, 3>(0, 0) = j1.block<3, 3>(0, 3);
  J1.block<3, 3>(0, 9) = j1.block<3, 3>(0, 6);
  return true;
}

// ---------------------------------------------------------------------------------------
// E1

L1NormalEnergy::L1NormalEnergy(const DeformMesh& mesh) : EnergyTerm(mesh) {
  for (int c = 0; c < 3; ++c) {
    Mat9 H = tri::area_vector_hessian(c);
    hess_pos_[static_cast<size_t>(c)] = psd_part(H);
    hess_neg_[static_cast<size_t>(c)] = psd_part(-H);
  }
  lagged_area_ = Eigen::VectorXd::Zero(mesh.B.rows());
}

void L1NormalEnergy::begin_iteration(const MatX3d& X) {
  Eigen::Vector3d x0, x1, x2;
  for (Index f = 0; f < mesh_.B.rows(); ++f) {
    face_vertices(X, mesh_.B, f, x0, x1, x2);
    lagged_area_(f) = tri::area_vector(x0, x1, x2).norm();
  }
}

double L1NormalEnergy::value(const MatX3d& X) const {
  double e = 0.0;
  Eigen::Vector3d x0, x1, x2;
  for (Index f = 0; f < mesh_.B.rows(); ++f) {
    face_vertices(X, mesh_.B, f, x0, x1, x2);
    Eigen::Vector3d a = tri::area_vector(x0, x1, x2);
    const double s = eps * lagged_area_(f) * lagged_area_(f);
    for (int c = 0; c < 3; ++c) e += std::sqrt(a(c) * a(c) + s);
  }
  return weight * e / mesh_.rest_area;
}

void L1NormalEnergy::add_gradient(const MatX3d& X, Eigen::VectorXd& g) const {
  const double w = weight / mesh_.rest_area;
  Eigen::Vector3d x0, x1, x2;
  for (Index f = 0; f < mesh_.B.rows(); ++f) {
    face_vertices(X, mesh_.B, f, x0, x1, x2);
    Eigen::Vector3d a = tri::area_vector(x0, x1, x2);
    tri::Mat39 J = tri::area_vector_jacobian(x0, x1, x2);
    const double s = eps * lagged_area_(f) * lagged_area_(f);
    Eigen::Vector3d coef;
    for (int c = 0; c < 3; ++c) coef(c) = a(c) / std::sqrt(a(c) * a(c) + s);
    add_face_gradient(mesh_.B, f, w * (J.transpose() * coef), g);
  }
}

void L1NormalEnergy::add_hessian(const MatX3d& X, BlockSparse& H) const {
  const double w = weight / mesh_.rest_area;
  Eigen::Vector3d x0, x1, x2;
  for (Index f = 0; f < mesh_.B.rows(); ++f) {
    face_vertices(X, mesh_.B, f, x0, x1, x2);
    Eigen::Vector3d a = tri::area_vector(x0, x1, x2);
    tri::Mat39 J = tri::area_vector_jacobian(x0, x1, x2);
    const double s = eps * lagged_area_(f) * lagged_area_(f);
    Mat9 He = Mat9::Zero();
    for (int c = 0; c < 3; ++c) {
      const double at = std::sqrt(a(c) * a(c) + s);
      const double ratio = a(c) / at;
      Vec9 gc = J.row(c).transpose();
      if (s > 0.0) He += (s / (at * at * at)) * gc * gc.transpose();
      if (ratio > 0.0)
        He += ratio * hess_pos_[static_cast<size_t>(c)];
      else if (ratio < 0.0)
        He += (-ratio) * hess_neg_[static_cast<size_t>(c)];
    }
    H.add_element<3>(face_nodes(mesh_.B, f), w * He);
  }
}

// ---------------------------------------------------------------------------------------
// ARAP

Eigen::VectorXd ArapEnergy::per_tet(const MatX3d& X) const {
  Eigen::VectorXd d(mesh_.T.rows());
  for (Index t = 0; t < mesh_.T.rows(); ++t) {
    Eigen::Matrix3d F = deformation_gradient(mesh_.grads, X, mesh_.T, t);
    d(t) = 0.5 * (F - closest_rotation(F)).squaredNorm();
  }
  return d;
}

double ArapEnergy::value(const MatX3d& X) const {
  double e = 0.0;
  for (Index t = 0; t < mesh_.T.rows(); ++t) {
    Eigen::Matrix3d F = deformation_gradient(mesh_.grads, X, mesh_.T, t);
    e += mesh_.grads.rest_volume(t) * 0.5 * (F - closest_rotation(F)).squaredNorm();
  }
  return weight * e / mesh_.rest_volume;
}

void ArapEnergy::add_gradient(const MatX3d& X, Eigen::VectorXd& g) const {
  for (Index t = 0; t < mesh_.T.rows(); ++t) {
    const auto& G = mesh_.grads.G[static_cast<size_t>(t)];
    Eigen::Matrix3d F = deformation_gradient(mesh_.grads, X, mesh_.T, t);
    Eigen::Matrix3d D = F - closest_rotation(F);
    const double c = weight * mesh_.grads.rest_volume(t) / mesh_.rest_volume;
    for (int i = 0; i < 4; ++i) {
      Eigen::Vector3d gi = c * (D * G.col(i));
      for (int k = 0; k < 3; ++k) g(3 * mesh_.T(t, i) + k) += gi(k);
    }
  }
}

void ArapEnergy::add_hessian(const MatX3d& /*X*/, BlockSparse& H) const {
  for (Index t = 0; t < mesh_.T.rows(); ++t) {
    const auto& G = mesh_.grads.G[static_cast<size_t>(t)];
    const double c = weight * mesh_.grads.rest_volume(t) / mesh_.rest_volume;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        H.add_diag(mesh_.T(t, i), mesh_.T(t, j), c * G.col(i).dot(G.col(j)));
  }
}

// ---------------------------------------------------------------------------------------
// E_eta

void ComplexityEnergy::begin_iteration(const MatX3d& X) {
  const Index nf = mesh_.B.rows();
  Eigen::VectorXd A(nf);
  Eigen::Vector3d x0, x1, x2;
  for (Index f = 0; f < nf; ++f) {
    face_vertices(X, mesh_.B, f, x0, x1, x2);
    A(f) = tri::area_vector(x0, x1, x2).norm();
  }
  Eigen::VectorXd S = Eigen::VectorXd::Zero(nf);
  for (Index f = 0; f < nf; ++f)
    for (int k = 0; k < 3; ++k) S(f) += A(mesh_.topo.face_nbrs[static_cast<size_t>(f)][static_cast<size_t>(k)]);
  edge_area_.resize(static_cast<Index>(mesh_.topo.pairs.size()));
  for (size_t e = 0; e < mesh_.topo.pairs.size(); ++e) {
    const FacePair& p = mesh_.topo.pairs[e];
    // A(e_ij) = gamma_ij A_i + gamma_ji A_j with gamma_ij = A_j / sum_{k in N(i)} A_k.
    double ai = A(p.f0), aj = A(p.f1);
    double v = 0.0;
    if (S(p.f0) > 0.0) v += aj / S(p.f0) * ai;
    if (S(p.f1) > 0.0) v += ai / S(p.f1) * aj;
    edge_area_(static_cast<Index>(e)) = v;
  }
}

double ComplexityEnergy::value(const MatX3d& X) const {
  double e = 0.0;
  Eigen::Vector3d n0, n1;
  Eigen::Matrix<double, 3, 12> J0, J1;
  for (size_t i = 0; i < mesh_.topo.pairs.size(); ++i) {
    if (!pair_normals(X, mesh_.topo.pairs[i], n0, n1, J0, J1)) continue;
    e += edge_area_(static_cast<Index>(i)) * (n0 - n1).squaredNorm();
  }
  return weight * e / mesh_.rest_area;
}

void ComplexityEnergy::add_gradient(const MatX3d& X, Eigen::VectorXd& g) const {
  Eigen::Vector3d n0, n1;
  Eigen::Matrix<double, 3, 12> J0, J1;
  for (size_t i = 0; i < mesh_.topo.pairs.size(); ++i) {
    const FacePair& p = mesh_.topo.pairs[i];
    if (!pair_normals(X, p, n0, n1, J0, J1)) continue;
    double w = 2.0 * weight * edge_area_(static_cast<Index>(i)) / mesh_.rest_area;
    add_pair_gradient(p, w * ((J0 - J1).transpose() * (n0 - n1)), g);
  }
}

void ComplexityEnergy::add_hessian(const MatX3d& X, BlockSparse& H) const {
  Eigen::Vector3d n0, n1;
  Eigen::Matrix<double, 3, 12> J0, J1;
  for (size_t i = 0; i < mesh_.topo.pairs.size(); ++i) {
    const FacePair& p = mesh_.topo.pairs[i];
    if (!pair_normals(X, p, n0, n1, J0, J1)) continue;
    double w = 2.0 * weight * edge_area_(static_cast<Index>(i)) / mesh_.rest_area;
    Eigen::Matrix<double, 3, 12> J = J0 - J1;
    Mat12 He = w * J.transpose() * J;
    H.add_element<4>(p.nodes, He);
  }
}

// ---------------------------------------------------------------------------------------
// E_e

double FeatureEdgeEnergy::value(const MatX3d& X) const {
  double e = 0.0;
  Eigen::Vector3d n0, n1;
  Eigen::Matrix<double, 3, 12> J0, J1;
  for (const Edge& ed : edges_) {
    if (!pair_normals(X, mesh_.topo.pairs[static_cast<size_t>(ed.pair)], n0, n1, J0, J1)) continue;
    double r = n0.dot(n1) - ed.target;
    e += r * r;
  }
  return weight * e;
}

void FeatureEdgeEnergy::add_gradient(const MatX3d& X, Eigen::VectorXd& g) const {
  Eigen::Vector3d n0, n1;
  Eigen::Matrix<double, 3, 12> J0, J1;
  for (const Edge& ed : edges_) {
    const FacePair& p = mesh_.topo.pairs[static_cast<size_t>(ed.pair)];
    if (!pair_normals(X, p, n0, n1, J0, J1)) continue;
    double r = n0.dot(n1) - ed.target;
    Vec12 grad = J0.transpose() * n1 + J1.transpose() * n0;
    add_pair_gradient(p, 2.0 * weight * r * grad, g);
  }
}

void FeatureEdgeEnergy::add_hessian(const MatX3d& X, BlockSparse& H) const {
  Eigen::Vector3d n0, n1;
  Eigen::Matrix<double, 3, 12> J0, J1;
  for (const Edge& ed : edges_) {
    const FacePair& p = mesh_.topo.pairs[static_cast<size_t>(ed.pair)];
    if (!pair_normals(X, p, n0, n1, J0, J1)) continue;
    Vec12 grad = J0.transpose() * n1 + J1.transpose() * n0;
    Mat12 He = 2.0 * weight * grad * grad.transpose();
    H.add_element<4>(p.nodes, He);
  }
}

// ---------------------------------------------------------------------------------------
// Barrier

double FoldoverBarrier::min_slack(const MatX3d& X) const {
  double mn = std::numeric_limits<double>::infinity();
  Eigen::Vector3d n0, n1;
  Eigen::Matrix<double, 3, 12> J0, J1;
  for (const FacePair& p : mesh_.topo.pairs) {
    if (!pair_normals(X, p, n0, n1, J0, J1)) return 0.0;
    mn = std::min(mn, n0.dot(n1) + 1.0);
  }
  return mn;
}

double FoldoverBarrier::value(const MatX3d& X) const {
  double e = 0.0;
  Eigen::Vector3d n0, n1;
  Eigen::Matrix<double, 3, 12> J0, J1;
  for (const FacePair& p : mesh_.topo.pairs) {
    if (!pair_normals(X, p, n0, n1, J0, J1)) return std::numeric_limits<double>::infinity();
    double s = n0.dot(n1) + 1.0;
    if (!(s > 0.0)) return std::numeric_limits<double>::infinity();
    e -= std::log(s);
  }
  return weight * e;
}

void FoldoverBarrier::add_gradient(const MatX3d& X, Eigen::VectorXd& g) const {
  Eigen::Vector3d n0, n1;
  Eigen::Matrix<double, 3, 12> J0, J1;
  for (const FacePair& p : mesh_.topo.pairs) {
    if (!pair_normals(X, p, n0, n1, J0, J1)) continue;
    double s = n0.dot(n1) + 1.0;
    if (!(s > 0.0)) continue;
    Vec12 grad = J0.transpose() * n1 + J1.transpose() * n0;
    add_pair_gradient(p, (-weight / s) * grad, g);
  }
}

void FoldoverBarrier::add_hessian(const MatX3d& X, BlockSparse& H) const {
  Eigen::Vector3d n0, n1;
  Eigen::Matrix<double, 3, 12> J0, J1;
  for (const FacePair& p : mesh_.topo.pairs) {
    if (!pair_normals(X, p, n0, n1, J0, J1)) continue;
    double s = n0.dot(n1) + 1.0;
    if (!(s > 0.0)) continue;
    Vec12 grad = J0.transpose() * n1 + J1.transpose() * n0;
    Mat12 He = (weight / (s * s)) * grad * grad.transpose();
    H.add_element<4>(p.nodes, He);
  }
}

// ---------------------------------------------------------------------------------------
// Tet inversion barrier

namespace {
// Volume of tet t and its gradient w.r.t. the 12 stacked node coordinates.
double tet_volume_gradient(const MatX3d& X, const MatX4i& T, Index t, Vec12* grad) {
  const Eigen::Vector3d a = X.row(T(t, 0)), b = X.row(T(t, 1)), c = X.row(T(t, 2)), d = X.row(T(t, 3));
  const double v = (b - a).dot((c - a).cross(d - a)) / 6.0;
  if (grad) {
    Eigen::Vector3d gb = (c - a).cross(d - a) / 6.0, gc = (d - a).cross(b - a) / 6.0,
                    gd = (b - a).cross(c - a) / 6.0;
    grad->segment<3>(0) = -(gb + gc + gd);
    grad->segment<3>(3) = gb;
    grad->segment<3>(6) = gc;
    grad->segment<3>(9) = gd;
  }
  return v;
}
}  // namespace

Index TetInversionBarrier::inverted(const MatX3d& X) const {
  Index n = 0;
  for (Index t = 0; t < mesh_.T.rows(); ++t)
    if (!(tet_volume_gradient(X, mesh_.T, t, nullptr) > 0.0)) ++n;
  return n;
}

double TetInversionBarrier::value(const MatX3d& X) const {
  double e = 0.0;
  for (Index t = 0; t < mesh_.T.rows(); ++t) {
    const double V = mesh_.grads.rest_volume(t);
    const double v = tet_volume_gradient(X, mesh_.T, t, nullptr);
    if (!(v > 0.0)) return std::numeric_limits<double>::infinity();
    e -= V * std::log(v / V);
  }
  return weight * e / mesh_.rest_volume;
}

void TetInversionBarrier::add_gradient(const MatX3d& X, Eigen::VectorXd& g) const {
  Vec12 gv;
  for (Index t = 0; t < mesh_.T.rows(); ++t) {
    const double V = mesh_.grads.rest_volume(t);
    const double v = tet_volume_gradient(X, mesh_.T, t, &gv);
    if (!(v > 0.0)) continue;
    const double c = -weight * V / (v * mesh_.rest_volume);
    for (int a = 0; a < 4; ++a)
      for (int k = 0; k < 3; ++k) g(3 * mesh_.T(t, a) + k) += c * gv(3 * a + k);
  }
}

void TetInversionBarrier::add_hessian(const MatX3d& X, BlockSparse& H) const {
  Vec12 gv;
  for (Index t = 0; t < mesh_.T.rows(); ++t) {
    const double V = mesh_.grads.rest_volume(t);
    const double v = tet_volume_gradient(X, mesh_.T, t, &gv);
    if (!(v > 0.0)) continue;
    const Mat12 He = (weight * V / (v * v * mesh_.rest_volume)) * gv * gv.transpose();
    H.add_element<4>({mesh_.T(t, 0), mesh_.T(t, 1), mesh_.T(t, 2), mesh_.T(t, 3)}, He);
  }
}

}  // namespace pcp
