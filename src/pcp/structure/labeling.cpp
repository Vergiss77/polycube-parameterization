#include "pcp/structure/labeling.h"

#include <algorithm>
#include <map>

#include <Eigen/Geometry>

namespace pcp {

nlohmann::json CleanupReport::to_json() const {
  return {{"zigzag_faces", zigzag_faces},
          {"degree1_merged", degree1_merged},
          {"degree2_merged", degree2_merged},
          {"conflict_merged", conflict_merged}};
}

Eigen::VectorXi initial_labels(const MatX3d& X, const MatX3i& B) {
  Eigen::VectorXi L(B.rows());
  for (Index f = 0; f < B.rows(); ++f) {
    Eigen::Vector3d x0 = X.row(B(f, 0)), x1 = X.row(B(f, 1)), x2 = X.row(B(f, 2));
    Eigen::Vector3d n = (x1 - x0).cross(x2 - x0);
    int axis = 0;
    n.cwiseAbs().maxCoeff(&axis);
    L(f) = make_label(axis, n(axis) >= 0.0 ? 1 : -1);
  }
  return L;
}

Index remove_zigzags(const BoundaryTopology& topo, Eigen::VectorXi& labels) {
  Index total = 0;
  for (int pass = 0; pass < 100; ++pass) {
    Index changed = 0;
    for (Index f = 0; f < labels.size(); ++f) {
      const auto& nb = topo.face_nbrs[static_cast<size_t>(f)];
      int a = labels(nb[0]), b = labels(nb[1]), c = labels(nb[2]);
      int me = labels(f);
      int target = -1;
      if (a == b && a != me) target = a;
      else if (a == c && a != me) target = a;
      else if (b == c && b != me) target = b;
      if (target >= 0) {
        labels(f) = target;
        ++changed;
      }
    }
    total += changed;
    if (changed == 0) break;
  }
  return total;
}

Patches compute_patches(const BoundaryTopology& topo, const Eigen::VectorXi& labels) {
  const Index nf = labels.size();
  Patches P;
  P.face_patch = Eigen::VectorXi::Constant(nf, -1);
  std::vector<int> stack;
  for (Index s = 0; s < nf; ++s) {
    if (P.face_patch(s) >= 0) continue;
    const int id = static_cast<int>(P.label.size());
    P.label.push_back(labels(s));
    P.faces.emplace_back();
    P.face_patch(s) = id;
    stack.assign(1, static_cast<int>(s));
    while (!stack.empty()) {
      int f = stack.back();
      stack.pop_back();
      P.faces.back().push_back(f);
      for (int g : topo.face_nbrs[static_cast<size_t>(f)])
        if (P.face_patch(g) < 0 && labels(g) == labels(s)) {
          P.face_patch(g) = id;
          stack.push_back(g);
        }
    }
    std::sort(P.faces.back().begin(), P.faces.back().end());
  }
  P.neighbors.assign(P.label.size(), {});
  for (const FacePair& p : topo.pairs) {
    int a = P.face_patch(p.f0), b = P.face_patch(p.f1);
    if (a != b) {
      P.neighbors[static_cast<size_t>(a)].push_back(b);
      P.neighbors[static_cast<size_t>(b)].push_back(a);
    }
  }
  for (auto& n : P.neighbors) {
    std::sort(n.begin(), n.end());
    n.erase(std::unique(n.begin(), n.end()), n.end());
  }
  return P;
}

CleanupReport cleanup_labels(const MatX3d& X, const MatX3i& B, const BoundaryTopology& topo,
                             Eigen::VectorXi& labels) {
  CleanupReport rep;
  rep.zigzag_faces = remove_zigzags(topo, labels);
  for (int round = 0; round < 10000; ++round) {
    Patches P = compute_patches(topo, labels);
    // Process the smallest offending patch first (deterministic, and small spurious patches are
    // the ones expected to go). Offending: fewer than 3 neighbors, or a larger neighbor on the
    // same axis.
    auto conflicting = [&](size_t i) {
      for (int q : P.neighbors[i])
        if (label_axis(P.label[static_cast<size_t>(q)]) == label_axis(P.label[i]) &&
            P.faces[static_cast<size_t>(q)].size() >= P.faces[i].size())
          return true;
      return false;
    };
    int victim = -1;
    size_t victim_size = 0;
    bool victim_conflict = false;
    for (size_t i = 0; i < P.label.size(); ++i) {
      const bool low_degree = P.neighbors[i].size() < 3;
      const bool conflict = !low_degree && conflicting(i);
      if (!low_degree && !conflict) continue;
      if (victim < 0 || P.faces[i].size() < victim_size) {
        victim = static_cast<int>(i);
        victim_size = P.faces[i].size();
        victim_conflict = conflict;
      }
    }
    if (victim < 0) break;
    const auto& nb = P.neighbors[static_cast<size_t>(victim)];
    if (nb.empty()) break;  // a single patch covering a component; nothing to merge into
    int target = nb[0];
    if (nb.size() >= 2) {
      // Neighbor sharing the longest boundary.
      std::map<int, double> len;
      for (const FacePair& p : topo.pairs) {
        int a = P.face_patch(p.f0), b = P.face_patch(p.f1);
        if (a == victim && b != victim) len[b] += (X.row(p.v0) - X.row(p.v1)).norm();
        if (b == victim && a != victim) len[a] += (X.row(p.v0) - X.row(p.v1)).norm();
      }
      double best = -1.0;
      for (int q : nb)
        if (len[q] > best) {
          best = len[q];
          target = q;
        }
    }
    if (victim_conflict)
      ++rep.conflict_merged;
    else if (nb.size() == 2)
      ++rep.degree2_merged;
    else
      ++rep.degree1_merged;
    for (int f : P.faces[static_cast<size_t>(victim)]) labels(f) = P.label[static_cast<size_t>(target)];
    rep.zigzag_faces += remove_zigzags(topo, labels);
  }
  return rep;
}

}  // namespace pcp
