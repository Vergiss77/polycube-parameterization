#pragma once

#include <filesystem>
#include <array>
#include <map>
#include <set>
#include <vector>
#include <string>

#include "pcp/core/constants.h"
#include "pcp/core/types.h"

namespace testutil {

inline std::filesystem::path models_dir() { return PCP_MODELS_DIR; }

inline std::filesystem::path temp_dir() {
  auto d = std::filesystem::temp_directory_path() / "pcp_tests";
  std::filesystem::create_directories(d);
  return d;
}

// Closed axis-aligned box surface [0,sx]x[0,sy]x[0,sz] with n x n quads (2 tris each) per face,
// outward orientation, welded vertices.
inline pcp::SurfaceMesh make_box(int n, double sx = 1, double sy = 1, double sz = 1) {
  std::map<std::array<int, 3>, int> id;
  std::vector<Eigen::Vector3d> verts;
  std::vector<Eigen::Vector3i> tris;
  auto vid = [&](std::array<int, 3> g) {
    auto it = id.find(g);
    if (it != id.end()) return it->second;
    int i = static_cast<int>(verts.size());
    id[g] = i;
    verts.emplace_back(sx * g[0] / n, sy * g[1] / n, sz * g[2] / n);
    return i;
  };
  for (int axis = 0; axis < 3; ++axis)
    for (int side = 0; side < 2; ++side) {
      int u = (axis + 1) % 3, v = (axis + 2) % 3;
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
          auto g = [&](int a, int b) {
            std::array<int, 3> p{};
            p[axis] = side * n;
            p[u] = a;
            p[v] = b;
            return vid(p);
          };
          int a = g(i, j), b = g(i + 1, j), c = g(i + 1, j + 1), d = g(i, j + 1);
          // (u, v, axis) is right-handed, so a->b->c faces +axis.
          if (side == 1) {
            tris.emplace_back(a, b, c);
            tris.emplace_back(a, c, d);
          } else {
            tris.emplace_back(a, c, b);
            tris.emplace_back(a, d, c);
          }
        }
    }
  pcp::SurfaceMesh m;
  m.V.resize(static_cast<pcp::Index>(verts.size()), 3);
  for (size_t i = 0; i < verts.size(); ++i) m.V.row(static_cast<pcp::Index>(i)) = verts[i];
  m.F.resize(static_cast<pcp::Index>(tris.size()), 3);
  for (size_t i = 0; i < tris.size(); ++i) m.F.row(static_cast<pcp::Index>(i)) = tris[i];
  return m;
}

// Torus with major radius R, minor radius r around the Z axis.
inline pcp::SurfaceMesh make_torus(int nu, int nv, double R = 2, double r = 0.7) {
  pcp::SurfaceMesh m;
  m.V.resize(nu * nv, 3);
  m.F.resize(2 * nu * nv, 3);
  for (int i = 0; i < nu; ++i)
    for (int j = 0; j < nv; ++j) {
      double u = 2 * pcp::kPi * i / nu, v = 2 * pcp::kPi * j / nv;
      m.V.row(i * nv + j) << (R + r * std::cos(v)) * std::cos(u),
          (R + r * std::cos(v)) * std::sin(u), r * std::sin(v);
    }
  int f = 0;
  for (int i = 0; i < nu; ++i)
    for (int j = 0; j < nv; ++j) {
      int a = i * nv + j, b = ((i + 1) % nu) * nv + j, c = ((i + 1) % nu) * nv + (j + 1) % nv,
          d = i * nv + (j + 1) % nv;
      m.F.row(f++) << a, b, c;
      m.F.row(f++) << a, c, d;
    }
  return m;
}

inline void write_text(const std::filesystem::path& p, const std::string& s) {
  FILE* f = std::fopen(p.string().c_str(), "wb");
  std::fwrite(s.data(), 1, s.size(), f);
  std::fclose(f);
}

}  // namespace testutil

namespace testutil {

// Icosphere of radius 1 with `levels` subdivisions.
inline pcp::SurfaceMesh make_icosphere(int levels) {
  const double t = (1.0 + std::sqrt(5.0)) / 2.0;
  std::vector<Eigen::Vector3d> v = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0},
                                    {0, -1, t}, {0, 1, t}, {0, -1, -t}, {0, 1, -t},
                                    {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
  std::vector<Eigen::Vector3i> f = {{0, 11, 5}, {0, 5, 1},  {0, 1, 7},   {0, 7, 10}, {0, 10, 11},
                                    {1, 5, 9},  {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                                    {3, 9, 4},  {3, 4, 2},  {3, 2, 6},   {3, 6, 8},  {3, 8, 9},
                                    {4, 9, 5},  {2, 4, 11}, {6, 2, 10},  {8, 6, 7},  {9, 8, 1}};
  for (auto& p : v) p.normalize();
  for (int l = 0; l < levels; ++l) {
    std::map<std::pair<int, int>, int> mid;
    auto m = [&](int a, int b) {
      auto key = std::make_pair(std::min(a, b), std::max(a, b));
      auto it = mid.find(key);
      if (it != mid.end()) return it->second;
      v.push_back((v[a] + v[b]).normalized());
      return mid[key] = static_cast<int>(v.size()) - 1;
    };
    std::vector<Eigen::Vector3i> nf;
    for (auto& tri : f) {
      int a = m(tri[0], tri[1]), b = m(tri[1], tri[2]), c = m(tri[2], tri[0]);
      nf.emplace_back(tri[0], a, c);
      nf.emplace_back(tri[1], b, a);
      nf.emplace_back(tri[2], c, b);
      nf.emplace_back(a, b, c);
    }
    f = nf;
  }
  pcp::SurfaceMesh s;
  s.V.resize(static_cast<pcp::Index>(v.size()), 3);
  for (size_t i = 0; i < v.size(); ++i) s.V.row(static_cast<pcp::Index>(i)) = v[i];
  s.F.resize(static_cast<pcp::Index>(f.size()), 3);
  for (size_t i = 0; i < f.size(); ++i) s.F.row(static_cast<pcp::Index>(i)) = f[i];
  return s;
}

// Closed surface of a union of unit voxels (integer coordinates), each voxel face split into
// n x n quads (2 triangles each), outward orientation, welded vertices. The voxels must not
// touch along edges or vertices only (the result would not be manifold).
inline pcp::SurfaceMesh make_voxels(const std::vector<std::array<int, 3>>& voxels, int n) {
  std::set<std::array<int, 3>> occ(voxels.begin(), voxels.end());
  std::map<std::array<int, 3>, int> id;
  std::vector<Eigen::Vector3d> verts;
  std::vector<Eigen::Vector3i> tris;
  auto vid = [&](std::array<int, 3> g) {
    auto it = id.find(g);
    if (it != id.end()) return it->second;
    int i = static_cast<int>(verts.size());
    id[g] = i;
    verts.emplace_back(double(g[0]) / n, double(g[1]) / n, double(g[2]) / n);
    return i;
  };
  for (const auto& v : voxels) {
    for (int a = 0; a < 3; ++a) {
      for (int sgn : {-1, 1}) {
        std::array<int, 3> nb = v;
        nb[a] += sgn;
        if (occ.count(nb)) continue;
        const int b = (a + 1) % 3, c = (a + 2) % 3;
        for (int i = 0; i < n; ++i)
          for (int j = 0; j < n; ++j) {
            std::array<int, 4> q;
            const int du[4] = {0, 1, 1, 0}, dv[4] = {0, 0, 1, 1};
            for (int k = 0; k < 4; ++k) {
              std::array<int, 3> g{};
              g[a] = (v[a] + (sgn > 0 ? 1 : 0)) * n;
              g[b] = v[b] * n + i + du[k];
              g[c] = v[c] * n + j + dv[k];
              q[k] = vid(g);
            }
            if (sgn > 0) {  // e_b x e_c = e_a: counter-clockwise seen from outside
              tris.emplace_back(q[0], q[1], q[2]);
              tris.emplace_back(q[0], q[2], q[3]);
            } else {
              tris.emplace_back(q[0], q[2], q[1]);
              tris.emplace_back(q[0], q[3], q[2]);
            }
          }
      }
    }
  }
  pcp::SurfaceMesh s;
  s.V.resize(static_cast<pcp::Index>(verts.size()), 3);
  for (size_t i = 0; i < verts.size(); ++i) s.V.row(static_cast<pcp::Index>(i)) = verts[i];
  s.F.resize(static_cast<pcp::Index>(tris.size()), 3);
  for (size_t i = 0; i < tris.size(); ++i) s.F.row(static_cast<pcp::Index>(i)) = tris[i];
  return s;
}

}  // namespace testutil
