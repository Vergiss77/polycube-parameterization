#pragma once

// Polycube height simplification (Chen et al. 2019, section 5): snap the heights of same-axis
// patch groups that are not ordered by any axis-parallel chain, so that fewer distinct heights
// (cuts) remain, without changing the polycube topology.

#include <array>
#include <vector>

#include <nlohmann/json.hpp>

#include "pcp/structure/final_map.h"
#include "pcp/structure/patch_graph.h"

namespace pcp {

struct SimplifyParams {
  bool enabled = true;
  // Every axis-parallel chain keeps its length ratio r / r_rest within [1/s, s].
  double max_stretch = 2.0;
  nlohmann::json to_json() const;
  void update_from_json(const nlohmann::json& j);  // throws pcp::Error on unknown keys
};

// A set of patches on one axis sharing one height.
struct HeightGroup {
  int axis = -1;
  double height = 0.0;
  std::vector<int> patches;
};

// Directed edge from a lower group to an upper group of the same axis: either one
// axis-parallel chain (Chen 2019) or a separation between two patches facing each other along
// the axis across material or across a gap (extension, see design D1b).
struct HeightEdge {
  int from = -1, to = -1;
  double rest = 0.0;  // height difference before simplification (> 0)
  bool separation = false;
  bool width = false;
};

// Mesh-independent input of the simplification.
struct HeightGraph {
  int n_patches = 0;
  std::vector<HeightGroup> groups;
  std::vector<HeightEdge> edges;
  int chains_without_ends = 0;  // chains lacking a lower or upper patch (no constraint)
  int separation_edges = 0;
  int width_edges = 0;
  nlohmann::json to_json() const;
  static HeightGraph from_json(const nlohmann::json& j);
};

// Groups: patches with the same axis and the same height (up to round-off, 1e-12 relative).
// Edges: one per chain with both lower and upper patch, rest length = height difference.
// Throws pcp::Error if a chain's upper patch is not strictly higher than its lower patch.
HeightGraph build_height_graph(const PatchGraph& g, const std::vector<double>& heights);

// Adds separation edges: from the centroid of every boundary triangle of the axis-aligned
// polycube X, rays along -/+ its axis; the first hit triangle of the same axis and opposite
// orientation gives an edge between the two patches (lower -> upper), one per patch pair.
// Returns the number of edges added.
int add_separation_edges(HeightGraph& hg, const PatchGraph& g, const std::vector<double>& heights,
                         const MatX3d& X, const MatX3i& B);

// Adds width edges: inside every planar patch P of the axis-aligned polycube X, from the
// centroid of each of its triangles, rays within P's plane along -/+ each of the two in-plane
// axes a; the first boundary edges hit on both sides border patches Q1 (below) and Q2 (above)
// of axis a, and the strip of P between them must keep a positive width: edge Q1 -> Q2.
// One edge per patch pair. Returns the number of edges added.
int add_width_edges(HeightGraph& hg, const PatchGraph& g, const std::vector<double>& heights,
                    const MatX3d& X, const MatX3i& B);

struct MergeRecord {
  int a = -1, b = -1;  // merged groups (initial ids; the merged group keeps id a < b)
  double height = 0.0;
  double cost = 0.0;   // increase of the total spring distortion
};

struct PairCost {
  bool feasible = false;
  double height = 0.0;
  double cost = 0.0;
};

// Current state of a simplification: representative and height of every initial group.
struct HeightState {
  std::vector<int> rep;         // representative group id (smallest id of the cluster)
  std::vector<double> height;   // valid for representatives
  static HeightState initial(const HeightGraph& g);
};

// Cost of merging representatives a and b (eq. 8 with the cost taken as the increase of the
// total distortion sum_l d(l), d = r_rest / r + r / r_rest), subject to the stretch bound.
PairCost merge_cost(const HeightGraph& g, const HeightState& st, int a, int b, double max_stretch);

struct SimplifyResult {
  std::vector<int> group_of_patch;   // final group (representative id) of every patch
  std::vector<double> patch_height;  // final height of every patch
  std::vector<int> group_rep;        // representative of every initial group
  std::vector<double> group_height;  // final height of every initial group
  std::vector<MergeRecord> merges;
  std::array<int, 3> cuts_before{}, cuts_after{};
  double min_ratio = 1.0, max_ratio = 1.0;  // over all edges, r / r_rest after simplification
  int chain_edges = 0, separation_edges = 0, width_edges = 0;
  int merges_proposed = 0;  // merges found by the greedy search (merges.size() are kept)
  int cuts_total_before() const { return cuts_before[0] + cuts_before[1] + cuts_before[2]; }
  int cuts_total_after() const { return cuts_after[0] + cuts_after[1] + cuts_after[2]; }
  nlohmann::json to_json() const;
};

// Greedy directed-graph simplification (Algorithm 1): repeatedly merges the feasible pair of
// same-axis, mutually unreachable groups with the smallest cost, until none is feasible.
SimplifyResult simplify_heights(const HeightGraph& g, double max_stretch);

// Result of keeping only the first k of the given merges (in order).
SimplifyResult apply_merges(const HeightGraph& g, const std::vector<MergeRecord>& merges, size_t k);

// Builds the height graph (chain, separation and width edges) of the exact polycube X_prev,
// then runs remap_with_rollback.
struct SimplifiedPolycube {
  HeightGraph graph;  // input of the simplification (chain and separation edges)
  SimplifyResult simplify;
  FinalMapResult map;
  TopologyCheck topology;
  double seconds_simplify = 0.0;
  // Rollback (design D1d): inversions of the fixed-height map (without untangling) for the
  // heights before simplification and for the kept merges, and the number of such maps solved.
  Index inverted_reference = 0;
  Index inverted_kept = 0;
  int probe_maps = 0;
};

// Greedy simplification of hg, then keeps the longest prefix of the merges whose re-mapping has
// no more inverted tets than the re-mapping with the original heights, and re-maps with it.
SimplifiedPolycube remap_with_rollback(const MatX3d& rest, const MatX3d& X_prev, const MatX4i& T,
                                       const MatX3i& B, const PatchGraph& g, const HeightGraph& hg,
                                       const SimplifyParams& sp, FinalMapParams fp);

SimplifiedPolycube simplify_and_remap(const MatX3d& rest, const MatX3d& X_prev, const MatX4i& T,
                                      const MatX3i& B, const PatchGraph& g,
                                      const std::vector<double>& heights,
                                      const SimplifyParams& sp, FinalMapParams fp);

}  // namespace pcp
