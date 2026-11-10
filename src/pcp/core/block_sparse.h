#pragma once

#include <array>
#include <vector>

#include <Eigen/Sparse>

#include "pcp/core/types.h"

namespace pcp {

// Symmetric 3N x 3N sparse matrix made of dense 3x3 node blocks, with a sparsity pattern
// fixed at construction (all node pairs that share an element). Values can be reset and
// accumulated block by block without reallocating, so a sparse factorization can reuse its
// symbolic analysis. Both triangles of the matrix are stored.
class BlockSparse {
 public:
  BlockSparse() = default;
  // elements: node lists; every pair of nodes in an element becomes a nonzero block.
  BlockSparse(Index n_nodes, const std::vector<std::vector<int>>& elements);

  Index n_nodes() const { return n_nodes_; }
  void set_zero();
  // Adds B to block (i, j). Both nodes must share an element of the pattern.
  void add_block(int i, int j, const Eigen::Matrix3d& B);
  // Adds s * I3 to block (i, j).
  void add_diag(int i, int j, double s);
  // Adds an element matrix (3k x 3k) for nodes `nodes` (size k).
  template <int K>
  void add_element(const std::array<int, K>& nodes,
                   const Eigen::Matrix<double, 3 * K, 3 * K>& H) {
    for (int a = 0; a < K; ++a)
      for (int b = 0; b < K; ++b)
        add_block(nodes[a], nodes[b], H.template block<3, 3>(3 * a, 3 * b));
  }

  const Eigen::SparseMatrix<double>& matrix() const { return M_; }
  Eigen::SparseMatrix<double>& matrix() { return M_; }

 private:
  // Position of node i among the sorted neighbors of node j; -1 if absent.
  int local(int i, int j) const;

  Index n_nodes_ = 0;
  std::vector<std::vector<int>> nbrs_;  // sorted, includes self
  Eigen::SparseMatrix<double> M_;
};

}  // namespace pcp
