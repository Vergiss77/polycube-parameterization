#include "pcp/core/block_sparse.h"

#include <algorithm>

namespace pcp {

BlockSparse::BlockSparse(Index n_nodes, const std::vector<std::vector<int>>& elements)
    : n_nodes_(n_nodes), nbrs_(static_cast<size_t>(n_nodes)) {
  for (Index i = 0; i < n_nodes; ++i) nbrs_[static_cast<size_t>(i)].push_back(static_cast<int>(i));
  for (const auto& e : elements)
    for (int a : e)
      for (int b : e)
        if (a != b) nbrs_[static_cast<size_t>(b)].push_back(a);
  Index nnz = 0;
  for (auto& n : nbrs_) {
    std::sort(n.begin(), n.end());
    n.erase(std::unique(n.begin(), n.end()), n.end());
    nnz += static_cast<Index>(n.size()) * 9;
  }
  const Index dim = 3 * n_nodes;
  M_.resize(dim, dim);
  M_.resizeNonZeros(nnz);
  int* outer = M_.outerIndexPtr();
  int* inner = M_.innerIndexPtr();
  Index pos = 0;
  for (Index j = 0; j < n_nodes; ++j) {
    const auto& n = nbrs_[static_cast<size_t>(j)];
    for (int q = 0; q < 3; ++q) {
      outer[3 * j + q] = static_cast<int>(pos);
      for (int i : n)
        for (int p = 0; p < 3; ++p) inner[pos++] = 3 * i + p;
    }
  }
  outer[dim] = static_cast<int>(pos);
  set_zero();
}

void BlockSparse::set_zero() {
  std::fill(M_.valuePtr(), M_.valuePtr() + M_.nonZeros(), 0.0);
}

int BlockSparse::local(int i, int j) const {
  const auto& n = nbrs_[static_cast<size_t>(j)];
  auto it = std::lower_bound(n.begin(), n.end(), i);
  if (it == n.end() || *it != i) return -1;
  return static_cast<int>(it - n.begin());
}

void BlockSparse::add_block(int i, int j, const Eigen::Matrix3d& B) {
  int l = local(i, j);
  if (l < 0) throw Error("BlockSparse: block outside the sparsity pattern");
  const int* outer = M_.outerIndexPtr();
  double* val = M_.valuePtr();
  for (int q = 0; q < 3; ++q) {
    int base = outer[3 * j + q] + 3 * l;
    for (int p = 0; p < 3; ++p) val[base + p] += B(p, q);
  }
}

void BlockSparse::add_diag(int i, int j, double s) {
  int l = local(i, j);
  if (l < 0) throw Error("BlockSparse: block outside the sparsity pattern");
  const int* outer = M_.outerIndexPtr();
  double* val = M_.valuePtr();
  for (int q = 0; q < 3; ++q) val[outer[3 * j + q] + 3 * l + q] += s;
}

}  // namespace pcp
