#include <doctest/doctest.h>
#include <Eigen/Sparse>
#include <nlohmann/json.hpp>
#include "pcp/core/version.h"

TEST_CASE("smoke: dependencies link") {
  Eigen::SparseMatrix<double> A(2, 2);
  A.insert(0, 0) = 2.0;
  A.insert(1, 1) = 4.0;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(A);
  Eigen::VectorXd x = ldlt.solve(Eigen::VectorXd::Ones(2));
  CHECK(x(1) == doctest::Approx(0.25));
  nlohmann::json j = {{"v", pcp::version()}};
  CHECK(j["v"] == "0.1.0");
}
