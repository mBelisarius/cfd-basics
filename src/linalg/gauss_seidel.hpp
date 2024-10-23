#ifndef CFD_BASICS_linalg_GAUSS_SEIDEL_H_
#define CFD_BASICS_linalg_GAUSS_SEIDEL_H_

#include "nuenv/core"

#include "solution.hpp"

namespace cfd_basics {
template<typename Scalar>
Solution<Scalar> GaussSeidel(const nuenv::MatrixSQX<Scalar>& a,
                             const nuenv::VectorX<Scalar>& b,
                             Scalar rtol = 1e-5,
                             nuenv::Index maxIter = 100) {
  nuenv::VectorX<Scalar> x0 = nuenv::VectorX<Scalar>::Zero(b.size());
  return GaussSeidel(a, b, x0, rtol, maxIter);
}

template<typename Scalar>
Solution<Scalar> gaussSeidel(const nuenv::MatrixSQX<Scalar>& a,
                             const nuenv::VectorX<Scalar>& b,
                             const nuenv::VectorX<Scalar>& x0,
                             Scalar rtol = 1e-5,
                             nuenv::size_t maxIter = 100) {
  bool status = false;
  string message;
  nuenv::Index iter;

  Scalar sum1;
  Scalar sum2;
  nuenv::VectorX<Scalar> x = x0;
  nuenv::VectorX<Scalar> x1;

  for (iter = 1; iter <= maxIter; iter++) {
    x1 = x;

    for (nuenv::Index i = 0; i < x0.size(); i++) {
      sum1 = 0;
      sum2 = 0;

      for (nuenv::Index j = 0; j < i; j++) {
        sum1 += a(i, j) * x(j);
      }

      for (nuenv::Index j = i + 1; j < x0.size(); j++) {
        sum2 += a(i, j) * x1(j);
      }

      x(i) = (-sum1 - sum2 + b(i)) / a(i, i);
    }

    if ((x - x1).norm() < x.norm() * rtol) {
      status = true;
      message = "Convergence achieved";
      break;
    }
  }

  if (!status) { message = "Convergence not achieved in Gauss-Seidel method"; }

  return Solution<Scalar>(status, message, iter, x);
}

} // cfd_basics

#endif //CFD_BASICS_linalg_GAUSS_SEIDEL_H_
