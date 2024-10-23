#ifndef CFD_BASICS_NUMCAL_INTEGRATION_H_
#define CFD_BASICS_NUMCAL_INTEGRATION_H_

#include "nuenv/core"

namespace cfd_basics {

template<typename Scalar>
Scalar Riemann(Scalar (* func)(Scalar), Scalar x0, Scalar x1, nuenv::Index n) {
  Scalar xp;
  Scalar sol = 0;

  // Memoization
  Scalar n_inv = (Scalar)1 / n;
  Scalar dx = (x1 - x0) / n;

  for (nuenv::Index p = 0; p < n; p++) {
    xp = (p + 0.5) * dx;
    sol += func(xp) * n_inv;
  }

  return sol;
}

} // cfd_basics

#endif //CFD_BASICS_NUMCAL_INTEGRATION_H_
