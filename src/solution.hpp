#ifndef CFD_BASICS_SOLUTION_H_
#define CFD_BASICS_SOLUTION_H_

#include "nuenv/core"

namespace cfd_basics {
template<typename Scalar>
struct Solution {
  Solution() : status(false), message("Not solved"), iterations(0), x {0} {}

  Solution(bool status, std::string message,
           nuenv::Index iterations, nuenv::VectorX<Scalar> x)
      : status(status), message(message),
        iterations(iterations), x(x) {}

  // Copy constructor
  Solution(const Solution<Scalar>& other)
      : status(other.status), message(other.message),
        iterations(other.iterations), x(other.x) {}

  bool status;
  std::string message;
  nuenv::Index iterations;
  nuenv::VectorX<Scalar> x;
};

} // cfd_basics

#endif /* CFD_BASICS_SOLUTION_H_ */
