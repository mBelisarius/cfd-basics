#ifndef CFD_BASICS_BUILDER_SYSTEM_HPP_
#define CFD_BASICS_BUILDER_SYSTEM_HPP_

#include "nuenv/core"

namespace cfd_basics {

template<typename Scalar>
struct System {
  System() = default;

  System(nuenv::Index nVolumes, nuenv::MatrixSQX<Scalar> coeffs, nuenv::VectorX<Scalar> constants)
      : NVolumes(nVolumes), Coeffs(coeffs), Constants(constants) {}

  nuenv::Index NVolumes {};
  nuenv::MatrixSQX<Scalar> Coeffs;
  nuenv::VectorX<Scalar> Constants;
};

} // cfd_basics

#endif // CFD_BASICS_BUILDER_SYSTEM_HPP_
