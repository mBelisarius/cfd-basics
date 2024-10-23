#ifndef CFD_BASICS_GENERATOR_BOUNDARY_FIELD_H_
#define CFD_BASICS_GENERATOR_BOUNDARY_FIELD_H_

#include "nuenv/core"

#include "mesh/poly_mesh.hpp"

namespace cfd_basics {
template<typename Scalar, typename ScalarField>
struct BoundaryField {
  BoundaryField(nuenv::Index index, std::function<Scalar(ScalarField)> condition)
      : index(index), condition(condition) {}

  nuenv::Index index;
  std::function<Scalar(ScalarField)> condition;
};

} // cfd_basics

#endif //CFD_BASICS_GENERATOR_BOUNDARY_FIELD_H_
