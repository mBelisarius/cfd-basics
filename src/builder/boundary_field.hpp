#ifndef CFD_BASICS_BOUNDARY_BOUNDARY_FIELD_HPP_
#define CFD_BASICS_BOUNDARY_BOUNDARY_FIELD_HPP_

#include "nuenv/core"

namespace cfd_basics {

template<typename Scalar>
class BoundaryField {
 public:
  BoundaryField(nuenv::Index index, nuenv::Index type, Scalar value)
      : index_(index), type_(type), value_(value) {}

 private:
  nuenv::Index index_;
  nuenv::Index type_;
  Scalar value_;
};

} // cfd_basics

#endif //CFD_BASICS_BOUNDARY_BOUNDARY_FIELD_HPP_
