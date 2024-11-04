#ifndef CFD_BASICS_BOUNDARY_BOUNDARY_FIELD_HPP_
#define CFD_BASICS_BOUNDARY_BOUNDARY_FIELD_HPP_

#include "nuenv/core"

namespace cfd_basics {

template<typename Scalar>
class BoundaryField {
 public:
  BoundaryField() = default;

  BoundaryField(nuenv::Index id, nuenv::Index type, Scalar value)
      : id_(id), type_(type), value_(value) {}

  nuenv::Index Id() { return id_; }

  nuenv::Index Type() { return type_; }

  Scalar Value() { return value_; }

 private:
  nuenv::Index id_;
  nuenv::Index type_;
  Scalar value_;
};

} // cfd_basics

#endif //CFD_BASICS_BOUNDARY_BOUNDARY_FIELD_HPP_
