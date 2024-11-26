#ifndef CFD_BASICS_BOUNDARY_BOUNDARY_FIELD_HPP_
#define CFD_BASICS_BOUNDARY_BOUNDARY_FIELD_HPP_

#include "nuenv/core"

namespace cfd_basics {

template<typename Scalar, typename ScalarField>
class BoundaryField {
 public:
  BoundaryField() = default;

  BoundaryField(nuenv::Index id, nuenv::Index type, nuenv::Lambda<Scalar(ScalarField)> value)
      : id_(id), type_(type), value_(value) {}

  nuenv::Index Id() { return id_; }

  nuenv::Index Type() { return type_; }

  Scalar Value(const ScalarField& x) { return value_(x); }

 private:
  nuenv::Index id_;
  nuenv::Index type_;
  nuenv::Lambda<Scalar(ScalarField)> value_;
};

} // cfd_basics

#endif //CFD_BASICS_BOUNDARY_BOUNDARY_FIELD_HPP_
