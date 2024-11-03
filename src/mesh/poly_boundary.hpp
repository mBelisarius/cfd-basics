#ifndef CFD_BASICS_MESH_POLY_BOUNDARY_HPP_
#define CFD_BASICS_MESH_POLY_BOUNDARY_HPP_

#include "nuenv/core"

namespace cfd_basics {

template<typename Scalar>
class PolyBoundary {
 public:
  PolyBoundary() = default;

  PolyBoundary(nuenv::Index index, nuenv::Index type, nuenv::Index nFaces, nuenv::Index startFace)
      : index_(index), type_(type), nFaces_(nFaces), startFace_(startFace) {}

  nuenv::Index Index() const { return index_; }

  nuenv::Index Type() const { return type_; }

  nuenv::Index NFaces() const { return nFaces_; }

  nuenv::Index StartFace() const { return startFace_; }

 private:
  nuenv::Index index_;
  nuenv::Index type_;
  nuenv::Index nFaces_;
  nuenv::Index startFace_;
};

} // cfd_basics

#endif // CFD_BASICS_MESH_POLY_BOUNDARY_HPP_
