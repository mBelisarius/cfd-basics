#ifndef CFD_BASICS_MESH_POLY_BOUNDARY_HPP_
#define CFD_BASICS_MESH_POLY_BOUNDARY_HPP_

#include "nuenv/core"

namespace cfd_basics {

template<typename Scalar>
class PolyBoundary {
 public:
  PolyBoundary() = default;

  PolyBoundary(nuenv::Index id, nuenv::Index type, nuenv::Index nFaces, nuenv::Index startFace)
      : id_(id), type_(type), nFaces_(nFaces), startFace_(startFace) {}

  nuenv::Index Id() const { return id_; }

  nuenv::Index Type() const { return type_; }

  nuenv::Index NFaces() const { return nFaces_; }

  nuenv::Index StartFace() const { return startFace_; }

 private:
  nuenv::Index id_;
  nuenv::Index type_;
  nuenv::Index nFaces_;
  nuenv::Index startFace_;
};

} // cfd_basics

#endif // CFD_BASICS_MESH_POLY_BOUNDARY_HPP_
