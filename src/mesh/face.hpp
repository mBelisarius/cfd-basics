#ifndef CFD_BASICS_MESH_FACE_HPP_
#define CFD_BASICS_MESH_FACE_HPP_

#include "nuenv/core"

namespace cfd_basics {

class Face {
 public:
  Face() = default;

  Face(nuenv::Index id, nuenv::Index owner, nuenv::Index neighbour,
       const nuenv::VectorX<nuenv::Index>& points)
      : order_(points.size()), id_(id), owner_(owner), neighbour_(neighbour),
        pointsId_(points) {}

  nuenv::Index Order() const {
    return order_;
  }

  nuenv::Index Id() const {
    return id_;
  }

  nuenv::Index Owner() const {
    return owner_;
  }

  nuenv::Index Neighbour() const {
    return neighbour_;
  }

  const nuenv::VectorX<nuenv::Index>& Points() const {
    return pointsId_;
  }

 private:
  nuenv::Index order_;
  nuenv::Index id_;
  nuenv::Index owner_;
  nuenv::Index neighbour_;
  nuenv::VectorX<nuenv::Index> pointsId_;
};

} // cfd_basics

#endif // CFD_BASICS_MESH_FACE_HPP_
