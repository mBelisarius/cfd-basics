#ifndef CFD_BASICS_MESH_FACE_HPP_
#define CFD_BASICS_MESH_FACE_HPP_

#include "nuenv/core"

namespace cfd_basics {

class Face {
 public:
  Face() = default;

  Face(nuenv::Index id, nuenv::Index ownerId, nuenv::Index neighbourId,
       const nuenv::VectorX<nuenv::Index>& pointsId)
      : order_(pointsId.size()), id_(id), ownerId_(ownerId), neighbourId_(neighbourId),
        pointsId_(pointsId) {}

  nuenv::Index Order() const {
    return order_;
  }

  nuenv::Index Id() const {
    return id_;
  }

  nuenv::Index OwnerId() const {
    return ownerId_;
  }

  nuenv::Index NeighbourId() const {
    return neighbourId_;
  }

  const nuenv::VectorX<nuenv::Index>& PointsId() const {
    return pointsId_;
  }

 private:
  nuenv::Index order_;
  nuenv::Index id_;
  nuenv::Index ownerId_;
  nuenv::Index neighbourId_;
  nuenv::VectorX<nuenv::Index> pointsId_;
};

} // cfd_basics

#endif // CFD_BASICS_MESH_FACE_HPP_
