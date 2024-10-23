#ifndef CFD_BASICS_MESH_POLYMESH_H_
#define CFD_BASICS_MESH_POLYMESH_H_

#include "nuenv/core"

#include "mesh/face.hpp"

namespace cfd_basics {

template<typename Scalar, typename ScalarField>
class PolyMesh {
 public:
  PolyMesh() = default;

  PolyMesh(const nuenv::VectorX<ScalarField>& points,
           const nuenv::VectorX<Face>& faces)
      : points_(points), faces_(faces) {}

  ScalarField FaceNormal(nuenv::Index faceId) {
    auto facePointsId = faces_[faceId].Points();

    ScalarField p1 = points_(facePointsId[0]);
    ScalarField p2 = points_(facePointsId[1]);
    ScalarField p3 = points_(facePointsId[2]);

    ScalarField v1 = p2 - p1;
    ScalarField v2 = p3 - p1;

    ScalarField normal = v1.cross(v2);
    normal.normalize();

    return normal;
  }

  ScalarField FaceCentre(nuenv::Index faceId) {
    ScalarField centroid = {0.0, 0.0, 0.0};
    auto face = faces_[faceId];

    for (const ScalarField& point : face.Points()) {
      centroid[0] += point[0];
      centroid[1] += point[1];
      centroid[2] += point[2];
    }

    Scalar order = face.Order();
    centroid[0] /= order;
    centroid[1] /= order;
    centroid[2] /= order;

    return centroid;
  }

 private:
  nuenv::VectorX<ScalarField> points_;
  nuenv::VectorX<Face> faces_;
};

} // cfd_basics

#endif //CFD_BASICS_MESH_POLYMESH_H_
