#ifndef CFD_BASICS_MESH_POLYMESH_H_
#define CFD_BASICS_MESH_POLYMESH_H_

#include "nuenv/core"

#include "mesh/cell.hpp"
#include "mesh/face.hpp"

namespace cfd_basics {

template<typename Scalar, typename ScalarField>
class PolyMesh {
 public:
  PolyMesh() = default;

  PolyMesh(nuenv::Index nCells,
           const nuenv::VectorX<ScalarField>& points,
           const nuenv::VectorX<Face>& faces);

  nuenv::Index NCells() const { return nCells_; }

  const nuenv::VectorX<ScalarField>& Points() const { return points_; }

  const nuenv::VectorX<Face>& Faces() const { return faces_; }

  const nuenv::VectorX<Cell>& Cells() const { return cells_; }

  ScalarField FaceNormal(nuenv::Index faceId);

  ScalarField FaceNormal3p(nuenv::Index faceId);

  ScalarField FaceCentre(nuenv::Index faceId);

  nuenv::VectorX<nuenv::Vector2X<Scalar>> FaceProject2D(nuenv::Index faceId);

  Scalar FaceArea(nuenv::Index faceId);

  ScalarField CellCentre(nuenv::Index cellId);

  Scalar CellVolume(nuenv::Index cellId);

 private:
  nuenv::Index nCells_;
  nuenv::VectorX<ScalarField> points_;
  nuenv::VectorX<Face> faces_;
  nuenv::VectorX<Cell> cells_;
};

template<typename Scalar, typename ScalarField>
PolyMesh<Scalar, ScalarField>::PolyMesh(
    nuenv::Index nCells,
    const nuenv::VectorX<ScalarField>& points,
    const nuenv::VectorX<Face>& faces)
    : nCells_(nCells), points_(points), faces_(faces), cells_(nCells) {
  for (const Face& face : faces_) {
    cells_[face.OwnerId()].SetId(face.OwnerId());
    cells_[face.OwnerId()].AddFace(face.Id());

    if (face.NeighbourId() > -1) {
      cells_[face.NeighbourId()].SetId(face.NeighbourId());
      cells_[face.NeighbourId()].AddFace(face.Id());
    }
  }
}

template<typename Scalar, typename ScalarField>
ScalarField PolyMesh<Scalar, ScalarField>::FaceNormal(nuenv::Index faceId) {
  auto facePointsId = faces_[faceId].PointsId();
  nuenv::Index nPoints = facePointsId.size();

  ScalarField normal {0.0, 0.0, 0.0};

  for (nuenv::Index i = 0; i < nPoints; ++i) {
    const ScalarField& current = points_(facePointsId[i]);
    const ScalarField& next = points_(facePointsId[(i + 1) % nPoints]);
    normal += (current - next).cross(current + next);
  }

  normal.normalize();

  return normal;
}

template<typename Scalar, typename ScalarField>
ScalarField PolyMesh<Scalar, ScalarField>::FaceNormal3p(nuenv::Index faceId) {
  auto facePointsId = faces_[faceId].PointsId();

  ScalarField p1 = points_(facePointsId[0]);
  ScalarField p2 = points_(facePointsId[1]);
  ScalarField p3 = points_(facePointsId[2]);

  ScalarField v1 = p2 - p1;
  ScalarField v2 = p3 - p1;

  ScalarField normal = v1.cross(v2);
  normal.normalize();

  return normal;
}

template<typename Scalar, typename ScalarField>
ScalarField PolyMesh<Scalar, ScalarField>::FaceCentre(nuenv::Index faceId) {
  auto face = faces_[faceId];
  auto facePointsId = faces_[faceId].PointsId();
  ScalarField centroid = {0.0, 0.0, 0.0};

  for (const auto& pointId : facePointsId) {
    centroid += points_[pointId];
  }

  nuenv::Index order = face.Order();
  centroid /= static_cast<Scalar>(order);

  return centroid;
}

template<typename Scalar, typename ScalarField>
nuenv::VectorX<nuenv::Vector2X<Scalar>> PolyMesh<Scalar, ScalarField>::FaceProject2D(nuenv::Index faceId) {
  auto facePointsId = faces_[faceId].PointsId();
  nuenv::Index order = faces_[faceId].Order();
  ScalarField normal = FaceNormal(faceId);
  ScalarField centre = FaceCentre(faceId);

  ScalarField e1 = (points_[facePointsId[1]] - points_[facePointsId[0]]);
  e1.normalize();
  ScalarField e2 = e1.cross(-normal);
  e2.normalize();

  nuenv::VectorX<nuenv::Vector2X<Scalar>> projected(order);
  for (nuenv::Index i = 0; i < order; ++i) {
    ScalarField vec = points_[facePointsId[i]] - centre;
    projected[i] = nuenv::Vector2X<Scalar> {vec.dot(e1), vec.dot(e2)};
  }

  return projected;
}

template<typename Scalar, typename ScalarField>
Scalar PolyMesh<Scalar, ScalarField>::FaceArea(nuenv::Index faceId) {
  // Shoelace Theorem
  auto points = FaceProject2D(faceId);
  nuenv::Index order = faces_[faceId].Order();

  Scalar area = 0.0;
  for (nuenv::Index i = 0; i < order; ++i) {
    nuenv::Index j = (i + 1) % order;
    area += points[i][0] * points[j][1] - points[i][1] * points[j][0];
  }

  return abs(area) / 2.0;
}

template<typename Scalar, typename ScalarField>
ScalarField PolyMesh<Scalar, ScalarField>::CellCentre(nuenv::Index cellId) {
  ScalarField centre {0.0, 0.0, 0.0};

  for (auto faceId : cells_[cellId].FacesId()) {
    centre += FaceCentre(faceId);
  }

  centre /= static_cast<Scalar>(cells_[cellId].NFaces());
  return centre;
}

template<typename Scalar, typename ScalarField>
Scalar PolyMesh<Scalar, ScalarField>::CellVolume(nuenv::Index cellId) {
  const Cell& cell = cells_[cellId];
  ScalarField cellCentre = CellCentre(cellId);
  Scalar volume = 0.0;

  for (auto faceId : cell.FacesId()) {
    ScalarField faceCentre = FaceCentre(faceId);
    ScalarField faceNormal = FaceNormal3p(faceId);
    Scalar faceArea = FaceArea(faceId);
    ScalarField heightVec = faceCentre - cellCentre;
    volume += (1.0 / 3.0) * faceArea * abs(heightVec.dot(faceNormal));
  }

  return volume;
}

} // cfd_basics

#endif //CFD_BASICS_MESH_POLYMESH_H_
