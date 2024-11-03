#ifndef CFD_BASICS_SCRIPTS_HPP_
#define CFD_BASICS_SCRIPTS_HPP_

#include "nuenv/core"

#include "mesh/poly_boundary.hpp"
#include "mesh/poly_mesh.hpp"

namespace cfd_basics {

template<typename Scalar>
nuenv::VectorX<PolyBoundary<Scalar>> MakePolyBoundaries1D(
    nuenv::Index volumes,
    nuenv::Index order
) {
  nuenv::VectorX<PolyBoundary<Scalar>> boundaries(3);

  // Left boudnary
  boundaries[0] = PolyBoundary<Scalar>(0, 1, 1, volumes - 1);

  // Right boundary
  boundaries[1] = PolyBoundary<Scalar>(1, 1, 1, volumes);

  // Adiabatic boundary
  boundaries[2] = PolyBoundary<Scalar>(2, 0, volumes * order, volumes + 1);

  return boundaries;
}

template<typename Scalar, typename ScalarField>
PolyMesh<Scalar, ScalarField> MakePolyMesh1D(
    Scalar length,
    Scalar area,
    nuenv::Index volumes,
    nuenv::Index order
) {
  nuenv::Index nPoints = (volumes + 1) * order;
  nuenv::VectorX<ScalarField> points(nPoints);
  Scalar theta0 = -nuenv::pi / static_cast<Scalar>(order);
  Scalar dTheta = 2.0 * nuenv::pi / static_cast<Scalar>(order);
  Scalar sideLength = sqrt(4.0 * tan(nuenv::pi / static_cast<Scalar>(order)) * area / static_cast<Scalar>(order));
  Scalar radius = sideLength / (2.0 * sin(nuenv::pi / static_cast<Scalar>(order)));

  for (nuenv::Index i = 0; i < volumes + 1; ++i) {
    for (nuenv::Index j = 0; j < order; ++j) {
      points[i * order + j][0] = length * static_cast<Scalar>(i) / static_cast<Scalar>(volumes);
      points[i * order + j][1] = radius * std::cos(theta0 + dTheta * static_cast<Scalar>(j));
      points[i * order + j][2] = radius * std::sin(theta0 + dTheta * static_cast<Scalar>(j));
    }
  }

  nuenv::VectorX<Face> faces((volumes + 1) + (volumes * order));
  nuenv::Index faceCount = 0;

  // Internal faces
  for (nuenv::Index i = 1; i < volumes; ++i) {
    nuenv::VectorX<nuenv::Index> face_points(order);
    for (nuenv::Index j = 0; j < order; ++j) {
      face_points[j] = i * order + j;
    }
    faces[faceCount] = Face(faceCount, i, i - 1, face_points);
    ++faceCount;
  }

  // Boundary faces
  {
    nuenv::VectorX<nuenv::Index> face_points(order);
    for (nuenv::Index j = 0; j < order; j++) {
      face_points[j] = 0 * order + j;
    }
    faces[faceCount] = Face(faceCount, 0, -1, face_points);
    ++faceCount;
  }
  {
    nuenv::VectorX<nuenv::Index> face_points(order);
    for (nuenv::Index j = 0; j < order; j++) {
      face_points[j] = volumes * order + j;
    }
    faces[faceCount] = Face(faceCount, volumes - 1, -1, face_points);
    ++faceCount;
  }

  // Adiabatic faces
  for (nuenv::Index i = 0; i < volumes; ++i) {
    for (nuenv::Index j = 0; j < order; ++j) {
      nuenv::VectorX<nuenv::Index> face_points(4);
      face_points[0] = i * order + j;
      face_points[1] = i * order + (j + 1) % order;
      face_points[2] = (i + 1) * order + (j + 1) % order;
      face_points[3] = (i + 1) * order + j;
      faces[faceCount] = Face(faceCount, i, -1, face_points);
      ++faceCount;
    }
  }

  return PolyMesh<Scalar, ScalarField>(volumes, points, faces);
}

} // cfd_basics

#endif // CFD_BASICS_SCRIPTS_HPP_
