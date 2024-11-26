#ifndef CFD_BASICS_SCRIPTS_HPP_
#define CFD_BASICS_SCRIPTS_HPP_

#include "nuenv/core"

#include "mesh/poly_boundary.hpp"
#include "mesh/poly_mesh.hpp"

namespace cfd_basics {

template<typename Scalar>
nuenv::VectorX<PolyBoundary<Scalar>> MakePolyBoundaries1d(
    nuenv::Index volumes,
    nuenv::Index order
) {
  nuenv::VectorX<PolyBoundary<Scalar>> boundaries(3);

  // Left boundary
  boundaries[0] = PolyBoundary<Scalar>(0, 1, 1, volumes - 1);

  // Right boundary
  boundaries[1] = PolyBoundary<Scalar>(1, 1, 1, volumes);

  // Adiabatic boundary
  boundaries[2] = PolyBoundary<Scalar>(2, 0, volumes * order, volumes + 1);

  return boundaries;
}

template<typename Scalar, typename ScalarField>
PolyMesh<Scalar, ScalarField> MakePolyMesh1d(
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

template<typename Scalar>
nuenv::VectorX<PolyBoundary<Scalar>> MakePolyBoundaries2dRegular(
    nuenv::Index volumesX,
    nuenv::Index volumesY
) {
//  nuenv::Index order = 4;
  nuenv::VectorX<PolyBoundary<Scalar>> boundaries(5);

  nuenv::Index internalFaces = (volumesX - 1) * volumesY + volumesX * (volumesY - 1);

  // Down boundary
  boundaries[0] = PolyBoundary<Scalar>(0, 1, volumesX, internalFaces);

  // Right boundary
  boundaries[1] = PolyBoundary<Scalar>(1, 1, volumesY, internalFaces + volumesX);

  // Up boundary
  boundaries[2] = PolyBoundary<Scalar>(2, 1, volumesX, internalFaces + volumesX + volumesY);

  // Left boundary
  boundaries[3] = PolyBoundary<Scalar>(3, 1, volumesY, internalFaces + 2 * volumesX + volumesY);

  // Adiabatic boundary
  boundaries[4] = PolyBoundary<Scalar>(4, 0, 2 * volumesX * volumesY, internalFaces + 2 * volumesX + 2 * volumesY);

  return boundaries;
}

template<typename Scalar, typename ScalarField>
PolyMesh<Scalar, ScalarField> MakePolyMesh2dRegular(
    Scalar lengthX,
    Scalar lengthY,
    nuenv::Index volumesX,
    nuenv::Index volumesY
) {
  nuenv::Index order = 4;

  nuenv::Index nPoints = 2 * (volumesX + 1) * (volumesY + 1);
  nuenv::VectorX<ScalarField> points(nPoints);

  for (nuenv::Index j = 0; j < volumesY + 1; ++j) {
    for (nuenv::Index i = 0; i < volumesX + 1; ++i) {
      points[i + (volumesX + 1) * j][0] = lengthX * static_cast<Scalar>(i) / static_cast<Scalar>(volumesX);
      points[i + (volumesX + 1) * j][1] = lengthY * static_cast<Scalar>(j) / static_cast<Scalar>(volumesY);
      points[i + (volumesX + 1) * j][2] = -0.5;

      points[i + (volumesX + 1) * j + nPoints / 2][0] =
          lengthX * static_cast<Scalar>(i) / static_cast<Scalar>(volumesX);
      points[i + (volumesX + 1) * j + nPoints / 2][1] =
          lengthY * static_cast<Scalar>(j) / static_cast<Scalar>(volumesY);
      points[i + (volumesX + 1) * j + nPoints / 2][2] = 0.5;
    }
  }

  nuenv::VectorX<Face> faces((volumesX + 1) * volumesY + volumesX * (volumesY + 1) + 2 * volumesX * volumesY);
  nuenv::Index faceCount = 0;

  // Internal faces - vertical
  for (nuenv::Index i = 0; i < (volumesX - 1) * volumesY; ++i) {
    nuenv::VectorX<nuenv::Index> face_points(order);
    nuenv::Index iPoint = i + 2 * (i / (volumesX - 1)) + 1;
    face_points[0] = iPoint;
    face_points[1] = iPoint + nPoints / 2;
    face_points[2] = iPoint + nPoints / 2 + (volumesX + 1);
    face_points[3] = iPoint + (volumesX + 1);

    nuenv::Index iFace = i + i / (volumesX - 1) + 1;
    faces[faceCount] = Face(faceCount, iFace, iFace - 1, face_points);
    ++faceCount;
  }

  // Internal faces - horizontal
  for (nuenv::Index i = 0; i < volumesX * (volumesY - 1); ++i) {
    nuenv::VectorX<nuenv::Index> face_points(order);
    nuenv::Index iPoint = i + (i / volumesX) + volumesX + 1;
    face_points[0] = iPoint + 1;
    face_points[1] = iPoint + nPoints / 2 + 1;
    face_points[2] = iPoint + nPoints / 2;
    face_points[3] = iPoint;

    nuenv::Index iFace = i + volumesX;
    faces[faceCount] = Face(faceCount, iFace, iFace - volumesX, face_points);
    ++faceCount;
  }

  // Boundary faces - Down
  for (nuenv::Index i = 0; i < volumesX; ++i) {
    nuenv::VectorX<nuenv::Index> face_points(order);
    nuenv::Index iPoint = i;
    face_points[0] = iPoint + 1;
    face_points[1] = iPoint + nPoints / 2 + 1;
    face_points[2] = iPoint + nPoints / 2;
    face_points[3] = iPoint;

    nuenv::Index iFace = i;
    faces[faceCount] = Face(faceCount, iFace, -1, face_points);
    ++faceCount;
  }
  // Boundary faces - Right
  for (nuenv::Index i = 0; i < volumesY; ++i) {
    nuenv::VectorX<nuenv::Index> face_points(order);
    nuenv::Index iPoint = (i + 1) * (volumesX + 1) - 1;
    face_points[0] = iPoint;
    face_points[1] = iPoint + nPoints / 2;
    face_points[2] = iPoint + nPoints / 2 + (volumesX + 1);
    face_points[3] = iPoint + (volumesX + 1);

    nuenv::Index iFace = (i + 1) * volumesX - 1;
    faces[faceCount] = Face(faceCount, iFace, -1, face_points);
    ++faceCount;
  }
  // Boundary faces - Up
  for (nuenv::Index i = 0; i < volumesX; ++i) {
    nuenv::VectorX<nuenv::Index> face_points(order);
    nuenv::Index iPoint = i + (volumesX + 1) * volumesY;
    face_points[0] = iPoint + 1;
    face_points[1] = iPoint + nPoints / 2 + 1;
    face_points[2] = iPoint + nPoints / 2;
    face_points[3] = iPoint;

    nuenv::Index iFace = i + volumesX * (volumesY - 1);
    faces[faceCount] = Face(faceCount, iFace, -1, face_points);
    ++faceCount;
  }
  // Boundary faces - Left
  for (nuenv::Index i = 0; i < volumesY; ++i) {
    nuenv::VectorX<nuenv::Index> face_points(order);
    nuenv::Index iPoint = i * (volumesX + 1);
    face_points[0] = iPoint;
    face_points[1] = iPoint + nPoints / 2;
    face_points[2] = iPoint + nPoints / 2 + (volumesX + 1);
    face_points[3] = iPoint + (volumesX + 1);

    nuenv::Index iFace = i * volumesX;
    faces[faceCount] = Face(faceCount, iFace, -1, face_points);
    ++faceCount;
  }

  // Adiabatic faces
  for (nuenv::Index i = 0; i < volumesX * volumesY; ++i) {
    nuenv::VectorX<nuenv::Index> face_points(4);
    face_points[0] = i + i / volumesX + 1;
    face_points[1] = i + i / volumesX + 1 + (volumesX + 1);
    face_points[2] = i + i / volumesX + (volumesX + 1);
    face_points[3] = i + i / volumesX;
    faces[faceCount] = Face(faceCount, i, -1, face_points);
    ++faceCount;
  }
  for (nuenv::Index i = 0; i < volumesX * volumesY; ++i) {
    nuenv::VectorX<nuenv::Index> face_points(4);
    face_points[0] = i + i / volumesX + 1 + nPoints / 2;
    face_points[1] = i + i / volumesX + 1 + (volumesX + 1) + nPoints / 2;
    face_points[2] = i + i / volumesX + (volumesX + 1) + nPoints / 2;
    face_points[3] = i + i / volumesX + nPoints / 2;
    faces[faceCount] = Face(faceCount, i, -1, face_points);
    ++faceCount;
  }

  return PolyMesh<Scalar, ScalarField>(volumesX * volumesY, points, faces);
}

} // cfd_basics

#endif // CFD_BASICS_SCRIPTS_HPP_
