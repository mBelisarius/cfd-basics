#ifndef CFD_BASICS_BUILDER_DIFFUSION_HPP_
#define CFD_BASICS_BUILDER_DIFFUSION_HPP_


#include "nuenv/core"

#include "mesh/poly_boundary.hpp"
#include "mesh/poly_mesh.hpp"

#include "boundary_field.hpp"
#include "cfd_utils.hpp"
#include "system.hpp"

namespace cfd_basics {

template<typename Scalar, typename ScalarField>
System<Scalar> BuilderDiffusion(
    PolyMesh<Scalar, ScalarField> polyMesh,
    nuenv::VectorX<PolyBoundary<Scalar>> polyBoundaries
//    nuenv::VectorX<BoundaryField<Scalar>> boundaryFields
) {
  nuenv::Index nVolumes = polyMesh.NCells();
  nuenv::MatrixSQX<Scalar> coeffs = nuenv::MatrixSQX<Scalar>::Zero(nVolumes, nVolumes);
  nuenv::VectorX<Scalar> constants = nuenv::VectorX<Scalar>::Zero(nVolumes);

  nuenv::Index boundariesStartFace = polyMesh.Faces().size();
  for (const auto& boundary : polyBoundaries) {
    if (boundary.StartFace() < boundariesStartFace) {
      boundariesStartFace = boundary.StartFace();
    }
  }

  for (nuenv::Index i = 0; i < boundariesStartFace; ++i) {
    const auto& face = polyMesh.Faces()[i];

    // TODO: kFace
    Scalar kFace = 1.0;
    Scalar faceArea = polyMesh.FaceArea(i);
    Scalar cellsDist = (polyMesh.CellCentre(face.OwnerId()) - polyMesh.CellCentre(face.NeighbourId())).norm();
    Scalar coeff_a = kFace * faceArea / cellsDist;

    coeffs(face.OwnerId(), face.OwnerId()) += coeff_a;
    coeffs(face.OwnerId(), face.NeighbourId()) += coeff_a;
    coeffs(face.NeighbourId(), face.NeighbourId()) += coeff_a;
    coeffs(face.NeighbourId(), face.OwnerId()) += coeff_a;

    // Todo: qdot
    Scalar qdot = 0.0;
    constants[face.OwnerId()] += qdot * polyMesh.CellVolume(face.OwnerId());
  }

  for (const auto& boundary : polyBoundaries) {
    switch (boundary.Type()) {
      case 0:  // Empty
        break;
      case 1:  // Dirichlet
        for (nuenv::Index faceId = boundary.StartFace(); faceId < boundary.StartFace() + boundary.NFaces(); ++faceId) {
          const auto& face = polyMesh.Faces()[faceId];

          // TODO: kFace, valueBound
          Scalar kFace = 1.0;
          Scalar faceArea = polyMesh.FaceArea(faceId);
          Scalar boundDist = (polyMesh.CellCentre(face.OwnerId()) - polyMesh.FaceCentre(faceId)).norm();
          Scalar valueBound = 100.0;

          coeffs(face.OwnerId(), face.OwnerId()) += kFace * faceArea / boundDist;
          constants[face.OwnerId()] += valueBound * kFace * faceArea / boundDist;
        }
        break;
      case 2:  // Neumann
        throw;
        break;
      case 3:  // Robin
        throw;
        break;
      default:throw;
    }
  }

  return System<Scalar>(nVolumes, coeffs, constants);
}

} // cfd_basics

#endif /* CFD_BASICS_BUILDER_DIFFUSION_HPP_ */
