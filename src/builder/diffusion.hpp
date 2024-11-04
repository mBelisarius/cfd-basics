#ifndef CFD_BASICS_BUILDER_DIFFUSION_HPP_
#define CFD_BASICS_BUILDER_DIFFUSION_HPP_


#include "nuenv/core"

#include "mesh/poly_boundary.hpp"
#include "mesh/poly_mesh.hpp"

#include "boundary_field.hpp"
#include "cfd_utils.hpp"
#include "properties.hpp"
#include "system.hpp"

namespace cfd_basics {

template<typename Scalar, typename ScalarField>
System<Scalar> BuilderDiffusion(
    PolyMesh<Scalar, ScalarField> polyMesh,
    nuenv::VectorX<PolyBoundary<Scalar>> polyBoundaries,
    nuenv::VectorX<BoundaryField<Scalar>> boundaryFields,
    PropertiesList<Scalar> properties
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

    // TODO: kFace with Patankar
    Scalar kOwner = properties[Property::kCondutivity];
    Scalar kNeighbour = properties[Property::kCondutivity];
    Scalar kFace = 0.5 * (kOwner + kNeighbour);
    Scalar faceArea = polyMesh.FaceArea(i);
    Scalar cellsDist = (polyMesh.CellCentre(face.OwnerId()) - polyMesh.CellCentre(face.NeighbourId())).norm();
    Scalar coeff_a = kFace * faceArea / cellsDist;

    coeffs(face.OwnerId(), face.OwnerId()) += coeff_a;
    coeffs(face.OwnerId(), face.NeighbourId()) -= coeff_a;
    coeffs(face.NeighbourId(), face.NeighbourId()) += coeff_a;
    coeffs(face.NeighbourId(), face.OwnerId()) -= coeff_a;
  }

  for (const Cell& cell : polyMesh.Cells()) {
    Scalar qdot = properties[Property::kHeatSource];
    constants[cell.Id()] += qdot * polyMesh.CellVolume(cell.Id());
  }

  for (const auto& boundary : polyBoundaries) {
    switch (boundary.Type()) {
      case 0:  // Empty
        break;
      case 1:  // Dirichlet
        for (nuenv::Index faceId = boundary.StartFace(); faceId < boundary.StartFace() + boundary.NFaces(); ++faceId) {
          const auto& face = polyMesh.Faces()[faceId];

          // TODO: kFace with Patankar
          Scalar kOwner = properties[Property::kCondutivity];
          Scalar kNeighbour = properties[Property::kCondutivity];
          Scalar kFace = 0.5 * (kOwner + kNeighbour);
          Scalar faceArea = polyMesh.FaceArea(faceId);
          Scalar boundDist = (polyMesh.CellCentre(face.OwnerId()) - polyMesh.FaceCentre(faceId)).norm();
          Scalar valueBound = boundaryFields[boundary.Id()].Value();

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
