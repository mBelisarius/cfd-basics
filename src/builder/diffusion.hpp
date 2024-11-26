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
    nuenv::VectorX<BoundaryField<Scalar, ScalarField>> boundaryFields,
    PropertiesList<Scalar> properties,
    nuenv::VectorX<Scalar> x0,
    Scalar dt,
    Scalar theta = 0.5
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

    // TODO: alphaFace with Patankar
    Scalar alphaOwner = properties[Property::kThermalDiffusivity];
    Scalar alphaNeighbour = properties[Property::kThermalDiffusivity];
    Scalar alphaFace = 0.5 * (alphaOwner + alphaNeighbour);
    Scalar faceArea = polyMesh.FaceArea(i);
    Scalar cellsDist = (polyMesh.CellCentre(face.OwnerId()) - polyMesh.CellCentre(face.NeighbourId())).norm();
    Scalar coeff_a = (alphaFace * faceArea / cellsDist) * theta * dt;

    coeffs(face.OwnerId(), face.OwnerId()) += coeff_a;
    coeffs(face.OwnerId(), face.NeighbourId()) -= coeff_a;
    coeffs(face.NeighbourId(), face.NeighbourId()) += coeff_a;
    coeffs(face.NeighbourId(), face.OwnerId()) -= coeff_a;

    Scalar const_b =
        (alphaFace * faceArea / cellsDist) * (1.0 - theta) * dt * (x0[face.OwnerId()] - x0[face.NeighbourId()]);
    constants[face.OwnerId()] -= const_b;
    constants[face.NeighbourId()] += const_b;
  }

  for (const Cell& cell : polyMesh.Cells()) {
    Scalar cellVolume = polyMesh.CellVolume(cell.Id());
    coeffs(cell.Id(), cell.Id()) += cellVolume;

    constants[cell.Id()] +=
        properties[Property::kHeatSource] * cellVolume * dt * properties[Property::kThermalDiffusivity]
            / properties[Property::kCondutivity];
    constants[cell.Id()] += cellVolume * x0[cell.Id()];
  }

  for (const auto& boundary : polyBoundaries) {
    switch (boundary.Type()) {
      case 0:  // Empty
        break;
      case 1:  // Dirichlet
        for (nuenv::Index faceId = boundary.StartFace(); faceId < boundary.StartFace() + boundary.NFaces(); ++faceId) {
          const auto& face = polyMesh.Faces()[faceId];

          // TODO: alphaFace with Patankar
          Scalar alphaOwner = properties[Property::kThermalDiffusivity];
          Scalar alphaFace = alphaOwner;
          Scalar faceArea = polyMesh.FaceArea(faceId);
          ScalarField faceCentre = polyMesh.FaceCentre(faceId);
          Scalar boundDist = (polyMesh.CellCentre(face.OwnerId()) - polyMesh.FaceCentre(faceId)).norm();
          Scalar valueBound = boundaryFields[boundary.Id()].Value(faceCentre);

          coeffs(face.OwnerId(), face.OwnerId()) += (alphaFace * faceArea / boundDist) * theta * dt;
          constants[face.OwnerId()] -= (alphaFace * faceArea / boundDist) * (1.0 - theta) * dt * x0[face.OwnerId()];
          constants[face.OwnerId()] += valueBound * (alphaFace * faceArea / boundDist) * dt;
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
