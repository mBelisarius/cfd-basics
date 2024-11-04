#include <iostream>
#include <map>

using std::cout;
using std::endl;
using std::map;

#include "nuenv/core"

#include "builder/boundary_field.hpp"
#include "builder/diffusion.hpp"
#include "builder/properties.hpp"
#include "linalg/gauss_seidel.hpp"
#include "scripts.hpp"

using namespace cfd_basics;

int main() {
  using Scalar = double;
  using ScalarField = nuenv::Vector3X<double>;

  nuenv::Index nVolumes = 10;

  auto mesh = MakePolyMesh1D<Scalar, ScalarField>(1.0, 0.1, nVolumes, 4);
  auto boundaries = MakePolyBoundaries1D<Scalar>(nVolumes, 4);

  nuenv::VectorX<BoundaryField<Scalar>> boundaryFields(3);
  boundaryFields[0] = BoundaryField<Scalar>(0, 1, 0.0);   // boundaryFieldLeft
  boundaryFields[1] = BoundaryField<Scalar>(0, 1, 100.0); // boundaryFieldRight
  boundaryFields[2] = BoundaryField<Scalar>(0, 0, 0.0);   // boundaryFieldEmpty

  PropertiesList<Scalar> properties(map<Property, Scalar> {
      {Property::kCondutivity, 400.0},
      {Property::kHeatSource, 5.0e5}
  });

  auto system = BuilderDiffusion(mesh, boundaries, boundaryFields, properties);

  auto sol = GaussSeidel<Scalar>(system.Coeffs, system.Constants, 1e-16, 10000);

  cout << "================================================================================" << endl;
  cout << system.Coeffs << endl;

  cout << "================================================================================" << endl;
  cout << system.Constants << endl;

  cout << "================================================================================" << endl;
  cout << sol.x << endl;

  return 0;
}
