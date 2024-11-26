#include <iomanip>
#include <iostream>
#include <map>

using std::cout;
using std::endl;
using std::map;

#include "nuenv/algorithm"
#include "nuenv/core"

#include "builder/boundary_field.hpp"
#include "builder/diffusion.hpp"
#include "builder/properties.hpp"
#include "linalg/gauss_seidel.hpp"
#include "scripts.hpp"

using namespace cfd_basics;

template<typename Scalar>
Scalar TrapezoidalIntegral(const nuenv::VectorX<Scalar>& x, const nuenv::VectorX<Scalar>& y) {
  Scalar integral = 0.0;

  for (nuenv::Index i = 0; i < x.size() - 1; ++i) {
    integral += 0.5 * (y[i + 1] + y[i]) * (x[i + 1] - x[i]);
  }

  return integral;
}

int main() {
  cout << std::scientific << std::setprecision(15);

  using Scalar = double;
  using ScalarField = nuenv::Vector3X<double>;

  nuenv::Index nVolumesX = 11;
  nuenv::Index nVolumesY = 11;

  PropertiesList<Scalar> properties(map<Property, Scalar> {
      {Property::kCondutivity, 1.0},
      {Property::kHeatSource, 0.0},
      {Property::kThermalDiffusivity, 1.0},
  });

  auto mesh = MakePolyMesh2dRegular<Scalar, ScalarField>(1.0, 1.0, nVolumesX, nVolumesY);
  auto boundaries = MakePolyBoundaries2dRegular<Scalar>(nVolumesX, nVolumesY);

  nuenv::VectorX<BoundaryField<Scalar, ScalarField>> boundaryFields(5);
  boundaryFields[0] = BoundaryField<Scalar, ScalarField>(0,
                                                         1,
                                                         [](const ScalarField& x) { return static_cast<Scalar>(0.0); });  // boundaryFieldDown
  boundaryFields[1] = BoundaryField<Scalar, ScalarField>(1,
                                                         1,
                                                         [](const ScalarField& x) { return static_cast<Scalar>(0.0); });  // boundaryFieldRight
  boundaryFields[2] = BoundaryField<Scalar, ScalarField>(2,
                                                         1,
                                                         [](const ScalarField& x) {
                                                           return static_cast<Scalar>(sin(nuenv::pi * x[0]));
                                                         });  // boundaryFieldUp
  boundaryFields[3] = BoundaryField<Scalar, ScalarField>(3,
                                                         1,
                                                         [](const ScalarField& x) { return static_cast<Scalar>(0.0); });  // boundaryFieldLeft
  boundaryFields[4] = BoundaryField<Scalar, ScalarField>(4,
                                                         0,
                                                         [](const ScalarField& x) { return static_cast<Scalar>(0.0); });  // boundaryFieldEmpty

  nuenv::VectorX<Scalar> x0 = nuenv::VectorX<Scalar>::Zero(nVolumesX * nVolumesY);

  auto system = BuilderDiffusion<Scalar, ScalarField>(mesh, boundaries, boundaryFields, properties, x0, 2.0, 0.5);
  auto sol = GaussSeidel<Scalar>(system.Coeffs, system.Constants, 1e-16, 1000);

  for (nuenv::Index i = 0; i < 1000; ++i) {
    system = BuilderDiffusion<Scalar, ScalarField>(mesh, boundaries, boundaryFields, properties, sol.x, 2.0, 0.5);
    sol = GaussSeidel<Scalar>(system.Coeffs, system.Constants, 1e-16, 1000);
  }

  cout << "================================================================================" << endl;
  cout << system.Coeffs << endl;

  cout << "================================================================================" << endl;
  cout << system.Constants << endl;

  cout << "================================================================================" << endl;
  cout << sol.x << endl;

  return 0;
}
