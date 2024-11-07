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

  nuenv::Index nVolumes = 10;

  auto mesh = MakePolyMesh1D<Scalar, ScalarField>(0.1, 1.0, nVolumes, 4);
  auto boundaries = MakePolyBoundaries1D<Scalar>(nVolumes, 4);

  nuenv::VectorX<BoundaryField<Scalar>> boundaryFields(3);
  boundaryFields[0] = BoundaryField<Scalar>(0, 1, 0.0);  // boundaryFieldLeft
  boundaryFields[1] = BoundaryField<Scalar>(0, 1, 0.0);  // boundaryFieldRight
  boundaryFields[2] = BoundaryField<Scalar>(0, 0, 0.0);  // boundaryFieldEmpty

  PropertiesList<Scalar> properties(map<Property, Scalar> {
      {Property::kCondutivity, 1.0},
      {Property::kHeatSource, 0.0},
      {Property::kThermalDiffusivity, 1.17e-4},
  });

  nuenv::VectorX<Scalar> xp = nuenv::VectorX<Scalar>::Zero(nVolumes + 2);
  xp[0] = 0.0;
  xp[xp.size() - 1] = 0.1;
  for (nuenv::Index i = 1; i < xp.size() - 1; ++i) {
    xp[i] +=
        0.1 * static_cast<Scalar>(i - 1) / static_cast<Scalar>(nVolumes) + 0.5 * 0.1 / static_cast<Scalar>(nVolumes);
  }

  nuenv::VectorX<Scalar> x0 = nuenv::LinearSpace(0.0, 0.1 - 0.1 / static_cast<Scalar>(nVolumes), nVolumes);;
  for (Scalar& x : x0) {
    x += 0.5 * 0.1 / static_cast<Scalar>(nVolumes);
    x = sin(x * nuenv::pi / 0.1);
  }

  auto system = BuilderDiffusion(mesh, boundaries, boundaryFields, properties, x0, 4.0, 0.5);
  auto sol = GaussSeidel<Scalar>(system.Coeffs, system.Constants, 1e-16, 10000);
  nuenv::VectorX<Scalar> sol_aux = nuenv::VectorX<Scalar>::Zero(nVolumes + 2);
  sol_aux[0] = 0.0;
  sol_aux[sol_aux.size() - 1] = 0.0;
  for (nuenv::Index i = 1; i < sol_aux.size() - 1; ++i) {
    sol_aux[i] = sol.x[i - 1];
  }
  cout << "================================================================================" << endl;
  cout << "T_avg = " << TrapezoidalIntegral(xp, sol_aux) << endl;

  system = BuilderDiffusion(mesh, boundaries, boundaryFields, properties, sol.x, 4.0, 0.5);
  sol = GaussSeidel<Scalar>(system.Coeffs, system.Constants, 1e-16, 10000);
  sol_aux[0] = 0.0;
  sol_aux[sol_aux.size() - 1] = 0.0;
  for (nuenv::Index i = 1; i < sol_aux.size() - 1; ++i) {
    sol_aux[i] = sol.x[i - 1];
  }
  cout << "================================================================================" << endl;
  cout << "T_avg = " << TrapezoidalIntegral(xp, sol_aux) << endl;

  system = BuilderDiffusion(mesh, boundaries, boundaryFields, properties, sol.x, 4.0, 0.5);
  sol = GaussSeidel<Scalar>(system.Coeffs, system.Constants, 1e-16, 10000);
  sol_aux[0] = 0.0;
  sol_aux[sol_aux.size() - 1] = 0.0;
  for (nuenv::Index i = 1; i < sol_aux.size() - 1; ++i) {
    sol_aux[i] = sol.x[i - 1];
  }
  cout << "================================================================================" << endl;
  cout << "T_avg = " << TrapezoidalIntegral(xp, sol_aux) << endl;

  system = BuilderDiffusion(mesh, boundaries, boundaryFields, properties, sol.x, 4.0, 0.5);
  sol = GaussSeidel<Scalar>(system.Coeffs, system.Constants, 1e-16, 10000);
  sol_aux[0] = 0.0;
  sol_aux[sol_aux.size() - 1] = 0.0;
  for (nuenv::Index i = 1; i < sol_aux.size() - 1; ++i) {
    sol_aux[i] = sol.x[i - 1];
  }
  cout << "================================================================================" << endl;
  cout << "T_avg = " << TrapezoidalIntegral(xp, sol_aux) << endl;

  system = BuilderDiffusion(mesh, boundaries, boundaryFields, properties, sol.x, 4.0, 0.5);
  sol = GaussSeidel<Scalar>(system.Coeffs, system.Constants, 1e-16, 10000);
  sol_aux[0] = 0.0;
  sol_aux[sol_aux.size() - 1] = 0.0;
  for (nuenv::Index i = 1; i < sol_aux.size() - 1; ++i) {
    sol_aux[i] = sol.x[i - 1];
  }
  cout << "================================================================================" << endl;
  cout << "T_avg = " << TrapezoidalIntegral(xp, sol_aux) << endl;

  cout << "================================================================================" << endl;
  cout << system.Coeffs << endl;

  cout << "================================================================================" << endl;
  cout << system.Constants << endl;

  cout << "================================================================================" << endl;
  cout << sol.x << endl;

  return 0;
}
