#include <iostream>

using std::cout;
using std::endl;

#include "nuenv/core"

#include "builder/diffusion.hpp"
#include "scripts.hpp"

using namespace cfd_basics;

int main() {
  using Scalar = double;
  using ScalarField = nuenv::Vector3X<double>;

  nuenv::Index nVolumes = 5;

  auto mesh = MakePolyMesh1D<Scalar, ScalarField>(1.0, 1.0, nVolumes, 4);
  auto boundaries = MakePolyBoundaries1D<Scalar>(nVolumes, 4);

  auto system = BuilderDiffusion(mesh, boundaries);

  cout << system.Coeffs << endl;

  return 0;
}
