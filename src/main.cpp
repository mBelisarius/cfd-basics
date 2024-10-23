#include <iostream>

#include "nuenv/core"

#include "mesh/poly_mesh.hpp"

using namespace cfd_basics;

int main() {
  using std::cout;
  using std::endl;

  using Scalar = double;
  using ScalarField = nuenv::Vector3X<double>;

  nuenv::VectorX_s<ScalarField, 9> points {
      ScalarField {0.0, -0.36602540379, -0.36602540379},
      ScalarField {0.0, -0.36602540379, 0.36602540379},
      ScalarField {0.0, 0.5, 0.0},
      ScalarField {0.5, -0.36602540379, -0.36602540379},
      ScalarField {0.5, -0.36602540379, 0.36602540379},
      ScalarField {0.5, 0.5, 0.0},
      ScalarField {1.0, -0.36602540379, -0.36602540379},
      ScalarField {1.0, -0.36602540379, 0.36602540379},
      ScalarField {1.0, 0.5, 0.0},
  };

  auto face_w = Face(0, 0, -1, nuenv::VectorX_s<nuenv::Index, 3> {0, 1, 2});

  auto face_w1 = Face(1, 0, -1, nuenv::VectorX_s<nuenv::Index, 4> {0, 1, 4, 3});
  auto face_w2 = Face(2, 0, -1, nuenv::VectorX_s<nuenv::Index, 4> {1, 2, 5, 4});
  auto face_w3 = Face(3, 0, -1, nuenv::VectorX_s<nuenv::Index, 4> {2, 0, 3, 5});

  auto face_o = Face(4, 0, 1, nuenv::VectorX_s<nuenv::Index, 3> {3, 4, 5});

  auto face_e1 = Face(5, 1, -1, nuenv::VectorX_s<nuenv::Index, 4> {3, 4, 7, 6});
  auto face_e2 = Face(6, 1, -1, nuenv::VectorX_s<nuenv::Index, 4> {4, 5, 8, 7});
  auto face_e3 = Face(7, 1, -1, nuenv::VectorX_s<nuenv::Index, 4> {5, 3, 6, 8});

  auto face_e = Face(8, 1, -1, nuenv::VectorX_s<nuenv::Index, 3> {6, 7, 8});

  nuenv::VectorX_s<Face, 9> faces {face_w, face_w1, face_w2, face_w3, face_o,
                                   face_e1, face_e2, face_e3, face_e};

  auto mesh = PolyMesh<Scalar, ScalarField>(points, faces);

  cout << mesh.FaceNormal(0) << endl;
  cout << mesh.FaceNormal(1) << endl;
  cout << mesh.FaceNormal(2) << endl;

  return 0;
}
