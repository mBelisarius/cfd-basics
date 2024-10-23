#ifndef CFD_BASICS_GENERATOR_DIFFUSION_H_
#define CFD_BASICS_GENERATOR_DIFFUSION_H_


#include "nuenv/core"

#include "linalg/gauss_seidel.hpp"
#include "mesh/poly_mesh.hpp"
#include "boundary_field.hpp"
#include "solution.hpp"

namespace cfd_basics {
template<typename Scalar, typename ScalarField>
class Diffusion2D {
 public:
  Diffusion2D() = default;

  Diffusion2D(const nuenv::MatrixSQX<Scalar>& a,
              const nuenv::VectorX<Scalar>& b);

  Diffusion2D(PolyMesh<Scalar, ScalarField> mesh,
              nuenv::VectorT<BoundaryField<Scalar>> boundaryField);

  nuenv::MatrixSQX<Scalar> a();

  nuenv::VectorX<Scalar> b();

  Solution<Scalar> sol();

  void generate(Scalar k, Scalar qdot);

  Solution<Scalar>
  solve(nuenv::VectorX<Scalar> x0, Scalar rtol, nuenv::Index maxIter);

 private:
  PolyMesh<Scalar, ScalarField> mesh_;
  nuenv::VectorT<BoundaryField<Scalar>> boundaryField_;
  nuenv::Index nt_;
  nuenv::MatrixSQX<Scalar> a_;
  nuenv::VectorX<Scalar> b_;
  Solution<Scalar> sol_;
};

template<typename Scalar, typename ScalarField>
Diffusion2D<Scalar, ScalarField>::Diffusion2D(const Matrix<Scalar>& a,
                                              const Vector<Scalar>& b)
    : a_(a), b_(b) {}

template<typename Scalar, typename ScalarField>
Diffusion2D<Scalar, ScalarField>::Diffusion2D(PolyMesh<Scalar, ScalarField> mesh,
                                              VectorSTL<BoundaryField<Scalar>> boundaryField)
    : mesh_(mesh), boundaryField_(boundaryField) {
  nt_ = mesh_.get_nFaces();
  for (auto& boundary : mesh_.get_boundaries()) {
    // TODO: index m_mesh.get_boundaries() instead
    if (boundary.index != 0) {
      nt_ += boundary.nFaces;
    }
  }

  a_ = Matrix<Scalar>::Zero(nt_, nt_);
  b_ = Vector<Scalar>::Zero(nt_);
}

template<typename Scalar, typename ScalarField>
Matrix<Scalar> Diffusion2D<Scalar, ScalarField>::a() {
  return a_;
}

template<typename Scalar, typename ScalarField>
Vector<Scalar> Diffusion2D<Scalar, ScalarField>::b() {
  return b_;
}

template<typename Scalar, typename ScalarField>
Solution<Scalar> Diffusion2D<Scalar, ScalarField>::sol() {
  return sol_;
}

template<typename Scalar, typename ScalarField>
void Diffusion2D<Scalar, ScalarField>::generate(Scalar k, Scalar qdot) {
  for (auto& face : mesh_.get_faces()) {
    b_(face.index) = qdot / k;
    for (auto& neighbourFace : mesh_.get_neighbours(face)) {
      Scalar value = k;
      a_(face.index, neighbourFace.index) = -value;
      a_(face.index, face.index) += value;
    }
  }

  nuenv::Index boundaryFaceIndex = mesh_.get_nFaces();

  for (auto& boundary : mesh_.get_boundaries()) {
    // TODO: index m_mesh.get_boundaries() instead
    if (boundary.index != 0) {
      for (auto& face : boundary.faces) {
        a_(face.index, boundaryFaceIndex) = -k;
        a_(boundaryFaceIndex, boundaryFaceIndex) = k;
        a_(boundaryFaceIndex, face.index) = k;
        b_(boundaryFaceIndex) = boundaryField_[boundary.index - 1]
            .condition(face.center);
        boundaryFaceIndex++;
      }
    }
  }
}

template<typename Scalar, typename ScalarField>
Solution<Scalar> Diffusion2D<Scalar, ScalarField>::solve(nuenv::VectorX<Scalar> x0,
                                                         Scalar rtol,
                                                         nuenv::Index maxIter) {
  sol_ = GaussSeidel<Scalar>(a_, b_, x0, rtol, maxIter);
  return sol_;
}

} // cfd_basics

#endif /* CFD_BASICS_GENERATOR_DIFFUSION_H_ */
