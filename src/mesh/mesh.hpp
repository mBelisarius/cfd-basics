#ifndef CFD_BASICS_MESH_MESH_H_
#define CFD_BASICS_MESH_MESH_H_

#include "nuenv/core"

namespace cfd_basics {

template<typename Scalar>
class Mesh1D {
 public:
  Mesh1D() = default;

  Mesh1D(const nuenv::VectorX<Scalar>& x, const nuenv::VectorX<Scalar>& k,
         const nuenv::VectorX<Scalar>& qdot);

  const nuenv::VectorX<Scalar>& x();

  Scalar x(nuenv::Index i);

  void set_x(const nuenv::VectorX<Scalar> x);

  const nuenv::VectorX<Scalar>& k();

  Scalar k(nuenv::Index i);

  void set_k(nuenv::VectorX<Scalar> k);

  const nuenv::VectorX<Scalar>& qdot();

  Scalar qdot(nuenv::Index i);

  void set_qdot(nuenv::VectorX<Scalar> qdot);

 private:
  nuenv::VectorX<Scalar> x_;
  nuenv::VectorX<Scalar> k_;
  nuenv::VectorX<Scalar> qdot_;
};

template<typename Scalar>
Mesh1D<Scalar>::Mesh1D(const nuenv::VectorX<Scalar>& x,
                       const nuenv::VectorX<Scalar>& k,
                       const nuenv::VectorX<Scalar>& qdot)
    : x_(x), k_(k), qdot_(qdot) {}

template<typename Scalar>
const nuenv::VectorX<Scalar>& Mesh1D<Scalar>::x() {
  return &x_;
}

template<typename Scalar>
Scalar Mesh1D<Scalar>::x(nuenv::Index i) {
  return x_[i];
}

template<typename Scalar>
void Mesh1D<Scalar>::set_x(const nuenv::VectorX<Scalar> x) {
  x_ = x;
}

template<typename Scalar>
const nuenv::VectorX<Scalar>& Mesh1D<Scalar>::k() {
  return &k_;
}

template<typename Scalar>
Scalar Mesh1D<Scalar>::k(nuenv::Index i) {
  return k_[i];
}

template<typename Scalar>
void Mesh1D<Scalar>::set_k(const nuenv::VectorX<Scalar> k) {
  k_ = k;
}

template<typename Scalar>
const nuenv::VectorX<Scalar>& Mesh1D<Scalar>::qdot() {
  return &qdot_;
}

template<typename Scalar>
Scalar Mesh1D<Scalar>::qdot(nuenv::Index i) {
  return qdot_[i];
}

template<typename Scalar>
void Mesh1D<Scalar>::set_qdot(const nuenv::VectorX<Scalar> qdot) {
  qdot_ = qdot;
}

} // cfd_basics

#endif //CFD_BASICS_MESH_MESH_H_
