#ifndef CFD_BASICS_LINALG_TRIDIAG_H_
#define CFD_BASICS_LINALG_TRIDIAG_H_

#include "nuenv/core"

namespace cfd_basics {
template<typename Scalar>
class TridiagMatrix {
 public:
  nuenv::Index Dimension();

  void SetVectors(const nuenv::VectorX<Scalar>& a,
                  const nuenv::VectorX<Scalar>& b,
                  const nuenv::VectorX<Scalar>& c,
                  const nuenv::VectorX<Scalar>& d);

  void SetMatrix(nuenv::MatrixSQX<Scalar>& m,
                 const nuenv::VectorX<Scalar>& constant);

  nuenv::MatrixSQX<Scalar> GetDense();

  nuenv::VectorX<Scalar> Solve();

 private:
  nuenv::Index dimension_;
  nuenv::VectorX<Scalar> a_, b_, c_, d_;

};

template<typename Scalar>
void TridiagMatrix<Scalar>::SetVectors(const nuenv::VectorX<Scalar>& a,
                                       const nuenv::VectorX<Scalar>& b,
                                       const nuenv::VectorX<Scalar>& c,
                                       const nuenv::VectorX<Scalar>& d) {
//  if (a.size() - 1 != b.size() != c.size() - 1 != d.size()) {
//    throw;
//  }

  a_ = a;
  b_ = b;
  c_ = c;
  d_ = d;

  dimension_ = d_.size();
}

template<typename Scalar>
nuenv::Index TridiagMatrix<Scalar>::Dimension() {
  return dimension_;
}

template<typename Scalar>
void TridiagMatrix<Scalar>::SetMatrix(nuenv::MatrixSQX<Scalar>& m,
                                      const nuenv::VectorX<Scalar>& constant) {
  a_ = m.getDiagonal(-1);
  b_ = m.getDiagonal(0);
  c_ = m.getDiagonal(1);
  d_ = constant;

  if (a_.size() - 1 != b_.size() != c_.size() - 1
      != d_.size()) {
    throw;
  }

  dimension_ = d_.size();
}

template<typename Scalar>
nuenv::MatrixSQX<Scalar> TridiagMatrix<Scalar>::GetDense() {
  nuenv::MatrixSQX<Scalar> matrix(dimension_, dimension_, 0);

  matrix(dimension_ - 1, dimension_ - 1) = b_(dimension_ - 1);

  for (nuenv::Index i = 0; i < dimension_ - 1; i++) {
    matrix(i, i) = b_(i);
    matrix(i + 1, i) = a_(i + 1);
    matrix(i, i + 1) = c_(i);
  }

  return matrix;
}

template<typename Scalar>
nuenv::VectorX<Scalar> TridiagMatrix<Scalar>::Solve() {
  if (dimension_ == 0) { throw; }

  nuenv::VectorX<Scalar> p = nuenv::VectorX<Scalar>::Zero(dimension_);
  nuenv::VectorX<Scalar> q = nuenv::VectorX<Scalar>::Zero(dimension_);
  nuenv::VectorX<Scalar> t = nuenv::VectorX<Scalar>::Zero(dimension_);

  Scalar aux;

  p(0) = c_(0) / b_(0);
  q(0) = d_(0) / b_(0);

  for (nuenv::Index i = 1; i < dimension_; i++) {
    aux = 1.0 / (b_(i) - a_(i) * p(i - 1));
    p(i) = c_(i) * aux;
    q(i) = (d_(i) - a_(i) * q(i - 1)) * aux;
  }

  t(dimension_ - 1) = q(dimension_ - 1);
  for (nuenv::Index i = dimension_ - 1; i > 0; i--) {
    t(i - 1) = q(i - 1) - p(i - 1) * t(i);
  }

  return t;
}

} // cfd_basics

#endif //CFD_BASICS_LINALG_TRIDIAG_H_
