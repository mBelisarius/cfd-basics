#ifndef CFD_BASICS_GENERATOR_GENERATOR_H_
#define CFD_BASICS_GENERATOR_GENERATOR_H_

#include "nuenv/core"

#include "linalg/tridiag.hpp"
#include "mesh/mesh.hpp"

namespace cfd_basics {

template<typename Scalar>
class DiscreteSystem {
 public:
  void SetMesh(Mesh1D<Scalar> mesh);

  void SetBoundaries(Scalar left, Scalar right);

  void CalculateMesh();

  TridiagMatrix<Scalar>
  GenerateSystem(const nuenv::VectorX<Scalar>& dxp,
                 const nuenv::VectorX<Scalar>& ke,
                 const nuenv::VectorX<Scalar>& dxe);

 private:
  Mesh1D<Scalar> mesh_;
  Scalar boundaryLeft_;
  Scalar boundaryRight_;
};

template<typename Scalar>
void DiscreteSystem<Scalar>::SetMesh(Mesh1D<Scalar> mesh) {
  this->mesh_ = mesh;
}

template<typename Scalar>
void DiscreteSystem<Scalar>::SetBoundaries(Scalar left, Scalar right) {
  boundaryLeft_ = left;
  boundaryRight_ = right;
}

template<typename Scalar>
void DiscreteSystem<Scalar>::CalculateMesh() {}

template<typename Scalar>
TridiagMatrix<Scalar>
DiscreteSystem<Scalar>::GenerateSystem(const nuenv::VectorX<Scalar>& dxp,
                                       const nuenv::VectorX<Scalar>& ke,
                                       const nuenv::VectorX<Scalar>& dxe) {
  nuenv::Index n = dxp.size();

  nuenv::VectorX<Scalar> aw(n);
  nuenv::VectorX<Scalar> ae(n);
  nuenv::VectorX<Scalar> ap(n);
  nuenv::VectorX<Scalar> bp(n);

  aw[0] = 0.0;
  ae[0] = -ke[1] / dxe[1];
  ap[0] = -ae[0] + 2.0 * ke[1] / dxp[0];
  bp[0] = mesh_.qdot(0) * dxp[0] + (2.0 * ke[1] / dxp[0]) * boundaryLeft_;

  aw[n - 1] = -ke[n] / dxe[n];
  ae[n - 1] = 0.0;
  ap[n - 1] = -aw[n - 1] + 2.0 * ke[n] / dxp[n - 1];
  bp[n - 1] = mesh_.qdot(n - 1) * dxp[n - 1] + (2.0 * ke[n] / dxp[n - 1]) * boundaryRight_;

  for (nuenv::Index i = 1; i < n - 1; i++) {
    aw[i] = -ke[i] / dxe[i];
    ae[i] = -ke[i + 1] / dxe[i + 1];
    ap[i] = -aw[i] - ae[i];
    bp[i] = mesh_.qdot(i) * dxp[i];
  }

  for (int it = 0; it < n; it++) {
    std::printf("%21.14e %21.14e %21.14e %21.14e \n",
                aw[it], ap[it], ae[it], bp[it]);
  }

  TridiagMatrix<Scalar> M;
  M.SetVectors(aw, ap, ae, bp);

  return M;
}

template<typename Scalar>
Scalar Patankar(Scalar kP, Scalar dxP, Scalar kE, Scalar dxE) {
  Scalar ki = ((dxP + dxE) * kP * kE) / (dxP * kP + dxE * kP);

  return ki;
}

//    template<typename Scalar>
//    inline Scalar ltmh(Scalar kp, Scalar tempP, Scalar tempE)
//    {
//        // TODO
//        // http://servidor.demec.ufpr.br/CFD/artigos_congressos/2013_Carvalho_Marchi_Perstchi_COBEM_2013.pdf
//
//        throw exception("Not implemented");
//        Scalar tempA = (3 * tempP + tempE) / 4;
//        Scalar tempB = (tempP + 3 * tempE) / 4;
//
//        Scalar ke = 2 * ;
//
//        return ke;
//    }
//
//    template<typename Scalar>
//    inline Vector<Scalar> dk1d(Vector<Scalar> kp)
//    {
//        // TODO
//        // http://servidor.demec.ufpr.br/CFD/artigos_congressos/2013_Carvalho_Marchi_Perstchi_COBEM_2013.pdf
//
//        throw exception("Not implemented");
//
//        size_t n = kp.size() - 1;
//        Vector<Scalar> ke(n + 1);
//
//        Scalar a, b, c;
//
//        ke[0] = kp[0];  // ke[0] = k(Ta);
//        ke[n] = kp[n];  // ke[n] = k(Tb);
//        for (size_t i = 1; i < n; i++)
//        {
//            a = 2 / ((temp[i + 1] - temp[i]) * (k1 + k2))
//            b = -k2;
//            c = k1;
//            ke[i] = a * (b + c);
//        }
//
//        return ke;
//    }

} // cfd_basics

#endif //CFD_BASICS_GENERATOR_GENERATOR_H_
