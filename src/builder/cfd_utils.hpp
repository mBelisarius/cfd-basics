#ifndef CFD_BASICS_BUILDER_CFD_UTILS_HPP_
#define CFD_BASICS_BUILDER_CFD_UTILS_HPP_

#include "nuenv/core"

namespace cfd_basics {

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

#endif //CFD_BASICS_BUILDER_CFD_UTILS_HPP_
