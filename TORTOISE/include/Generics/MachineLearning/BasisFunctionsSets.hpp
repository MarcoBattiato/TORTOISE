//
//  BasisFunctionsSets.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 30/9/23.
//

#ifndef BasisFunctionsSets_hpp
#define BasisFunctionsSets_hpp

#include <Eigen/Dense>

namespace Tortoise {
namespace SparseFunctionals {

template <int order> struct TaylorBasis{
    template <typename Derived> static auto basisF(const Eigen::ArrayBase<Derived>& x){
        return x * TaylorBasis<order-1>::basisF(x);
    }
};
template <> struct TaylorBasis<0>{
    template <typename Derived> static auto basisF(const Eigen::ArrayBase<Derived>& x){
        return Eigen::ArrayBase<Derived>::Constant(x.rows(), x.cols(), 1);
    }
};
template <> struct TaylorBasis<1>{
    template <typename Derived> static auto basisF(const Eigen::ArrayBase<Derived>& x){
        return x.matrix().array();
    }
};

template <int order> struct FourierBasis{
    template <typename Derived> static auto basisF(const Eigen::ArrayBase<Derived>& x){
        return sin(3.14159265358979323846 * x.matrix().array()) * FourierBasis<order-1>::basisF(x);
    }
};

template <> struct FourierBasis<0>{
    template <typename Derived> static auto basisF(const Eigen::ArrayBase<Derived>& x){
        return Eigen::ArrayBase<Derived>::Constant(x.rows(), x.cols(), 1);
    }
};

} // namespace SparseFunctionals
} // namespace Tortoise

#endif /* BasisFunctionsSets_hpp */
