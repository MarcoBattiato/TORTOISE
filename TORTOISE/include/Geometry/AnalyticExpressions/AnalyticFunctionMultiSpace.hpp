// Copyright © 2019, Marco Battiato <marco.battiato@ntu.edu.sg; battiato.marco@gmail.com>, All rights reserved.
//
// Licensed under the GNU GENERAL PUBLIC LICENSE Version 3 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   https://www.gnu.org/licenses/
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
// implied. See the License for the specific language governing
// permissions and limitations under the License.
//
//
// DISCLAIMER: This is a version under active development and testing.
// Not all features have been sufficiently tested, and several features
// are only partially implemented. Moreover both the core and the interface
// may change. You are discouraged from using this version to publish
// results without the supervision of the developer.
//
// If you want, you are welcome to act as a beta tester. In that case
// please contact Marco Battiato at:
// marco.battiato@ntu.edu.sg or battiato.marco@gmail.com
//
// Check if a newer, stable and tested version has, in the meanwhile,
// been made available at:
// https://github.com/MarcoBattiato/TORTOISE
//
//
//  AnalyticFunctionMultiSpace.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 24/9/21.
//

#ifndef AnalyticFunctionMultiSpace_hpp
#define AnalyticFunctionMultiSpace_hpp

#include <Geometry/GeometryCore/Geometry.hpp>

#include <cassert>
#include <vector>
#include <iostream>

namespace Tortoise {


// ============
// AnalyticFunction base class
// ============
template <typename AnFunc> class AnalyticFunctionMultiSpace{
public:
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return static_cast<AnFunc const&>(*this)(arrayOfPoints);
    }
};

// ============
// Most atomic analytic functions
// ============
class AnalyticFunctionMultiSpacePoint : public AnalyticFunctionMultiSpace<AnalyticFunctionMultiSpacePoint> {
    const int _which;  // which point
public:
    AnalyticFunctionMultiSpacePoint(int which): _which(which){}
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        std::cout << ">" << _which << "<";
        return arrayOfPoints.col(_which);
    }
};
class AnalyticFunctionMultiSpacePointCoordinate : public AnalyticFunctionMultiSpace<AnalyticFunctionMultiSpacePointCoordinate> {
    const int _which;  // which point
    const int _dir;    // direction
public:
    AnalyticFunctionMultiSpacePointCoordinate(int which, int dir): _which(which), _dir(dir){}
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return arrayOfPoints(_dir,_which);
    }
};


// ============
// Unary Operations
// ============
template <typename AnFunc> class AnalyticFunctionMultiSpaceNegative :
public AnalyticFunctionMultiSpace<AnalyticFunctionMultiSpaceNegative<AnFunc>> {
    const AnFunc& _u;
public:
    AnalyticFunctionMultiSpaceNegative(const AnFunc& u): _u(u){}
    template <typename Derived> auto
    operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return -_u(arrayOfPoints);
    }
};
template <typename AnFunc> AnalyticFunctionMultiSpaceNegative<AnFunc>
operator-(const AnalyticFunctionMultiSpace<AnFunc>& u) {
   return AnalyticFunctionMultiSpaceNegative<AnFunc>(*static_cast<const AnFunc*>(&u));
}

// ============
// Binary Operations
// ============
#define ANALYTICFUNCTIONMULTISPACE_OPERATOR_OVERLOAD(_NAME_CLASS_,_OPER_)               \
template <typename AnFunc1, typename AnFunc2> class _NAME_CLASS_ :                      \
public AnalyticFunctionMultiSpace<_NAME_CLASS_<AnFunc1,AnFunc2>> {                      \
    const AnFunc1& _u;                                                                  \
    const AnFunc2& _v;                                                                  \
public:                                                                                 \
    _NAME_CLASS_(const AnFunc1& u, const AnFunc2& v): _u(u), _v(v) {}                   \
    template <typename Derived> auto operator()                                         \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                           \
        return _u(arrayOfPoints) _OPER_ _v(arrayOfPoints);                              \
    }                                                                                   \
};                                                                                      \
template <typename AnFunc1, typename AnFunc2> _NAME_CLASS_<AnFunc1, AnFunc2>            \
operator _OPER_(const AnalyticFunctionMultiSpace<AnFunc1>& u,                           \
            const AnalyticFunctionMultiSpace<AnFunc2>& v) {                             \
   return _NAME_CLASS_<AnFunc1, AnFunc2>(*static_cast<const AnFunc1*>(&u),              \
                                         *static_cast<const AnFunc2*>(&v));             \
}

ANALYTICFUNCTIONMULTISPACE_OPERATOR_OVERLOAD(AnalyticFunctionMultiSpaceSum, +)
ANALYTICFUNCTIONMULTISPACE_OPERATOR_OVERLOAD(AnalyticFunctionMultiSpaceSub, -)
ANALYTICFUNCTIONMULTISPACE_OPERATOR_OVERLOAD(AnalyticFunctionMultiSpaceMult,*)
ANALYTICFUNCTIONMULTISPACE_OPERATOR_OVERLOAD(AnalyticFunctionMultiSpaceDiv, /)

// ============
// Operations with scalars (symmetric)
// ============
#define ANALYTICFUNCTIONMULTISPACE_OPERATOR_SCALAR_OVERLOAD(_NAME_CLASS_,_OPER_)        \
template <typename AnFunc> class _NAME_CLASS_ :                                         \
public AnalyticFunctionMultiSpace<_NAME_CLASS_<AnFunc>> {                               \
    const AnFunc& _u;                                                                   \
    const Real _v;                                                             \
public:                                                                                 \
    _NAME_CLASS_(const AnFunc& u, const Real v): _u(u), _v(v) {}               \
    template <typename Derived> auto operator()                                         \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                           \
        return _u(arrayOfPoints) _OPER_ _v;                                             \
    }                                                                                   \
};                                                                                      \
template <typename AnFunc> _NAME_CLASS_<AnFunc>                                         \
operator _OPER_(const AnalyticFunctionMultiSpace<AnFunc>& u,                            \
                const Real v) {                                                \
   return _NAME_CLASS_<AnFunc>(*static_cast<const AnFunc*>(&u), v);                     \
}                                                                                       \
template <typename AnFunc> _NAME_CLASS_<AnFunc>                                         \
operator _OPER_(const Real v,                                                  \
                const AnalyticFunctionMultiSpace<AnFunc>& u) {                          \
   return _NAME_CLASS_<AnFunc>(*static_cast<const AnFunc*>(&u), v);                     \
}

ANALYTICFUNCTIONMULTISPACE_OPERATOR_SCALAR_OVERLOAD(AnalyticFunctionMultiSpaceSumScal,+)
ANALYTICFUNCTIONMULTISPACE_OPERATOR_SCALAR_OVERLOAD(AnalyticFunctionMultiSpaceMultScal,*)

// ============
// Operations with scalars (asymmetric)
// ============
#define ANALYTICFUNCTIONMULTISPACE_OPERATOR_SCALARASYMMETRIC_OVERLOAD(_NAME_CLASS_,_OPER_)  \
template <typename AnFunc> class _NAME_CLASS_ :                                         \
public AnalyticFunctionMultiSpace<_NAME_CLASS_<AnFunc>> {                               \
    const AnFunc& _u;                                                                   \
    const Real _v;                                                             \
public:                                                                                 \
    _NAME_CLASS_(const AnFunc& u, const Real v): _u(u), _v(v) {}               \
    template <typename Derived> auto operator()                                         \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                           \
        return _u(arrayOfPoints) _OPER_ _v;                                             \
    }                                                                                   \
};                                                                                      \
template <typename AnFunc> _NAME_CLASS_<AnFunc>                                         \
operator _OPER_(const AnalyticFunctionMultiSpace<AnFunc>& u,                            \
                const Real v) {                                                \
   return _NAME_CLASS_<AnFunc>(*static_cast<const AnFunc*>(&u), v);                     \
}                                                                                       \
template <typename AnFunc> class _NAME_CLASS_ ## 1 :                                    \
public AnalyticFunctionMultiSpace<_NAME_CLASS_ ## 1<AnFunc>> {                          \
    const AnFunc& _u;                                                                   \
    const Real _v;                                                             \
public:                                                                                 \
    _NAME_CLASS_ ## 1(const AnFunc& u, const Real v): _u(u), _v(v) {}          \
    template <typename Derived> auto operator()                                         \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                           \
        return _v _OPER_ _u(arrayOfPoints);                                             \
    }                                                                                   \
};                                                                                      \
template <typename AnFunc> _NAME_CLASS_ ## 1<AnFunc>                                    \
operator _OPER_(const Real v,                                                  \
                const AnalyticFunctionMultiSpace<AnFunc>& u) {                          \
   return _NAME_CLASS_ ## 1<AnFunc>(*static_cast<const AnFunc*>(&u), v);                \
}

ANALYTICFUNCTIONMULTISPACE_OPERATOR_SCALARASYMMETRIC_OVERLOAD(AnalyticFunctionMultiSpaceSubScal,-)
ANALYTICFUNCTIONMULTISPACE_OPERATOR_SCALARASYMMETRIC_OVERLOAD(AnalyticFunctionMultiSpaceDivScal,/)

// ============
// Operations with vectors (symmetric)
// ============
#define ANALYTICFUNCTIONMULTISPACE_OPERATOR_VECTOR_OVERLOAD(_NAME_CLASS_,_OPER_)        \
template <typename AnFunc> class _NAME_CLASS_ :                                         \
public AnalyticFunctionMultiSpace<_NAME_CLASS_<AnFunc>> {                               \
    const AnFunc& _u;                                                                   \
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> _v;              \
public:                                                                                 \
    template <typename Derived> _NAME_CLASS_(                                           \
        const AnFunc& u, const Eigen::MatrixBase<Derived>& v): _u(u), _v(v) {}          \
    template <typename Derived> auto operator()                                         \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                           \
        return _u(arrayOfPoints) _OPER_ _v;                                             \
    }                                                                                   \
};                                                                                      \
template <typename AnFunc, typename Derived> _NAME_CLASS_<AnFunc>                       \
    operator _OPER_(const AnalyticFunctionMultiSpace<AnFunc>& u,                        \
                    const Eigen::MatrixBase<Derived>& v) {                              \
   return _NAME_CLASS_<AnFunc>(*static_cast<const AnFunc*>(&u), v);                     \
}                                                                                       \
template <typename AnFunc, typename Derived> _NAME_CLASS_<AnFunc>                       \
    operator _OPER_(const Eigen::MatrixBase<Derived>& v,                                \
                const AnalyticFunctionMultiSpace<AnFunc>& u) {                          \
   return _NAME_CLASS_<AnFunc>(*static_cast<const AnFunc*>(&u), v);                     \
}

ANALYTICFUNCTIONMULTISPACE_OPERATOR_VECTOR_OVERLOAD(AnalyticFunctionMultiSpaceSumVect,+)
ANALYTICFUNCTIONMULTISPACE_OPERATOR_VECTOR_OVERLOAD(AnalyticFunctionMultiSpaceMultVect,*)

// ============
// Operations with vectors (asymmetric)
// ============
#define ANALYTICFUNCTIONMULTISPACE_OPERATOR_VECTORASYMMETRIC_OVERLOAD(_NAME_CLASS_,_OPER_)  \
template <typename AnFunc> class _NAME_CLASS_ :                                         \
public AnalyticFunctionMultiSpace<_NAME_CLASS_<AnFunc>> {                               \
    const AnFunc& _u;                                                                   \
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> _v;              \
public:                                                                                 \
    template <typename Derived>_NAME_CLASS_                                             \
        (const AnFunc& u, const Eigen::MatrixBase<Derived>& v): _u(u), _v(v) {}         \
    template <typename Derived> auto operator()                                         \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                           \
        return _u(arrayOfPoints) _OPER_ _v;                                             \
    }                                                                                   \
};                                                                                      \
template <typename AnFunc, typename Derived> _NAME_CLASS_<AnFunc>                       \
operator _OPER_(const AnalyticFunctionMultiSpace<AnFunc>& u,                            \
                const Eigen::MatrixBase<Derived>& v) {                                  \
   return _NAME_CLASS_<AnFunc>(*static_cast<const AnFunc*>(&u), v);                     \
}                                                                                       \
template <typename AnFunc> class _NAME_CLASS_ ## 1 :                                    \
public AnalyticFunctionMultiSpace<_NAME_CLASS_ ## 1<AnFunc>> {                          \
    const AnFunc& _u;                                                                   \
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> _v;              \
public:                                                                                 \
    template <typename Derived> _NAME_CLASS_ ## 1                                       \
        (const AnFunc& u, const Eigen::MatrixBase<Derived>& v): _u(u), _v(v) {}         \
    template <typename Derived> auto operator()                                         \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                           \
        return _v _OPER_ _u(arrayOfPoints);                                             \
    }                                                                                   \
};                                                                                      \
template <typename AnFunc, typename Derived> _NAME_CLASS_ ## 1<AnFunc>                  \
operator _OPER_(const Eigen::MatrixBase<Derived>& v,                                    \
                const AnalyticFunctionMultiSpace<AnFunc>& u) {                          \
   return _NAME_CLASS_ ## 1<AnFunc>(*static_cast<const AnFunc*>(&u), v);                \
}

ANALYTICFUNCTIONMULTISPACE_OPERATOR_VECTORASYMMETRIC_OVERLOAD(AnalyticFunctionMultiSpaceSubVect,-)
ANALYTICFUNCTIONMULTISPACE_OPERATOR_VECTORASYMMETRIC_OVERLOAD(AnalyticFunctionMultiSpaceDivVect,/)

// ============
// Function of one Analytic Function
// ============
#define ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(_NAME_FUN_)                  \
template <typename AnFunc> class AnalyticFunctionMultiSpaceComp ## _NAME_FUN_ :             \
public AnalyticFunctionMultiSpace<AnalyticFunctionMultiSpaceComp ## _NAME_FUN_<AnFunc>> {   \
    const AnFunc& _u;                                                                       \
static auto internalFunc(const Real u){                                            \
            return std::_NAME_FUN_(u);};                                                    \
template<class Derived> static auto internalFunc(const Derived& vec){                       \
            return _NAME_FUN_(vec.array()).matrix();};                                      \
public:                                                                                     \
    AnalyticFunctionMultiSpaceComp ## _NAME_FUN_(const AnFunc& u): _u(u) {}                 \
    template <typename Derived> auto operator()                                             \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                               \
        return internalFunc(_u(arrayOfPoints));                                             \
    }                                                                                       \
};                                                                                          \
template <typename AnFunc> AnalyticFunctionMultiSpaceComp ## _NAME_FUN_<AnFunc>             \
        _NAME_FUN_(const AnalyticFunctionMultiSpace<AnFunc>& u){                            \
        return AnalyticFunctionMultiSpaceComp ## _NAME_FUN_<AnFunc>(                        \
            *static_cast<const AnFunc*>(&u)); }                                             \

ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(abs)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(exp)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(log)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(sqrt)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(sin)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(cos)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(tan)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(asin)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(acos)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(atan)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(sinh)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(cosh)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(tanh)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(asinh)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(acosh)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(atanh)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(ceil)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(floor)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(round)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(erf)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(erfc)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLE_OVERLOAD(lgamma)

// ============
// Function of one Analytic Function (vector only)
// ============
#define ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLEVEC_OVERLOAD(_NAME_FUN_)               \
template <typename AnFunc> class AnalyticFunctionMultiSpaceComp ## _NAME_FUN_ :             \
public AnalyticFunctionMultiSpace<AnalyticFunctionMultiSpaceComp ## _NAME_FUN_<AnFunc>> {   \
    const AnFunc& _u;                                                                       \
public:                                                                                     \
    AnalyticFunctionMultiSpaceComp ## _NAME_FUN_(const AnFunc& u): _u(u) {}                 \
    template <typename Derived> auto operator()                                             \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                               \
        return _u(arrayOfPoints)._NAME_FUN_();                                              \
    }                                                                                       \
};                                                                                          \
template <typename AnFunc> AnalyticFunctionMultiSpaceComp ## _NAME_FUN_<AnFunc>             \
        _NAME_FUN_(const AnalyticFunctionMultiSpace<AnFunc>& u){                            \
        return AnalyticFunctionMultiSpaceComp ## _NAME_FUN_<AnFunc>(                        \
            *static_cast<const AnFunc*>(&u)); }                                             \

//ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLEVEC_OVERLOAD(conj)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLEVEC_OVERLOAD(sum)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLEVEC_OVERLOAD(prod)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLEVEC_OVERLOAD(mean)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLEVEC_OVERLOAD(minCoeff)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLEVEC_OVERLOAD(maxCoeff)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLEVEC_OVERLOAD(transpose)
ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLEVEC_OVERLOAD(norm)

// ============
// Function of one Analytic Function with parameter
// ============
#define ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLEPARAM_OVERLOAD(_NAME_FUN_)                 \
template <typename AnFunc> class AnalyticFunctionMultiSpaceComp ## _NAME_FUN_ :                 \
public AnalyticFunctionMultiSpace<AnalyticFunctionMultiSpaceComp ## _NAME_FUN_<AnFunc>> {       \
    const AnFunc& _u;                                                                           \
    const Real _param;                                                                 \
static auto internalFunc(const Real u, const Real param){                     \
            return std::_NAME_FUN_(u, param);};                                                 \
template<class Derived> static auto internalFunc(const Derived& u, const Real param){  \
            return _NAME_FUN_(u.array(), param).matrix();};                                     \
public:                                                                                         \
    AnalyticFunctionMultiSpaceComp ## _NAME_FUN_                                                \
            (const AnFunc& u, const Real param): _u(u), _param(param) {}               \
    template <typename Derived> auto operator()                                                 \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                                   \
        return internalFunc(_u(arrayOfPoints), _param);                                         \
    }                                                                                           \
};                                                                                              \
template <typename AnFunc> AnalyticFunctionMultiSpaceComp ## _NAME_FUN_<AnFunc>                 \
        _NAME_FUN_(const AnalyticFunctionMultiSpace<AnFunc>& u, Real v){               \
        return AnalyticFunctionMultiSpaceComp ## _NAME_FUN_<AnFunc>(                            \
            *static_cast<const AnFunc*>(&u),v); }                                               \

ANALYTICFUNCTIONMULTISPACE_COMPOSITION_SINGLEPARAM_OVERLOAD(pow)

// ============
// Function of two Analytic Functions
// ============
template <typename AnFunc1, typename AnFunc2> class AnalyticFunctionMultiSpaceCompDot :
public AnalyticFunctionMultiSpace <AnalyticFunctionMultiSpaceCompDot<AnFunc1,AnFunc2>> {
    const AnFunc1& _u;
    const AnFunc2& _v;
public:
    AnalyticFunctionMultiSpaceCompDot(const AnFunc1& u, const AnFunc2& v): _u(u), _v(v) {}
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return _u(arrayOfPoints).dot(_v(arrayOfPoints));
    }
};
template <typename AnFunc1, typename AnFunc2> AnalyticFunctionMultiSpaceCompDot<AnFunc1, AnFunc2>
    dot(const AnalyticFunctionMultiSpace<AnFunc1>& u, const AnalyticFunctionMultiSpace<AnFunc2>& v) {
   return AnalyticFunctionMultiSpaceCompDot<AnFunc1, AnFunc2>(
        *static_cast<const AnFunc1*>(&u),*static_cast<const AnFunc2*>(&v));
}

template <typename AnFunc1> class AnalyticFunctionMultiSpaceCompDotConst :
public AnalyticFunctionMultiSpace <AnalyticFunctionMultiSpaceCompDotConst<AnFunc1>> {
    const AnFunc1& _u;
    const Eigen::Matrix<Real, Eigen::Dynamic, 1> _v;
public:
    template <typename Derived> AnalyticFunctionMultiSpaceCompDotConst(const AnFunc1& u, const Eigen::MatrixBase<Derived>& v): _u(u), _v(v) {}
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return _u(arrayOfPoints).dot(_v);
    }
};
template <typename AnFunc1, typename Derived> AnalyticFunctionMultiSpaceCompDotConst<AnFunc1>
    dot(const AnalyticFunctionMultiSpace<AnFunc1>& u, const Eigen::MatrixBase<Derived>& v) {
   return AnalyticFunctionMultiSpaceCompDotConst<AnFunc1>(
        *static_cast<const AnFunc1*>(&u),v);
}
template <typename AnFunc1, typename Derived> AnalyticFunctionMultiSpaceCompDotConst<AnFunc1>
    dot(const Eigen::MatrixBase<Derived>& v, const AnalyticFunctionMultiSpace<AnFunc1>& u) {
   return AnalyticFunctionMultiSpaceCompDotConst<AnFunc1>(
        *static_cast<const AnFunc1*>(&u),v);
}

template <typename AnFunc1, typename AnFunc2> class AnalyticFunctionMultiSpaceCompCross :
public AnalyticFunctionMultiSpace <AnalyticFunctionMultiSpaceCompCross<AnFunc1,AnFunc2>> {
    const AnFunc1& _u;
    const AnFunc2& _v;
public:
    AnalyticFunctionMultiSpaceCompCross(const AnFunc1& u, const AnFunc2& v): _u(u), _v(v) {}
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return (_u(arrayOfPoints).cross(_v(arrayOfPoints)));
    }
};
template <typename AnFunc1, typename AnFunc2> AnalyticFunctionMultiSpaceCompCross<AnFunc1, AnFunc2>
    cross(const AnalyticFunctionMultiSpace<AnFunc1>& u, const AnalyticFunctionMultiSpace<AnFunc2>& v) {
   return AnalyticFunctionMultiSpaceCompCross<AnFunc1, AnFunc2>(
        *static_cast<const AnFunc1*>(&u),*static_cast<const AnFunc2*>(&v));
}

template <typename AnFunc1> class AnalyticFunctionMultiSpaceCompCrossConst :
public AnalyticFunctionMultiSpace <AnalyticFunctionMultiSpaceCompCrossConst<AnFunc1>> {
    const AnFunc1& _u;
    const Eigen::Matrix<Real, Eigen::Dynamic, 1> _v;
public:
    template <typename Derived> AnalyticFunctionMultiSpaceCompCrossConst(const AnFunc1& u, const Eigen::MatrixBase<Derived>& v): _u(u), _v(v) {}
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return _u(arrayOfPoints).cross(_v);
    }
};
template <typename AnFunc1, typename Derived> AnalyticFunctionMultiSpaceCompCrossConst<AnFunc1>
    cross(const AnalyticFunctionMultiSpace<AnFunc1>& u, const Eigen::MatrixBase<Derived>& v) {
   return AnalyticFunctionMultiSpaceCompCrossConst<AnFunc1>(
        *static_cast<const AnFunc1*>(&u),v);
}
template <typename AnFunc1, typename Derived> AnalyticFunctionMultiSpaceCompCrossConst<AnFunc1>
    cross(const Eigen::MatrixBase<Derived>& v, const AnalyticFunctionMultiSpace<AnFunc1>& u) {
    return AnalyticFunctionMultiSpaceCompCrossConst<AnFunc1>(
        *static_cast<const AnFunc1*>(&u), -v);
}




// Initialisation of coordinates and points
extern AnalyticFunctionMultiSpacePoint k0;
extern AnalyticFunctionMultiSpacePoint k1;
extern AnalyticFunctionMultiSpacePoint k2;
extern AnalyticFunctionMultiSpacePoint k3;
extern AnalyticFunctionMultiSpacePoint k4;

extern AnalyticFunctionMultiSpacePointCoordinate k0x;
extern AnalyticFunctionMultiSpacePointCoordinate k0y;
extern AnalyticFunctionMultiSpacePointCoordinate k0z;
extern AnalyticFunctionMultiSpacePointCoordinate k1x;
extern AnalyticFunctionMultiSpacePointCoordinate k1y;
extern AnalyticFunctionMultiSpacePointCoordinate k1z;
extern AnalyticFunctionMultiSpacePointCoordinate k2x;
extern AnalyticFunctionMultiSpacePointCoordinate k2y;
extern AnalyticFunctionMultiSpacePointCoordinate k2z;
extern AnalyticFunctionMultiSpacePointCoordinate k3x;
extern AnalyticFunctionMultiSpacePointCoordinate k3y;
extern AnalyticFunctionMultiSpacePointCoordinate k3z;
extern AnalyticFunctionMultiSpacePointCoordinate k4x;
extern AnalyticFunctionMultiSpacePointCoordinate k4y;
extern AnalyticFunctionMultiSpacePointCoordinate k4z;


} // namespace Tortoise


#endif /* AnalyticFunctionMultiSpace_hpp */
