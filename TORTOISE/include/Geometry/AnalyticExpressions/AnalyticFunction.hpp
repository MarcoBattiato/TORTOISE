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
//  AnalyticFunction.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 22/9/21.
//
//  This class uses expression templates on functors to provide some functionalities to express and use analytic functions

#ifndef AnalyticFunction_hpp
#define AnalyticFunction_hpp

#include <Geometry/GeometryCore/Geometry.hpp>

#include <cassert>
#include <vector>


namespace Tortoise {

namespace AnalyticExpression {

// ============
// AnalyticFunction base class
// ============
template <typename AnFunc> class AnalyticFunction{
public:
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& variable) const {
        return static_cast<AnFunc const&>(*this)(variable);
    }
};

// ============
// Most atomic analytic functions
// ============
class AnalyticFunctionPoint : public AnalyticFunction<AnalyticFunctionPoint> {
public:
//    template<typename Derived> Point<Derived::RowsAtCompileTime> operator()(const Eigen::MatrixBase<Derived>& variable) const {
//        return variable;
//    }
    template<typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& variable) const {
        return variable.col(0);
    }
};
class AnalyticFunctionCoordinate : public AnalyticFunction<AnalyticFunctionCoordinate> {
    const int _dir;    // direction
public:
    AnalyticFunctionCoordinate(int dir): _dir(dir){}
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return arrayOfPoints(_dir);
    }
};

// ============
// Unary Operations
// ============
template <typename AnFunc> class AnalyticFunctionNegative :
public AnalyticFunction<AnalyticFunctionNegative<AnFunc>> {
    const AnFunc& _u;
public:
    AnalyticFunctionNegative(const AnFunc& u): _u(u){}
    template <typename Derived> auto
    operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return -_u(arrayOfPoints);
    }
};
template <typename AnFunc> AnalyticFunctionNegative<AnFunc>
operator-(const AnalyticFunction<AnFunc>& u) {
   return AnalyticFunctionNegative<AnFunc>(*static_cast<const AnFunc*>(&u));
}

// ============
// Binary Operations
// ============
#define ANALYTICFUNCTION_OPERATOR_OVERLOAD(_NAME_CLASS_,_OPER_)                         \
template <typename AnFunc1, typename AnFunc2> class _NAME_CLASS_ :                      \
public AnalyticFunction<_NAME_CLASS_<AnFunc1,AnFunc2>> {                                \
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
operator _OPER_(const AnalyticFunction<AnFunc1>& u,                                     \
            const AnalyticFunction<AnFunc2>& v) {                                       \
   return _NAME_CLASS_<AnFunc1, AnFunc2>(*static_cast<const AnFunc1*>(&u),              \
                                         *static_cast<const AnFunc2*>(&v));             \
}

ANALYTICFUNCTION_OPERATOR_OVERLOAD(AnalyticFunctionSum, +)
ANALYTICFUNCTION_OPERATOR_OVERLOAD(AnalyticFunctionSub, -)
ANALYTICFUNCTION_OPERATOR_OVERLOAD(AnalyticFunctionMult,*)
ANALYTICFUNCTION_OPERATOR_OVERLOAD(AnalyticFunctionDiv, /)

// ============
// Operations with scalars (symmetric)
// ============
#define ANALYTICFUNCTION_OPERATOR_SCALAR_OVERLOAD(_NAME_CLASS_,_OPER_)                  \
template <typename AnFunc> class _NAME_CLASS_ :                                         \
public AnalyticFunction<_NAME_CLASS_<AnFunc>> {                                         \
    const AnFunc& _u;                                                                   \
    const Real _v;                                                                      \
public:                                                                                 \
    _NAME_CLASS_(const AnFunc& u, const Real v): _u(u), _v(v) {}                        \
    template <typename Derived> auto operator()                                         \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                           \
        return _u(arrayOfPoints) _OPER_ _v;                                             \
    }                                                                                   \
};                                                                                      \
template <typename AnFunc> _NAME_CLASS_<AnFunc>                                         \
operator _OPER_(const AnalyticFunction<AnFunc>& u,                                      \
                const Real v) {                                                         \
   return _NAME_CLASS_<AnFunc>(*static_cast<const AnFunc*>(&u), v);                     \
}                                                                                       \
template <typename AnFunc> _NAME_CLASS_<AnFunc>                                         \
operator _OPER_(const Real v,                                                           \
                const AnalyticFunction<AnFunc>& u) {                                    \
   return _NAME_CLASS_<AnFunc>(*static_cast<const AnFunc*>(&u), v);                     \
}

ANALYTICFUNCTION_OPERATOR_SCALAR_OVERLOAD(AnalyticFunctionSumScal,+)
ANALYTICFUNCTION_OPERATOR_SCALAR_OVERLOAD(AnalyticFunctionMultScal,*)

// ============
// Operations with scalars (asymmetric)
// ============
#define ANALYTICFUNCTION_OPERATOR_SCALARASYMMETRIC_OVERLOAD(_NAME_CLASS_,_OPER_)        \
template <typename AnFunc> class _NAME_CLASS_ :                                         \
public AnalyticFunction<_NAME_CLASS_<AnFunc>> {                                         \
    const AnFunc& _u;                                                                   \
    const Real _v;                                                                      \
public:                                                                                 \
    _NAME_CLASS_(const AnFunc& u, const Real v): _u(u), _v(v) {}                        \
    template <typename Derived> auto operator()                                         \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                           \
        return _u(arrayOfPoints) _OPER_ _v;                                             \
    }                                                                                   \
};                                                                                      \
template <typename AnFunc> _NAME_CLASS_<AnFunc>                                         \
operator _OPER_(const AnalyticFunction<AnFunc>& u,                                      \
                const Real v) {                                                         \
   return _NAME_CLASS_<AnFunc>(*static_cast<const AnFunc*>(&u), v);                     \
}                                                                                       \
template <typename AnFunc> class _NAME_CLASS_ ## 1 :                                    \
public AnalyticFunction<_NAME_CLASS_ ## 1<AnFunc>> {                                    \
    const AnFunc& _u;                                                                   \
    const Real _v;                                                                      \
public:                                                                                 \
    _NAME_CLASS_ ## 1(const AnFunc& u, const Real v): _u(u), _v(v) {}                   \
    template <typename Derived> auto operator()                                         \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                           \
        return _v _OPER_ _u(arrayOfPoints);                                             \
    }                                                                                   \
};                                                                                      \
template <typename AnFunc> _NAME_CLASS_ ## 1<AnFunc>                                    \
operator _OPER_(const Real v,                                                           \
                const AnalyticFunction<AnFunc>& u) {                                    \
   return _NAME_CLASS_ ## 1<AnFunc>(*static_cast<const AnFunc*>(&u), v);                \
}

ANALYTICFUNCTION_OPERATOR_SCALARASYMMETRIC_OVERLOAD(AnalyticFunctionSubScal,-)
ANALYTICFUNCTION_OPERATOR_SCALARASYMMETRIC_OVERLOAD(AnalyticFunctionDivScal,/)

// ============
// Operations with vectors (symmetric)
// ============
#define ANALYTICFUNCTION_OPERATOR_VECTOR_OVERLOAD(_NAME_CLASS_,_OPER_)                  \
template <typename AnFunc> class _NAME_CLASS_ :                                         \
public AnalyticFunction<_NAME_CLASS_<AnFunc>> {                                         \
    const AnFunc& _u;                                                                   \
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> _v;                       \
public:                                                                                 \
    template <typename Derived> _NAME_CLASS_(                                           \
        const AnFunc& u, const Eigen::MatrixBase<Derived>& v): _u(u), _v(v) {}          \
    template <typename Derived> auto operator()                                         \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                           \
        return _u(arrayOfPoints) _OPER_ _v;                                             \
    }                                                                                   \
};                                                                                      \
template <typename AnFunc, typename Derived> _NAME_CLASS_<AnFunc>                       \
    operator _OPER_(const AnalyticFunction<AnFunc>& u,                                  \
                    const Eigen::MatrixBase<Derived>& v) {                              \
   return _NAME_CLASS_<AnFunc>(*static_cast<const AnFunc*>(&u), v);                     \
}                                                                                       \
template <typename AnFunc, typename Derived> _NAME_CLASS_<AnFunc>                       \
    operator _OPER_(const Eigen::MatrixBase<Derived>& v,                                \
                const AnalyticFunction<AnFunc>& u) {                                    \
   return _NAME_CLASS_<AnFunc>(*static_cast<const AnFunc*>(&u), v);                     \
}

ANALYTICFUNCTION_OPERATOR_VECTOR_OVERLOAD(AnalyticFunctionSumVect,+)
ANALYTICFUNCTION_OPERATOR_VECTOR_OVERLOAD(AnalyticFunctionMultVect,*)

// ============
// Operations with vectors (asymmetric)
// ============
#define ANALYTICFUNCTION_OPERATOR_VECTORASYMMETRIC_OVERLOAD(_NAME_CLASS_,_OPER_)        \
template <typename AnFunc> class _NAME_CLASS_ :                                         \
public AnalyticFunction<_NAME_CLASS_<AnFunc>> {                                         \
    const AnFunc& _u;                                                                   \
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> _v;                       \
public:                                                                                 \
    template <typename Derived>_NAME_CLASS_                                             \
        (const AnFunc& u, const Eigen::MatrixBase<Derived>& v): _u(u), _v(v) {}         \
    template <typename Derived> auto operator()                                         \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                           \
        return _u(arrayOfPoints) _OPER_ _v;                                             \
    }                                                                                   \
};                                                                                      \
template <typename AnFunc, typename Derived> _NAME_CLASS_<AnFunc>                       \
operator _OPER_(const AnalyticFunction<AnFunc>& u,                                      \
                const Eigen::MatrixBase<Derived>& v) {                                  \
   return _NAME_CLASS_<AnFunc>(*static_cast<const AnFunc*>(&u), v);                     \
}                                                                                       \
template <typename AnFunc> class _NAME_CLASS_ ## 1 :                                    \
public AnalyticFunction<_NAME_CLASS_ ## 1<AnFunc>> {                                    \
    const AnFunc& _u;                                                                   \
    const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> _v;                       \
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
                const AnalyticFunction<AnFunc>& u) {                                    \
   return _NAME_CLASS_ ## 1<AnFunc>(*static_cast<const AnFunc*>(&u), v);                \
}

ANALYTICFUNCTION_OPERATOR_VECTORASYMMETRIC_OVERLOAD(AnalyticFunctionSubVect,-)
ANALYTICFUNCTION_OPERATOR_VECTORASYMMETRIC_OVERLOAD(AnalyticFunctionDivVect,/)

// ============
// Function of one Analytic Function
// ============
#define ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(_NAME_FUN_)                            \
template <typename AnFunc> class AnalyticFunctionComp ## _NAME_FUN_ :                       \
public AnalyticFunction<AnalyticFunctionComp ## _NAME_FUN_<AnFunc>> {                       \
    const AnFunc& _u;                                                                       \
static auto internalFunc(const Real u){                                                     \
            return std::_NAME_FUN_(u);};                                                    \
template<class Derived> static auto internalFunc(const Derived& vec){                       \
            return _NAME_FUN_(vec.array()).matrix();};                                      \
public:                                                                                     \
    AnalyticFunctionComp ## _NAME_FUN_(const AnFunc& u): _u(u) {}                           \
    template <typename Derived> auto operator()                                             \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                               \
        return internalFunc(_u(arrayOfPoints));                                             \
    }                                                                                       \
};                                                                                          \
template <typename AnFunc> AnalyticFunctionComp ## _NAME_FUN_<AnFunc>                       \
        _NAME_FUN_(const AnalyticFunction<AnFunc>& u){                                      \
        return AnalyticFunctionComp ## _NAME_FUN_<AnFunc>(                                  \
            *static_cast<const AnFunc*>(&u)); }                                             \

ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(abs)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(exp)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(log)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(sqrt)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(sin)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(cos)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(tan)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(asin)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(acos)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(atan)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(sinh)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(cosh)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(tanh)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(asinh)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(acosh)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(atanh)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(ceil)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(floor)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(round)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(erf)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(erfc)
ANALYTICFUNCTION_COMPOSITION_SINGLE_OVERLOAD(lgamma)

// ============
// Function of one Analytic Function (vector only)
// ============
#define ANALYTICFUNCTION_COMPOSITION_SINGLEVEC_OVERLOAD(_NAME_FUN_)                         \
template <typename AnFunc> class AnalyticFunctionComp ## _NAME_FUN_ :                       \
public AnalyticFunction<AnalyticFunctionComp ## _NAME_FUN_<AnFunc>> {                       \
    const AnFunc& _u;                                                                       \
public:                                                                                     \
    AnalyticFunctionComp ## _NAME_FUN_(const AnFunc& u): _u(u) {}                           \
    template <typename Derived> auto operator()                                             \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                               \
        return _u(arrayOfPoints)._NAME_FUN_();                                              \
    }                                                                                       \
};                                                                                          \
template <typename AnFunc> AnalyticFunctionComp ## _NAME_FUN_<AnFunc>                       \
        _NAME_FUN_(const AnalyticFunction<AnFunc>& u){                                      \
        return AnalyticFunctionComp ## _NAME_FUN_<AnFunc>(                                  \
            *static_cast<const AnFunc*>(&u)); }                                             \

//ANALYTICFUNCTION_COMPOSITION_SINGLEVEC_OVERLOAD(conj)
ANALYTICFUNCTION_COMPOSITION_SINGLEVEC_OVERLOAD(sum)
ANALYTICFUNCTION_COMPOSITION_SINGLEVEC_OVERLOAD(prod)
ANALYTICFUNCTION_COMPOSITION_SINGLEVEC_OVERLOAD(mean)
ANALYTICFUNCTION_COMPOSITION_SINGLEVEC_OVERLOAD(minCoeff)
ANALYTICFUNCTION_COMPOSITION_SINGLEVEC_OVERLOAD(maxCoeff)
ANALYTICFUNCTION_COMPOSITION_SINGLEVEC_OVERLOAD(transpose)
ANALYTICFUNCTION_COMPOSITION_SINGLEVEC_OVERLOAD(norm)

// ============
// Function of one Analytic Function with parameter
// ============
#define ANALYTICFUNCTION_COMPOSITION_SINGLEPARAM_OVERLOAD(_NAME_FUN_)                           \
template <typename AnFunc> class AnalyticFunctionComp ## _NAME_FUN_ :                           \
public AnalyticFunction<AnalyticFunctionComp ## _NAME_FUN_<AnFunc>> {                           \
    const AnFunc& _u;                                                                           \
    const Real _param;                                                                          \
static auto internalFunc(const Real u, const Real param){                                       \
            return std::_NAME_FUN_(u, param);};                                                 \
template<class Derived> static auto internalFunc(const Derived& u, const Real param){           \
            return _NAME_FUN_(u.array(), param).matrix();};                                     \
public:                                                                                         \
    AnalyticFunctionComp ## _NAME_FUN_                                                          \
            (const AnFunc& u, const Real param): _u(u), _param(param) {}                        \
    template <typename Derived> auto operator()                                                 \
    (const Eigen::MatrixBase<Derived>& arrayOfPoints) const {                                   \
        return internalFunc(_u(arrayOfPoints), _param);                                         \
    }                                                                                           \
};                                                                                              \
template <typename AnFunc> AnalyticFunctionComp ## _NAME_FUN_<AnFunc>                           \
        _NAME_FUN_(const AnalyticFunction<AnFunc>& u, Real v){                                  \
        return AnalyticFunctionComp ## _NAME_FUN_<AnFunc>(                                      \
            *static_cast<const AnFunc*>(&u),v); }                                               \

ANALYTICFUNCTION_COMPOSITION_SINGLEPARAM_OVERLOAD(pow)

// ============
// Function of two Analytic Functions
// ============
template <typename AnFunc1, typename AnFunc2> class AnalyticFunctionCompDot :
public AnalyticFunction <AnalyticFunctionCompDot<AnFunc1,AnFunc2>> {
    const AnFunc1& _u;
    const AnFunc2& _v;
public:
    AnalyticFunctionCompDot(const AnFunc1& u, const AnFunc2& v): _u(u), _v(v) {}
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return _u(arrayOfPoints).dot(_v(arrayOfPoints));
    }
};
template <typename AnFunc1, typename AnFunc2> AnalyticFunctionCompDot<AnFunc1, AnFunc2>
    dot(const AnalyticFunction<AnFunc1>& u, const AnalyticFunction<AnFunc2>& v) {
   return AnalyticFunctionCompDot<AnFunc1, AnFunc2>(
        *static_cast<const AnFunc1*>(&u),*static_cast<const AnFunc2*>(&v));
}

template <typename AnFunc1> class AnalyticFunctionCompDotConst :
public AnalyticFunction <AnalyticFunctionCompDotConst<AnFunc1>> {
    const AnFunc1& _u;
    const Eigen::Matrix<Real, Eigen::Dynamic, 1> _v;
public:
    template <typename Derived> AnalyticFunctionCompDotConst(const AnFunc1& u, const Eigen::MatrixBase<Derived>& v): _u(u), _v(v) {}
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return _u(arrayOfPoints).dot(_v);
    }
};
template <typename AnFunc1, typename Derived> AnalyticFunctionCompDotConst<AnFunc1>
    dot(const AnalyticFunction<AnFunc1>& u, const Eigen::MatrixBase<Derived>& v) {
   return AnalyticFunctionCompDotConst<AnFunc1>(
        *static_cast<const AnFunc1*>(&u),v);
}
template <typename AnFunc1, typename Derived> AnalyticFunctionCompDotConst<AnFunc1>
    dot(const Eigen::MatrixBase<Derived>& v, const AnalyticFunction<AnFunc1>& u) {
   return AnalyticFunctionCompDotConst<AnFunc1>(
        *static_cast<const AnFunc1*>(&u),v);
}

template <typename AnFunc1, typename AnFunc2> class AnalyticFunctionCompCross :
public AnalyticFunction <AnalyticFunctionCompCross<AnFunc1,AnFunc2>> {
    const AnFunc1& _u;
    const AnFunc2& _v;
public:
    AnalyticFunctionCompCross(const AnFunc1& u, const AnFunc2& v): _u(u), _v(v) {}
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return (_u(arrayOfPoints).cross(_v(arrayOfPoints)));
    }
};
template <typename AnFunc1, typename AnFunc2> AnalyticFunctionCompCross<AnFunc1, AnFunc2>
    cross(const AnalyticFunction<AnFunc1>& u, const AnalyticFunction<AnFunc2>& v) {
   return AnalyticFunctionCompCross<AnFunc1, AnFunc2>(
        *static_cast<const AnFunc1*>(&u),*static_cast<const AnFunc2*>(&v));
}

template <typename AnFunc1> class AnalyticFunctionCompCrossConst :
public AnalyticFunction <AnalyticFunctionCompCrossConst<AnFunc1>> {
    const AnFunc1& _u;
    const Eigen::Matrix<Real, Eigen::Dynamic, 1> _v;
public:
    template <typename Derived> AnalyticFunctionCompCrossConst(const AnFunc1& u, const Eigen::MatrixBase<Derived>& v): _u(u), _v(v) {}
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& arrayOfPoints) const {
        return _u(arrayOfPoints).cross(_v);
    }
};
template <typename AnFunc1, typename Derived> AnalyticFunctionCompCrossConst<AnFunc1>
    cross(const AnalyticFunction<AnFunc1>& u, const Eigen::MatrixBase<Derived>& v) {
   return AnalyticFunctionCompCrossConst<AnFunc1>(
        *static_cast<const AnFunc1*>(&u),v);
}
template <typename AnFunc1, typename Derived> AnalyticFunctionCompCrossConst<AnFunc1>
    cross(const Eigen::MatrixBase<Derived>& v, const AnalyticFunction<AnFunc1>& u) {
    return AnalyticFunctionCompCrossConst<AnFunc1>(
        *static_cast<const AnFunc1*>(&u), -v);
}

extern AnalyticFunctionCoordinate x;
extern AnalyticFunctionCoordinate y;
extern AnalyticFunctionCoordinate z;

extern AnalyticFunctionPoint      k;
extern AnalyticFunctionCoordinate kx;
extern AnalyticFunctionCoordinate ky;
extern AnalyticFunctionCoordinate kz;

} // namespace AnalyticExpression

} // namespace Tortoise

#endif /* AnalyticFunction_hpp */
