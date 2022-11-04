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
//  Geometry.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 15/5/20.
//
// NOTE: I wrote a while ago that the code worked only for 2D and 3D. I do not remember anymore why I wrote it. I checked (not extensively) and it seems
// that it works for 1D. Please be careful and let me know if something breaks.
//
// This file is a header-only and provides functionalities to operate on a vector space and the space of its non-homogeneous linear forms.
// By non-homogeneous linear forms I mean a linear form plus a constant. For instance in 3D it would be an expression a*x+b*y+c*z+d.
//
// The vectors of the space are of type Point<NDim>, where NDim is the number of dimensions.
// As it will be common to treat many points and execute the same operation on all of them, two types are provided:  VectorPoint<NDim> and ArrayPoint<NDim>.
// ArrayPoint<NDim> is designed to hold a small number of points known at compile time.
// VectorPoint<NDim> is designed to hold a large number of points. The dimension related to the number of points is left
//  to be Dynamic since statically sized matrices in Eigen have a hard limit on their size.
//
// The file provides functions for the construction of these objects as well as for applying forms to points
//
//
//
//
// ================================================================================
// Points
// ================================================================================
// Points in NDim dimensional space
//  > TYPES
// Point<NDim>              : type storing a single NDim point
// VectorPoint<NDim>        : type storing a dynamic number of NDim points (designed to hold a large number of points)
// ArrayPoint<NDim,NPoints> : type storing a static number of NDim points (designed to hold a small number of points)
//
//  > FUNCTIONS
// VectorPoint<NDim> randomPoints(const int numberPoints)          : Constructs numberPoints random points in the reference cube [0,1]^NDim
// VectorPoint<NDim> randomPointsReference(const int numberPoints) : Constructs numberPoints random points in the reference tetrahedron
//
//
// ================================================================================
// Elements of space (Tetrahedra) and Reference tetrahedron
// ================================================================================
// The most fundamental polytope is the tetrahedron (called element when talking about meshes)
//  > TYPES
// Element<NDim>           : Stores the vertices of a given element (the i-th vertex is obtained with refNodes<NDim>.col(i) )
//
//  > CONSTANTS
// Element<NDim> refNodes  : Reference tetrahedron (the i-th vertex is obtained with refNodes<NDim>.col(i) )
//
//
// ================================================================================
// LinearForm<NDim>, VectorLinearForm<NDim>, ArrayLinearForm<NDim,NPoints>
// ================================================================================
// Non-homogeneous linear forms. For instance in 3D it would be an expression a*x+b*y+c*z+d.
//  > TYPES
// LinearForm<NDim>              : type storing a single NDim point
// VectorLinearForm<NDim>        : type storing a dynamic number of NDim points (designed to hold a large number of points)
// ArrayLinearForm<NDim,NPoints> : type storing a static number of NDim points (designed to hold a small number of points)
//


//
// ========================================
//  INTERNAL IMPLEMENTATION
// ========================================
//
// The internal implementation is powered by Eigen (http://eigen.tuxfamily.org)
// IMPORTANT: The use of templated functions to pass eigen arguments has been preferred, since Eigen::Ref does not work with block

#ifndef Geometry_hpp
#define Geometry_hpp

#include <Eigen/Dense>
#include <iostream>
#include <array>
#include <vector>
#include <math.h>
#include <functional>
//#include "DisableStupidWarnings.h"

namespace Tortoise {

enum scattLegDirection { in = 1, out = -1};
using Real = double;

constexpr inline int factorial(int t_n){
    if (t_n <= 0)
        return 1;
    return t_n * factorial(t_n - 1);
};

// The following function is taken from https://stackoverflow.com/questions/1505675/power-of-an-integer-in-c
// The author is Matthieu M. https://stackoverflow.com/users/147192/matthieu-m
template <int p> int constexpr IntPower(const int x) {
  if constexpr (p == 0) return 1;
  if constexpr (p == 1) return x;

  int tmp = IntPower<p / 2>(x);
  if constexpr ((p % 2) == 0) { return tmp * tmp; }
  else { return x * tmp * tmp; }
}


// Flexible data type used to store the representation of entire functions, or data per element, or per section
// The dimension has to be set accordingly
using DataVector                            = Eigen::Matrix<Real, 1, Eigen::Dynamic>;

//*******************************
// Point   P
//*******************************
// ===  Types
// ================
template <int NDim> using CartIndex                         = Eigen::Matrix<int, NDim, 1> ;
template <int NDim> using Point                             = Eigen::Matrix<Real, NDim, 1>;
// A point can be initialised easily     Point<2> exmplPoint = {1.2,-2.5};
template <int NDim> using VectorPoint                       = Eigen::Matrix<Real, NDim, Eigen::Dynamic>;
// VectorPoint is designed to hold a large number of points usually generated randomly by the functions defined below)
// In VectorPoint the dimension related to the number of points is left to be Dynamic since there is a hard limit to the size of statically stored matrices
template <int NDim, int NPoints> using ArrayPoint           = Eigen::Matrix<Real, NDim, NPoints>;
// A vector of points with statically assigned dimension (to be used for small number of points)

// Specialised types with special meaning
template <int NDim> using ElemNodalData             = Eigen::Matrix<Real, 1, NDim+1>;
// Values at nodes of a tetrahedron
template <int NDim> using Element                   = ArrayPoint<NDim, NDim+1>;
// Stores the vertices of a given element


// ===  Constants
// ================
template <int NDim> Element<NDim> refNodes = [] {
    ArrayPoint<NDim, NDim+1> tmp = ArrayPoint<NDim, NDim+1>::Zero();
    for (int i=0; i<NDim; ++i){ tmp(i,i+1) = 1.0; }
    return tmp; }();

// ===  Creators
// ================
//template <int NDim> Point<NDim>       origin()


template <int NDim> VectorPoint<NDim> randomPoints(const int numberPoints){
    return 0.5*Eigen::Array<Real, NDim, Eigen::Dynamic>::Random(NDim,numberPoints)+0.5;
}                // Constructs random coordinates in the range [0,1]
template <int NDim> VectorPoint<NDim> randomPointsReference(const int numberPoints){
    VectorPoint<NDim> toreturn(NDim,numberPoints);
    toreturn = 0.5*Eigen::Array<Real, NDim, Eigen::Dynamic>::Random(NDim,numberPoints)+0.5;
    switch (NDim) {
        case 2: {
    // cut'n fold the square into the reference triangle
            for (int i= 0; i<numberPoints; ++i){
                if (toreturn(0,i)+toreturn(1,i)>1.0) {toreturn(0,i) = 1.0 -toreturn(0,i); toreturn(1,i) = 1.0 -toreturn(1,i);}
            }
            break;
        }
        case 3: {
    // cut'n fold the square into the reference tetrahedron   see http://vcg.isti.cnr.it/jgt/tetra.htm
            for (int i= 0; i<numberPoints; ++i){
                if (toreturn(0,i)+toreturn(1,i)>1.0) {toreturn(0,i) = 1.0 -toreturn(0,i);toreturn(1,i) = 1.0 -toreturn(1,i);}                                       // cut'n fold the cube into a prism
                if(toreturn(1,i)+toreturn(2,i)>1.0) {Real tmp = toreturn(2,i); toreturn(2,i) = 1.0 - toreturn(0,i) - toreturn(1,i); toreturn(1,i) = 1.0 - tmp;    // cut'n fold the prism into a tetrahedron
                } else if(toreturn(0,i)+toreturn(1,i)+toreturn(2,i)>1.0) {Real tmp = toreturn(2,i); toreturn(2,i) = toreturn(0,i) + toreturn(1,i) + toreturn(2,i) - 1.0; toreturn(0,i) = 1 - toreturn(1,i) - tmp;}
            }
            break;
        }
        default:
            break;
    }
    return  toreturn;
}       // Constructs random points in the reference tetrahedron

// *******
// TO BE DEVELOPED: Gauss points by splitting tetrahedron in smaller tetrahedra

//template <int NDim, int order> VectorPoint<NDim> splitGaussPointsReference(){
//
//}



//*******************************
// Non-homogeneous linear forms  ( Functionals of Points   f(P) )
//*******************************
// ===  Types
// ================
template <int NDim> using LinearForm                        = Eigen::Matrix<Real, 1, NDim+1>;
// It is defines such that
// LinearForm<2> exmplLinForm = {-0.3,2.3,1.5};
// applying exmplLinForm to exmplPoint returns   exmplLinForm(0) + exmplLinForm(1)+exmplPoint(0) + exmplLinForm(2)+exmplPoint(1)
template <int NDim> using VectorLinearForm                  = Eigen::Matrix<Real, Eigen::Dynamic, NDim+1>;
// In VectorLinearForm the dimension related to the number of linear forms is left to be Dynamic since there is a hard limit to the size of statically stored Matrices
template <int NDim, int NForms> using ArrayLinearForm       = Eigen::Matrix<Real, NForms, NDim+1>;
// A vector of linear forms with statically assigned dimension (to be used for small number of linear forms)

// ===  Creators
// ================
template <int NDim> LinearForm<NDim> createLinearForm(const std::array<Real,NDim+1> nodevalue){
    LinearForm<NDim> toreturn;
    toreturn(0) = nodevalue[0];
    for (int i=1; i<NDim+1; ++i){ toreturn(i) = nodevalue[i]-nodevalue[0];}
    return toreturn;
};
// Creates a linear form that takes the given values at the reference tetrahedron's node coordinate. The nodes are ordered {0,0,0}, {1,0,0}, {0,1,0}. {0,0,1}
// nodevalue is in the form Eigen::Matrix<Real, nLinForms, NDim+1>
template <typename Derived> auto createLinearForm(const Eigen::MatrixBase<Derived>& nodevalue){
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> toreturn = nodevalue;
    toreturn.rightCols(toreturn.cols()-1).colwise() -= toreturn.col(0);
    return toreturn;
};
// Creates a linear form that takes the given values at the reference tetrahedron's node coordinate. The nodes are ordered {0,0,0}, {1,0,0}, {0,1,0}. {0,0,1}
template <int NDim> LinearForm<NDim> createLinearForm(const std::function<Real(Point<NDim>)>& t_f){
    LinearForm<NDim> values;
    for (int i = 0; i<NDim+1; ++i){ values(i) = t_f(refNodes<NDim>.col(i));}
    return createLinearForm(values);
}
// Creates a linear form that approximates the given function. The linear form has the same values of the function at the reference tetrahedron's node coordinate.

// ===  Operations
// ================
// This function conveniently works for all the following combinations
//  ->  apply(LinForm, Point)           ->  apply(LinForm, PointVector)             ->  apply(LinForm, ArrayPoint)
//  ->  apply(VectorLinForm, Point)     ->  apply(VectorLinForm, VectorPoint)       ->  apply(VectorLinForm, ArrayPoint)
//  ->  apply(ArrayLinearForm, Point)   ->  apply(ArrayLinearForm, VectorPoint)     ->  apply(ArrayLinearForm, ArrayPoint)
template <typename DerivedA,typename DerivedB> auto apply(const Eigen::MatrixBase<DerivedA>& linearForm, const Eigen::MatrixBase<DerivedB>& point){
    return (linearForm.rightCols(linearForm.cols()-1)*point).colwise()+linearForm.col(0);
}

//*******************************
// Linear operators on the non-homogeneous linear forms
//*******************************
template <int NDim> using LinTransform   = Eigen::Matrix<Real, NDim+1, NDim+1>;





//*******************************
// Gauss quadrature
//*******************************
// The gauss points and weights are for reference elements
// The reference element in 1D is the [0,1] interval
// The reference element in 2D is the {[0,0],[1,0],[0,1]} triangle
// The reference element in 3D is the {[0,0,0],[1,0,0],[0,1,0],[0,0,1]} tetrahedron
constexpr int ngausspoint(const int dim, const int ord) {
    switch (dim) {
        case 1: return ord; break;
        case 2: switch (ord) {
                    case 1: return 1; break;
                    case 2: return 3; break;
                    case 3: return 4; break;
                    case 4: return 6; break;
                    case 5: return 7; break;
                    case 6: return 12; break;
                    case 7: return 13; break;
                    default: return -1; break;
                }
        case 3: switch (ord) {
                    case 1: return 1; break;
                    case 2: return 4; break;
                    case 3: return 5; break;
                    case 4: return 10; break;
                    case 5: return 11; break;
                    case 6: return 16; break;
                    case 7: return 13; break;
                    default: return -1; break;
            }
        default: return -1; break;
    }
}

template <int NDim, int order> const ArrayPoint<NDim,ngausspoint(NDim,order)> gaussPoints;
template <int NDim, int order> const ArrayPoint<1,ngausspoint(NDim,order)> gaussWeigths;

// 1D
// ==============
template<> const inline ArrayPoint<1,1> gaussPoints<1,1> {0.5};
template<> const inline ArrayPoint<1,1> gaussWeigths<1,1> {1.0};
// -----
template<> const inline ArrayPoint<1,2> gaussPoints<1,2> {0.2113248654051871,0.7886751345948129};
template<> const inline ArrayPoint<1,2> gaussWeigths<1,2> {0.5,0.5};
// -----
template<> const inline ArrayPoint<1,3> gaussPoints<1,3> {0.1127016653792583,0.5,0.8872983346207417};
template<> const inline ArrayPoint<1,3> gaussWeigths<1,3> {0.2777777777777778,0.4444444444444444,0.2777777777777778};
// -----
template<> const inline ArrayPoint<1,4> gaussPoints<1,4> {0.0694318442029737,0.3300094782075719,0.6699905217924281,0.9305681557970263};
template<> const inline ArrayPoint<1,4> gaussWeigths<1,4> {0.1739274225687269,0.3260725774312731,0.3260725774312731,0.1739274225687269};
// -----
template<> const inline ArrayPoint<1,5> gaussPoints<1,5> = [] {ArrayPoint<1,5> tmp;
    tmp << 0.0469100770306680,0.0230765344947158,0.5,0.7692346550528415,0.9530899229693319;
    return tmp; }();
template<> const inline ArrayPoint<1,5> gaussWeigths<1,5> = [] {ArrayPoint<1,5> tmp;
    tmp << 0.1184634425280945,0.2393143352496832,0.2844444444444444,0.2393143352496832,0.1184634425280945;
    return tmp;}();
// 2D
// ==============
template<> const inline ArrayPoint<2,ngausspoint(2,1)> gaussPoints<2,1> = [] {ArrayPoint<2,ngausspoint(2,1)> tmp;
tmp << 0.33333333333333, 0.33333333333333;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(2,1)> gaussWeigths<2,1> = [] {ArrayPoint<1,ngausspoint(2,1)> tmp;
tmp << 1.00000000000000;
return tmp; }();
// -----
template<> const inline ArrayPoint<2,ngausspoint(2,2)> gaussPoints<2,2> = [] {ArrayPoint<2,ngausspoint(2,2)> tmp;
tmp << 0.16666666666667, 0.66666666666667, 0.16666666666667,
    0.16666666666667, 0.16666666666667, 0.66666666666667;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(2,2)> gaussWeigths<2,2> = [] {ArrayPoint<1,ngausspoint(2,2)> tmp;
tmp << 0.33333333333333, 0.33333333333333, 0.33333333333333;
return tmp; }();
// -----
template<> const inline ArrayPoint<2,ngausspoint(2,3)> gaussPoints<2,3> = [] {ArrayPoint<2,ngausspoint(2,3)> tmp;
tmp << 0.20000000000000, 0.60000000000000, 0.20000000000000, 0.33333333333333,
    0.20000000000000, 0.20000000000000, 0.60000000000000, 0.33333333333333;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(2,3)> gaussWeigths<2,3> = [] {ArrayPoint<1,ngausspoint(2,3)> tmp;
tmp << 0.52083333333333, 0.52083333333333, 0.52083333333333, -0.5625000000000;
return tmp; }();
// -----
template<> const inline ArrayPoint<2,ngausspoint(2,4)> gaussPoints<2,4> = [] {ArrayPoint<2,ngausspoint(2,4)> tmp;
tmp << 0.09157621350977, 0.44594849091597, 0.81684757298046,  0.44594849091597, 0.09157621350977, 0.10810301816807,
    0.09157621350977, 0.10810301816807, 0.09157621350977, 0.44594849091597, 0.81684757298046, 0.44594849091597;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(2,4)> gaussWeigths<2,4> = [] {ArrayPoint<1,ngausspoint(2,4)> tmp;
    tmp << 0.10995174365532, 0.22338158967801, 0.10995174365532, 0.22338158967801, 0.10995174365532, 0.22338158967801;
return tmp; }();
// -----
template<> const inline ArrayPoint<2,ngausspoint(2,5)> gaussPoints<2,5> = [] {ArrayPoint<2,ngausspoint(2,5)> tmp;
tmp << 0.10128650732346, 0.47014206410511, 0.79742698535309, 0.47014206410511, 0.10128650732346, 0.05971587178977, 0.33333333333333,
    0.10128650732346, 0.05971587178977, 0.10128650732346, 0.47014206410511, 0.79742698535309, 0.47014206410511, 0.33333333333333;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(2,5)> gaussWeigths<2,5> = [] {ArrayPoint<1,ngausspoint(2,5)> tmp;
    tmp << 0.12593918054483, 0.13239415278851, 0.12593918054483, 0.13239415278851, 0.12593918054483, 0.13239415278851, 0.22500000000000;
return tmp; }();
// 3D
// ==============
template<> const inline ArrayPoint<3,ngausspoint(3,1)> gaussPoints<3,1> = [] {ArrayPoint<3,ngausspoint(3,1)> tmp;
tmp << 0.250000000000,
    0.2500000000000,
    0.2500000000000;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(3,1)> gaussWeigths<3,1> = [] {ArrayPoint<1,ngausspoint(3,1)> tmp;
tmp << 1.00000000000000;
return tmp; }();
// -----
template<> const inline ArrayPoint<3,ngausspoint(3,2)> gaussPoints<3,2> = [] {ArrayPoint<3,ngausspoint(3,2)> tmp;
tmp << 0.585410196624, 0.138196601125, 0.138196601125, 0.138196601125,
    0.1381966011250, 0.5854101966249, 0.1381966011250, 0.1381966011250,
    0.1381966011250, 0.1381966011250, 0.5854101966249, 0.1381966011250;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(3,2)> gaussWeigths<3,2> = [] {ArrayPoint<1,ngausspoint(3,2)> tmp;
    tmp << 0.25000000000000, 0.25000000000000, 0.25000000000000, 0.25000000000000;
return tmp; }();
// -----
template<> const inline ArrayPoint<3,ngausspoint(3,3)> gaussPoints<3,3> = [] {ArrayPoint<3,ngausspoint(3,3)> tmp;
    tmp << 0.250000000000, 0.500000000000, 0.166666666666, 0.166666666666, 0.166666666666,
    0.2500000000000, 0.1666666666666, 0.1666666666666, 0.1666666666666, 0.5000000000000,
    0.2500000000000, 0.1666666666666, 0.1666666666666, 0.5000000000000, 0.1666666666666;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(3,3)> gaussWeigths<3,3> = [] {ArrayPoint<1,ngausspoint(3,3)> tmp;
    tmp << -0.8000000000000, 0.45000000000000, 0.45000000000000, 0.45000000000000, 0.45000000000000;
return tmp; }();
// -----
template<> const inline ArrayPoint<3,ngausspoint(3,4)> gaussPoints<3,4> = [] {ArrayPoint<3,ngausspoint(3,4)> tmp;
    tmp << 0.568430584196, 0.143856471934, 0.143856471934, 0.143856471934, 0.000000000000, 0.500000000000, 0.500000000000, 0.500000000000, 0.000000000000, 0.000000000000,
    0.1438564719343, 0.1438564719343, 0.1438564719343, 0.5684305841968, 0.5000000000000, 0.0000000000000, 0.5000000000000, 0.0000000000000, 0.5000000000000, 0.0000000000000,
    0.1438564719343, 0.1438564719343, 0.5684305841968, 0.1438564719343, 0.5000000000000, 0.5000000000000, 0.0000000000000, 0.0000000000000, 0.0000000000000, 0.5000000000000;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(3,4)> gaussWeigths<3,4> = [] {ArrayPoint<1,ngausspoint(3,4)> tmp;
    tmp << 0.21776506988041, 0.21776506988041, 0.21776506988041, 0.21776506988041, 0.02148995341306, 0.02148995341306, 0.02148995341306, 0.02148995341306, 0.02148995341306, 0.02148995341306;
return tmp; }();
// -----
template<> const inline ArrayPoint<3,ngausspoint(3,5)> gaussPoints<3,5> = [] {ArrayPoint<3,ngausspoint(3,5)> tmp;
    tmp << 0.250000000000, 0.785714285714, 0.071428571428, 0.071428571428, 0.071428571428, 0.100596423833, 0.399403576166, 0.399403576166, 0.399403576166, 0.100596423833, 0.100596423833,
    0.2500000000000, 0.0714285714285, 0.0714285714285, 0.0714285714285, 0.7857142857142, 0.3994035761668, 0.1005964238332, 0.3994035761668, 0.1005964238332, 0.3994035761668, 0.1005964238330,
    0.2500000000000, 0.0714285714285, 0.0714285714285, 0.7857142857142, 0.0714285714285, 0.3994035761668, 0.3994035761668, 0.1005964238332, 0.1005964238332, 0.1005964238332, 0.3994035761668;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(3,5)> gaussWeigths<3,5> = [] {ArrayPoint<1,ngausspoint(3,5)> tmp;
    tmp << -0.0789333333333, 0.04573333333333, 0.04573333333333, 0.04573333333333, 0.04573333333333, 0.14933333333333, 0.14933333333333, 0.14933333333333, 0.14933333333333, 0.14933333333333, 0.14933333333333;
return tmp; }();

} // namespace Tortoise
#endif /* Geometry_hpp */
