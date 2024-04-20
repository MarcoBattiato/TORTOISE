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

#include <Generics/Utilities/UsefulMathFunctions.hpp>
#include <Geometry/GeometryCore/Point.hpp>
#include <Geometry/GeometryCore/ReferenceElement.hpp>
#include <Geometry/GeometryCore/GaussQuadratureConst.hpp>
#include <Geometry/GeometryCore/GaussQuadratureCartConst.hpp>
#include <Geometry/GeometryCore/LinearForm.hpp>



namespace Tortoise {

enum scattLegDirection { in = 1, out = -1};

// Flexible data type used to store the representation of entire functions, or data per element, or per section
// The dimension has to be set accordingly
using DataVector                            = Eigen::Matrix<Real, 1, Eigen::Dynamic>;


} // namespace Tortoise
#endif /* Geometry_hpp */
