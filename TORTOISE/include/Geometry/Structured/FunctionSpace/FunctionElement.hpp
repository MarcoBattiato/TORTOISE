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
//  FunctionElement.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//
// Function expressed in linear form (not nodal or modal)

#ifndef FunctionElement_hpp
#define FunctionElement_hpp

#include <Geometry/GeometryCore/Geometry.hpp>
#include <Generics/Features/MathFieldSpace.hpp>
#include <Geometry/Structured/DirectSpace/Mesh.hpp>
#include <Geometry/Structured/DirectSpace/MeshIterators.hpp>

#include <Eigen/Dense>

namespace Tortoise {

template <int NDim> class Mesh;

template<int NDim> class FunctionElement :                  // Inherits all the operations of a field from MathFieldSpace and AsymmetricMathFieldSpace
public Features::MathFieldSpace<FunctionElement<NDim>>,
public Features::AsymmetricMathFieldSpace<FunctionElement<NDim>, Real> { //,
    //public AsymmetricMathFieldSpace<FunctionElement<NDim>, std::function<Real(Point<NDim>)>>{
        
public:
    Mesh<NDim> const * const            mesh;
    int const                           elemID;
    GeometryCore::LinearForm<NDim>      vec;

public:
    //=======================================================
    // Constructors
    //===================
    // ATTENTION!!! The constructors are executed by linear form value, since the most frequent constructors will be from Function.
    // Nodal value constructors (the most intuitive for the final user) are not offered as the final user should not use this class
    // The only exception is the constant value constructor
    FunctionElement(const Mesh<NDim>& t_mesh, const MeshElementIterator<NDim>& elID);
    FunctionElement(const Mesh<NDim>& t_mesh, const MeshElementIterator<NDim>& elID, const Real t_value);
    template<typename Derived> FunctionElement(const Mesh<NDim>& t_mesh, const MeshElementIterator<NDim>& elID, const Eigen::MatrixBase<Derived>& functRepr): mesh(&t_mesh), elemID(elID.currentElementID), vec(functRepr) {};
    FunctionElement(const Mesh<NDim>& t_mesh, const int elID);
    FunctionElement(const Mesh<NDim>& t_mesh, const int elID, const Real t_value);
    template<typename Derived> FunctionElement(const Mesh<NDim>& t_mesh, const int elID, const Eigen::MatrixBase<Derived>& functRepr): mesh(&t_mesh), elemID(elID), vec(functRepr) {};
    FunctionElement(const Mesh<NDim>* t_mesh, const MeshElementIterator<NDim>& elID);
    FunctionElement(const Mesh<NDim>* t_mesh, const MeshElementIterator<NDim>& elID, const Real t_value);
    template<typename Derived> FunctionElement(const Mesh<NDim>* t_mesh, const MeshElementIterator<NDim>& elID, const Eigen::MatrixBase<Derived>& functRepr): mesh(t_mesh), elemID(elID.currentElementID), vec(functRepr) {};
    FunctionElement(const Mesh<NDim>* t_mesh, const int elID);
    FunctionElement(const Mesh<NDim>* t_mesh, const int elID, const Real t_value);
    template<typename Derived> FunctionElement(const Mesh<NDim>* t_mesh, const int elID, const Eigen::MatrixBase<Derived>& functRepr): mesh(t_mesh), elemID(elID), vec(functRepr) {};
    FunctionElement(const FunctionElement<NDim>& other);
    FunctionElement(FunctionElement<NDim>&& other);
        
    //********************************
    //* Arithmetic
    //********************************
    // ATTENTION!!!!!!! Multiplications and divisions by anything except a scalar do not conserve precision!
    // If you don't know what this means, do not use those operations! You will get unexpected results in most of the cases
    FunctionElement<NDim>& operator=(FunctionElement<NDim> other);
    FunctionElement<NDim>& operator=(FunctionElement<NDim>&& other);
    FunctionElement<NDim>& operator=(const Real& scalar);

    FunctionElement<NDim>& operator+=(const FunctionElement<NDim>& other);
    FunctionElement<NDim>& operator-=(const FunctionElement<NDim>& other);
    FunctionElement<NDim>& operator*=(const FunctionElement<NDim>& other);
    FunctionElement<NDim>& operator/=(const FunctionElement<NDim>& other);
    
    FunctionElement<NDim> operator-() const;
        
    FunctionElement<NDim>& operator+=(const Real& other);
    FunctionElement<NDim>& operator-=(const Real& other);
    FunctionElement<NDim>& operator*=(const Real& other);
    FunctionElement<NDim>& operator/=(const Real& other);
//    friend FunctionElement<NDim> operator/(const Real& lhs, FunctionElement<NDim> rhs);
    
    //********************************
    //* Local Evaluation
    //********************************
    Real operator()(const Point<NDim>& t_point_ref) const;                                      // Value at a given point in reference coordinates.
    template <typename Derived> auto operator() (const Eigen::MatrixBase<Derived>& points_ref) const;   // Value at several given points in reference coordinates.
        
    //********************************
    //* Other operations
    //********************************
    
    Real max() const;       // Max
    Real min() const;       // Min
    Real average() const;   // Average within an element
    Real integrate() const;
    Real integrate(const FunctionElement<NDim>& rhs) const;                                 // Calculates the integral of the function multiplied another one
    Real integrate(const FunctionElement<NDim>& rhs1, const FunctionElement<NDim>& rhs2) const;
    //    Function<NDim> Derivative();                                                              // Calculates the derivative
            
    //********************************
    //* I/O
    //********************************

        
public:
// Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
// See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    // Remove the possibility of passing temporaries to the constructor
    FunctionElement(const Mesh<NDim>&& t_mesh, const MeshElementIterator<NDim>& elID) = delete;
    FunctionElement(const Mesh<NDim>&& t_mesh, const MeshElementIterator<NDim>& elID, const Real t_value) = delete;
    template<typename Derived> FunctionElement(const Mesh<NDim>&& t_mesh, const MeshElementIterator<NDim>& elID, const Eigen::MatrixBase<Derived>& functRepr) = delete;
    FunctionElement(const Mesh<NDim>&& t_mesh, const int elID) = delete;
    FunctionElement(const Mesh<NDim>&& t_mesh, const int elID, const Real t_value) = delete;
    template<typename Derived> FunctionElement(const Mesh<NDim>&& t_mesh, const int elID, const Eigen::MatrixBase<Derived>& functRepr) = delete;
    
};

template<int NDim> std::ostream &operator<<(std::ostream &t_os, const FunctionElement<NDim>& funct);
template<int NDim> FunctionElement<NDim> operator/(const Real& lhs, FunctionElement<NDim> rhs);

} // namespace Tortoise
#endif /* FunctionElement_hpp */
