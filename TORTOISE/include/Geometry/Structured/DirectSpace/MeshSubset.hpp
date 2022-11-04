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
//  MeshSubset.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 23/6/21.
//

#ifndef MeshSubset_hpp
#define MeshSubset_hpp

#include <Geometry/Structured/DirectSpace/Mesh.hpp>
#include <Geometry/Structured/DirectSpace/MeshIterators.hpp>
#include <Generics/Features/VectorSpace.hpp>

#include <cassert>

namespace Tortoise {

template <int NDim> class TrapezoidalSubset;

template <int NDim> class MeshSubset {
public:
    const Mesh<NDim>*       mesh;
    CartIndex<NDim>         lSecCoord;   // First included section coordinates
    CartIndex<NDim>         hSecCoord;   // First excluded section coordinates
    
public:
    //=======================================================
    // Constructors
    //===================
    // Subset covering the whole mesh
    MeshSubset(const Mesh<NDim>& t_mesh);
    MeshSubset(const Mesh<NDim>* t_mesh);
    
    // Subset with specified coverage
    MeshSubset(const Mesh<NDim>& t_mesh, const CartIndex<NDim>& t_lSecCoord, const CartIndex<NDim>& t_hSecCoord);
    MeshSubset(const Mesh<NDim>* t_mesh, const CartIndex<NDim>& t_lSecCoord, const CartIndex<NDim>& t_hSecCoord);
    
    // Subset covering the section pointed at by the SectionIterator
    MeshSubset(const Mesh<NDim>& t_mesh, const SectionIterator<NDim>& section);
    MeshSubset(const Mesh<NDim>& t_mesh, const CartIndex<NDim>& section);
    MeshSubset(const Mesh<NDim>* t_mesh, const SectionIterator<NDim>& section);
    MeshSubset(const Mesh<NDim>* t_mesh, const CartIndex<NDim>& section);

    // Subset from the TrapezoidalSubset
    MeshSubset(const Mesh<NDim>& t_mesh, const TrapezoidalSubset<NDim>& subset);
    MeshSubset(const Mesh<NDim>* t_mesh, const TrapezoidalSubset<NDim>& subset);
    
    //=======================================================
    // Converter to UnitSquareSubset
    //===================
    // Gives the representation of the subset in direct space
    TrapezoidalSubset<NDim> trapezoidalSubset() const;

    //=======================================================
    // Iterators
    //===================
    // Constructs an initialised iterator that allows for the iteration over all the sections included in the MeshSubset
    SectionIterator<NDim>           sectionIterator() const;
    MeshSubsetElementIterator<NDim> elementIterator() const;

    CartIndex<NDim>                 randomSection() const;
    MeshSubsetElementIterator<NDim> randomElement() const;
    
    int                             numberCoveredSections() const;
    
    //********************************
    //* Arithmetic
    //********************************
    // Allows for interval summation and accounts for umklapp
    // Notice that the operations do not make sense if the Region of definition of the MeshSubset are different
    MeshSubset<NDim>& operator = (const MeshSubset<NDim>& other);
    MeshSubset<NDim>& operator += (const MeshSubset<NDim>& other) { *this = trapezoidalSubset() + other.trapezoidalSubset(); return *this; }
    MeshSubset<NDim>& operator -= (const MeshSubset<NDim>& other) { *this = trapezoidalSubset() - other.trapezoidalSubset(); return *this; }
    MeshSubset<NDim>& operator += (const TrapezoidalSubset<NDim>& other) { *this = trapezoidalSubset() + other; return *this; }
    MeshSubset<NDim>& operator -= (const TrapezoidalSubset<NDim>& other) { *this = trapezoidalSubset() - other; return *this; }

    // !!!WARNING!!! This does not check if the operation is consistent: all the involved MeshSubset are defined on the same Region. It is up to the user to use this conscientiously
    MeshSubset& operator = (const TrapezoidalSubset<NDim>& subset);
    TrapezoidalSubset<NDim> operator-() const { return -trapezoidalSubset();} ;
    friend TrapezoidalSubset<NDim> operator+(const MeshSubset<NDim>& lhs, const MeshSubset<NDim>& rhs){ return lhs.trapezoidalSubset() + rhs.trapezoidalSubset();}
    friend TrapezoidalSubset<NDim> operator+(const MeshSubset<NDim>& lhs, const Point<NDim>& rhs){ return lhs.trapezoidalSubset() + rhs;}
    friend TrapezoidalSubset<NDim> operator+(const Point<NDim>& rhs, const MeshSubset<NDim>& lhs){ return lhs.trapezoidalSubset() + rhs;}
    friend TrapezoidalSubset<NDim> operator+(const MeshSubset<NDim>& lhs, TrapezoidalSubset<NDim> rhs){ return lhs.trapezoidalSubset() + rhs;}
    friend TrapezoidalSubset<NDim> operator+(TrapezoidalSubset<NDim> lhs, const MeshSubset<NDim>& rhs){ return lhs + rhs.trapezoidalSubset();}
    friend TrapezoidalSubset<NDim> operator-(const MeshSubset<NDim>& lhs, const MeshSubset<NDim>& rhs){ return lhs.trapezoidalSubset() - rhs.trapezoidalSubset();}
    friend TrapezoidalSubset<NDim> operator-(const MeshSubset<NDim>& lhs, TrapezoidalSubset<NDim> rhs){ return lhs.trapezoidalSubset() - rhs;}
    friend TrapezoidalSubset<NDim> operator-(const MeshSubset<NDim>& lhs, const Point<NDim>& rhs){ return lhs.trapezoidalSubset() - rhs;}
    friend TrapezoidalSubset<NDim> operator-(TrapezoidalSubset<NDim> lhs, const MeshSubset<NDim>& rhs){ return lhs - rhs.trapezoidalSubset();}
    friend TrapezoidalSubset<NDim> operator*(const Real scalar, const MeshSubset<NDim>& rhs){ return scalar * rhs.trapezoidalSubset();}
    friend TrapezoidalSubset<NDim> operator*(const MeshSubset<NDim>& lhs, const Real scalar){ return scalar * lhs.trapezoidalSubset();}
 
    
    
    //=======================================================
    // Technicalities
    //=======================================================
public:
// Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
// See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    // Remove the possibility of passing temporaries to the constructor
    MeshSubset(const Mesh<NDim>&& t_mesh) = delete;
    MeshSubset(const Mesh<NDim>&& t_mesh, const CartIndex<NDim>& t_lSecCoord, const CartIndex<NDim>& t_hSecCoord) = delete;
    MeshSubset(const Mesh<NDim>&& t_mesh, const SectionIterator<NDim>& section) = delete;
    MeshSubset(const Mesh<NDim>&& t_mesh, const CartIndex<NDim>& section) = delete;
    MeshSubset(const Mesh<NDim>&& t_mesh, const TrapezoidalSubset<NDim>& subset) = delete;
    
};

template<int NDim> std::ostream &operator<<(std::ostream &t_os, MeshSubset<NDim> const& meshSubset);
template<int NDim> std::ostream &operator<<(std::ostream &t_os, TrapezoidalSubset<NDim> const& meshSubset);
Plotter3D& operator << (Plotter3D& plotter, const MeshSubset<2>& subset);
Plotter3D& operator << (Plotter3D& plotter, const TrapezoidalSubset<2>& t_subsetToPlot);
// Example plotter3d << UnitSquareSubsetToPlot<2>(mesh, unitSquareSubset) << LABEL("UnitSquareSubset");






//=======================================================
//=======================================================
// UnitSquareSubset (Not to be used by user)
//===================
//===================
// Square subset within the unit square
// This will be used as all the subsets within the region will be treated in relative coordinates
template <int NDim> class TrapezoidalSubset : public VectorSpace<TrapezoidalSubset<NDim>,Real>, AsymmetricVectorSpace<TrapezoidalSubset<NDim>, Point<NDim>> {
public:
    Point<NDim>    lcoord, hcoord;
public:
    TrapezoidalSubset(const Point<NDim>& t_lcoord, const Point<NDim>& t_hcoord);
    TrapezoidalSubset operator-() const;
    TrapezoidalSubset& operator=(TrapezoidalSubset<NDim> other);
    TrapezoidalSubset& operator+=(const TrapezoidalSubset<NDim>& other);
    TrapezoidalSubset& operator-=(const TrapezoidalSubset<NDim>& other) ;
    TrapezoidalSubset& operator*=(const Real scalar) ;
    
    TrapezoidalSubset& operator+=(const Point<NDim>& other);
    TrapezoidalSubset& operator-=(const Point<NDim>& other);

    public:
    // Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
    // See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};




} // namespace Tortoise

#endif /* MeshSubset_hpp */
