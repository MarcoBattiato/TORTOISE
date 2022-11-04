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
//  MeshIterators.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 23/6/21.
//

#ifndef MeshIterators_hpp
#define MeshIterators_hpp

#include <Geometry/GeometryCore/Geometry.hpp>

#include <cassert>
#include <stdio.h>

namespace Tortoise {

//*******************************
// Section Iterator
//*******************************
template <int NDim> class SectionIterator {
public:
    CartIndex<NDim>   maxcoord;
    CartIndex<NDim>   lcoord, hcoord;
    CartIndex<NDim>   currentcoord;
    bool              unfinished;
public:
    // Constructor from mesh
//    SectionIterator(const Mesh<NDim>& mesh);
//    SectionIterator(const Mesh<NDim>* mesh);
    
    // Explicit constructor
    SectionIterator(const CartIndex<NDim> t_maxcoord, const CartIndex<NDim> t_lcoord, const CartIndex<NDim> t_hcoord);
    
    // Copy and move constructors
    SectionIterator (const SectionIterator<NDim>& other);
    SectionIterator (SectionIterator<NDim>&& other);
    
    //=======================================================
    // Increment and assignment
    //===================
    // Prefix increment
    SectionIterator& operator++();
    // Assignment
    SectionIterator& operator=(const SectionIterator<NDim>& other);
    SectionIterator& operator=(SectionIterator<NDim>&& other);
    
public:
// Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
// See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
template<int NDim> std::ostream &operator<<(std::ostream &t_os, SectionIterator<NDim> const& iter);





//*******************************
// Element Iterators
//*******************************
// Element Iterators: meshes
//*******************************
// They do not follow the section's coordinate and simply iterate over an integer index
// This iterator is actually slower than simply iterating over the index, but it is introduced for consistency with the other iterators and
// especially to have overloaded iterators for mesh and for mesh subsets (which instead need more complicated iterators)
template <int NDim> class MeshElementIterator {
public:
    int     maxElemID;
    int     currentElementID;
    bool    unfinished;
public:
    // Constructor from mesh
//    MeshElementIterator(const Mesh<NDim>& mesh);
//    MeshElementIterator(const Mesh<NDim>* mesh);
    
    // Uninitialised iterator
    MeshElementIterator();
    
    // Explicit constructor
    MeshElementIterator(const int t_maxElemID, const int t_currentElementID);
    
    // Copy and move constructor
    MeshElementIterator(const MeshElementIterator<NDim>& other);
    MeshElementIterator(MeshElementIterator<NDim>&& other);
    
    //=======================================================
    // Increment and assignment
    //===================
    // Prefix increment
    MeshElementIterator& operator++();
    // Assignment
    MeshElementIterator& operator=(const MeshElementIterator<NDim>& other);
    MeshElementIterator& operator=(MeshElementIterator<NDim>&& other);
};
template<int NDim> std::ostream &operator<<(std::ostream &t_os, MeshElementIterator<NDim> const& iter);


//*******************************
// Element Iterators: mesh subsets
//*******************************
// More complex iterators that iterates only over the appropriate sections
template <int NDim> class MeshSubsetElementIterator {
public:
    CartIndex<NDim>     maxcoord;
    CartIndex<NDim>     lcoord, hcoord;
    CartIndex<NDim>     currentSeccoord;
    int                 currentElementinSectionID;
    int                 currentElementID;
    bool                unfinished;
public:
    // Constructor from mesh subsets
//    MeshSubsetElementIterator(const Mesh<NDim>& mesh);
//    MeshSubsetElementIterator(const Mesh<NDim>* mesh);
    
    // Explicit constructor
    MeshSubsetElementIterator(const CartIndex<NDim> t_maxcoord, const CartIndex<NDim> t_lcoord, const CartIndex<NDim> t_hcoord);
    
    // Copy and move constructors
    MeshSubsetElementIterator(const MeshSubsetElementIterator<NDim>& other);
    MeshSubsetElementIterator(MeshSubsetElementIterator<NDim>&& other);
    
    //=======================================================
    // Increment and assignment
    //===================
    // Prefix increment
    MeshSubsetElementIterator& operator++();
    // Assignment
    MeshSubsetElementIterator& operator=(const MeshSubsetElementIterator<NDim>& other);
    MeshSubsetElementIterator& operator=(MeshSubsetElementIterator<NDim>&& other);
    
    //
//    int currentElementID() const;
    
public:
// Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
// See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};
template<int NDim> std::ostream &operator<<(std::ostream &t_os, MeshSubsetElementIterator<NDim> const& iter);




//*******************************
// Vertex Iterator
//*******************************

template <int NDim> class VertexInElemIterator {
public:
    static constexpr int    maxcoord = NDim+1;
    int                     currentVertID;
    bool                    unfinished;
public:
    VertexInElemIterator ();
    VertexInElemIterator (const VertexInElemIterator<NDim>& other);
    VertexInElemIterator (VertexInElemIterator<NDim>&& other);
    VertexInElemIterator& operator++();
    VertexInElemIterator& operator=(const VertexInElemIterator& other);
    VertexInElemIterator& operator=(VertexInElemIterator&& other);
};
template<int NDim> std::ostream &operator<<(std::ostream &t_os, VertexInElemIterator<NDim> const& iter);

} // namespace Tortoise

#endif /* MeshIterators_hpp */
