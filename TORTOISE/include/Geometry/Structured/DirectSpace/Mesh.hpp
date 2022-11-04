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
//  Mesh.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//
// Class that allows for the construction of meshes.
// A Mesh has to be created to be within a Region. The mesh can cover any subset of a Region, but its sides must be aligned to the Region sides.
// The reason for this limitation is to ensure an important optimisation for the scattering calculations, but it is too technical to be explained here.
//


#ifndef Mesh_hpp
#define Mesh_hpp

#include <Geometry/Structured/DirectSpace/Region.hpp>
#include <Geometry/Structured/DirectSpace/MeshIterators.hpp>
#include <Geometry/Plotter/Plotter.hpp>

#include <array>
#include <initializer_list>

namespace Tortoise {

//*******************************
// Types definitions
//*******************************

// Data type used to store the nodal representation of functions in each element.
// WARNING: ElemNodalData and LinearForm are not distinguished by the compiler but they are used to store different
// representations. They should not be used interchangeably (as they are intepreted differently), but it is up to the
// user not to use them improperly. The LinearForm representation is contructed by    createLinearForm(ElemNodalData)

template <int NDim> using VectorElement     = std::vector<Eigen::Matrix<Real, NDim, NDim+1>, Eigen::aligned_allocator<Eigen::Matrix<Real, NDim, NDim+1> > > ;
//
using VecOfDataVector                       = std::vector<Eigen::Matrix<Real, 1, Eigen::Dynamic>>;

template<int NDim> class ElementID;

//*******************************
// Main class definition
//*******************************

template<int NDim> class Mesh {         // NDim = Number of Dimensions of the space
    // Geometric Attributes
public:
    //*******************************
    // Dimensionality constants
    //*******************************
    static constexpr int                            nVertElement = NDim+1;                  // Number of vertices in Element
    static constexpr int                            nElementSection = factorial(NDim);      // Number of element per reference Section
    static const VectorElement<NDim>                vertElemInRefSec;                       // Vertices of Ref Elements (They all share the first vertex at the origin)

public: // DO NOT CHANCE THE ORDER!!! IT WILL BREAK THE CONSTRUCTOR!!!
    //*******************************
    // Mesh features
    //*******************************
    const Region<NDim>*                             region;                                 // Region the mesh is part of

    const Point<NDim>                               origin;                                 // origin of the Mesh
    const ArrayPoint<NDim,NDim>                     gVec;                                   // Mesh sides vectors
    const ArrayPoint<NDim,NDim>                     invgVec;                                // Inverse of Matrix formed by Mesh sides vectors

    const CartIndex<NDim>                           nSecSplits;                             // Number of sectioning per dimensions

    const int                                       numberSections;                         // Total number of sections
    const int                                       numberElements;                         // Total number of elements
    const int                                       dataVectorDim;                          // Total number of basis functions

    const ArrayPoint<NDim,NDim>                     dGVec;                                  // Sides of the sections
    const ArrayPoint<NDim,NDim>                     invdGVec;                               // Inverse of Matrix formed by sides of the sections

    const Real                                      elemVolume;                             // Area of the element
    const LinTransform<NDim>                        massMatrix;                             // Mass matrix for linear functional representation
    const LinTransform<NDim>                        invMassMatrix;                          // Inverse of Mass matrix for linear functional representation

    const Point<NDim>                               relativeOrigin;                         // origin of the Mesh in fraction of Region
    const Point<NDim>                               relativegVec;                           // Mesh sides vectors in fraction of Region
    const Point<NDim>                               relativedGVec;                          // Mesh sides vectors in fraction of Region

    const VecOfDataVector                           verticesCoord;                          // Stores the coordinates of the vertices

    const VectorLinearForm<NDim>                    basisFuncLF;                            // Linear forms corresponding to the basis functions
    const std::vector<VectorLinearForm<NDim>>       basisFuncDerivLF;                       // Linear forms corresponding to the derivatives of the basis functions
    const VecOfDataVector                           basisFuncDerivLF0thOrd;                 // Linear forms corresponding to the derivatives of the basis functions
    
    // Methods
public:
    //=======================================================
    // Constructors
    // !!! Notice !!! gVec are the lattice vectors and NOT the end points of the mesh.
    // The endpoints of the mesh are origin+gVec[1], origin+gVec[2], ...
    // The origin has to be passed in relative position compared to the region's origin and gVec's
    //===================
    Mesh(const Region<NDim>& t_region, const Point<NDim>& t_relativeOrigin, const Point<NDim>& t_relativegVec, const CartIndex<NDim>& t_nSecSplits);
    Mesh(const Region<NDim>* t_region, const Point<NDim>& t_relativeOrigin, const Point<NDim>& t_relativegVec, const CartIndex<NDim>& t_nSecSplits);
    // In 1D a constructor accepting Real numbers is provided
    Mesh(const Region<NDim>& t_region, const Real t_relativeOrigin, const Real t_relativegVec, const int t_nSecSplits) requires (NDim == 1);
    Mesh(const Region<NDim>* t_region, const Real t_relativeOrigin, const Real t_relativegVec, const int t_nSecSplits) requires (NDim == 1);
    //===================
    
    //********************************
    //* Mesh structure
    //********************************
    //=======================================================
    // Mesh Iterators and accessor
    // Iterating through the mesh should be done as
    // >>  for (auto elem = meshexmpl.elementIterator(); elem.unfinished; ++elem) {
    // >>       something[elem]
    // >>  }
    // This allows for abstraction and also works for MeshSubsets that run only over some of the elements
    // Notice that the iterator is always initialised to the first section, element, or vertex
    //===================
    // Indices conversions (more technical, the final user will not need these)
    CartIndex<NDim> sectionID(const int t_elementID) const;
    int             sectionID(const CartIndex<NDim>& SectionCartID) const;
    int             elemInSectionID(const int t_elementID) const;
    int             elemID(const CartIndex<NDim>& SectionCartID, const int elemInSect) const;
    int             elemID(const std::array<int, NDim>& SectionCartID, const int elemInSect) const;

    // Mesh iterators
    SectionIterator<NDim>       sectionIterator() const;
    MeshElementIterator<NDim>   elementIterator() const;
    VertexInElemIterator<NDim>  vertexInElemIterator() const;
    
    // Random accessors
    CartIndex<NDim>             randomSection() const;
    SectionIterator<NDim>       randomSectionIterator() const;
    int                         randomElement() const;
    
    
    //=======================================================
    // Mesh geometry
    //===================
    // Mesh Points from section ID
    Point<NDim> operator()(const std::array<int, NDim>& SectionCartID) const;                                            // Mesh point in cartesian coordinates
    Point<NDim> operator()(const SectionIterator<NDim>& t_sectionID) const;                                              // Mesh point from section ID
    Point<NDim> operator()(const CartIndex<NDim>& t_sectionID) const;                                                    // Mesh point from section ID
    // Vertex from element ID and vertex in element ID
    Point<NDim> operator()(const MeshElementIterator<NDim>& t_elementID, const VertexInElemIterator<NDim>& t_nVertex) const;
    Point<NDim> operator()(const int t_elementID, const int t_nVertex) const;
    Point<NDim> operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const VertexInElemIterator<NDim>& t_nVertex) const;
    // Point from element ID and coordinates in reference element
    Point<NDim> operator()(const MeshElementIterator<NDim>& t_elementID, const Point<NDim>& refCoord) const;
    Point<NDim> operator()(const MeshElementIterator<NDim>& t_elementID, const Real refCoord) const requires (NDim == 1);
    Point<NDim> operator()(const int t_elementID, const Point<NDim>& refCoord) const;
    Point<NDim> operator()(const int t_elementID, const Real refCoord) const requires (NDim == 1);
    Point<NDim> operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const Point<NDim>& refCoord) const;
    Point<NDim> operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const Real refCoord) const requires (NDim == 1);
    // Baricentre of a given element
    Point<NDim> elemCentre(const MeshElementIterator<NDim>& t_elementID) const;
    Point<NDim> elemCentre(const MeshSubsetElementIterator<NDim>& t_elementID) const;
    Point<NDim> elemCentre(const int t_elementID) const;
    
    
    //=======================================================
    // Adjacency Map and Faces
    //===================
    int             adjacentElement(const int elementID, const int face) const;
    int             adjacentFace(const int elementID, const int face) const;
    Point<NDim>     faceNormal(const int elementID, const int face) const;
    Real            faceArea(const int elementID, const int face) const;
    
    //********************************
    //* Associated DG space
    //********************************
    //=======================================================
    // DataVectors
    //===================
    DataVector representationVector() const { return DataVector::Zero(dataVectorDim);}
    DataVector elemDataVector() const { return DataVector::Zero(numberElements);}
    DataVector sectDataVector() const { return DataVector::Zero(numberSections);}
    
    //===================
    
    //=======================================================
    // I/O
    //===================
    void plot(const std::string& t_title = "") const;
    void plotDetailed(const std::string& t_title = "") const;
    void writeToTxtFile(const std::string &FileName) const;
    
    //=======================================================
    // Technicalities
    //=======================================================
    public:
    // Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
    // See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // Removes the possibility of passing temporaries to the constructor
    Mesh(const Region<NDim>&& t_region, const Point<NDim>& t_relativeOrigin, const Point<NDim>& t_relativegVec, const CartIndex<NDim>& t_nSecSplits) = delete;
    Mesh(const Region<NDim>&& t_region, const Real t_relativeOrigin, const Real t_relativegVec, const int t_nSecSplits) requires (NDim ==1) = delete;


};

template<int NDim> std::ostream &operator<<(std::ostream &t_os, Mesh<NDim> const& t_mesh);
Plotter3D& operator << (Plotter3D& plotter, const Mesh<3>& t_mesh);
Plotter3D& operator << (Plotter3D& plotter, const Mesh<2>& t_mesh);
Plotter2D& operator << (Plotter2D& plotter, const Mesh<2>& t_mesh);
Plotter2D& operator << (Plotter2D& plotter, const Mesh<1>& t_mesh);

template<int NDim> bool operator==(const Mesh<NDim>& lhs, const Mesh<NDim>& rhs);

} // namespace Tortoise

#endif /* Mesh_hpp */
