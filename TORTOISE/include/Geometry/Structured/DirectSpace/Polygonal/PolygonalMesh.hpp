//
//  PolygonalMesh.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 16/9/22.
//

#ifndef PolygonalMesh_hpp
#define PolygonalMesh_hpp

#include <Geometry/GeometryCore/Geometry.hpp>
#include <Geometry/Plotter/Plotter.hpp>
#include <Geometry/Structured/DirectSpace/Mesh.hpp>
#include <Geometry/Structured/FunctionSpace/Function.hpp>
#include <map>
#include <array>

namespace Tortoise {

template <int NDim> class PolygonalSectionIterator;
template <int NDim> class PolygonalElementIterator;

template <int NDim> using VectorCartCoord          = Eigen::Matrix<int, NDim, Eigen::Dynamic>;

// The structure of a PolygonalMesh, for instance in 2D,
//  x \ y | 0  1  2  3  4  5  6  7
// ----------------------------------
//    0   | 0  1  1  0  1  1  1  0
//    1   | 0  1  1  0  0  1  1  0
//    2   | 0  0  0  0  0  0  0  0
//    3   | 1  1  1  1  1  0  0  0
//    4   | 1  1  1  0  1  1  1  1
//    5   | 0  0  0  0  0  0  0  0
//
// is saved as a number of horizontal stripes of stored length
// 0 :
//      1  : (2, 0)
//      4  : (3, 2)
// 1 :
//      1  : (2, 5)
//      5  : (2, 7)
// 3 :
//      0  : (5, 9)
// 4 :
//      0  : (3, 14)
//      4  : (4, 17)
//
// The first index is the index of the section row
// The second index is the  index of the section column where the stripe begins
// The first stored integer is the length of the stripe
// The second stored index is the position of the cumulative number of sections until then
// (it will be used to calculate the position of data in vectors storing information on the
// PolygonalMesh, for instance the vector stored in a Function)

template <int NDim> class PolygonalMeshStructure;
template <> class PolygonalMeshStructure<1>: public std::map<int,std::array<int,2>> {};
template <> class PolygonalMeshStructure<2>: public std::map<int,std::map<int,std::array<int,2>>> {};
template <> class PolygonalMeshStructure<3>: public std::map<int,std::map<int,std::map<int,std::array<int,2>>>>{};


template<int NDim> class PolygonalMesh {         // NDim = Number of Dimensions of the space
// A PolygonalMesh is constructed as a series of stripes
    
    // Geometric Attributes
public:
    //*******************************
    // Dimensionality constants
    //*******************************
    static constexpr int                            nVertElement = NDim+1;                  // Number of vertices in Element
    static constexpr int                            nElementSection = Utilities::factorial(NDim);      // Number of element per reference Section
    static const VectorElement<NDim>                vertElemInRefSec ;  // Vertices of Ref Elements (They all share the first vertex at the origin)

public: // DO NOT CHANCE THE ORDER!!! IT WILL BREAK THE CONSTRUCTOR!!!
    //*******************************
    // Mesh features
    //*******************************
    const Region<NDim>*                             region;                                 // Region the polygonal mesh is part of
    const Mesh<NDim>*                               meshFull;                                   // Mesh the polygonal mesh is derived from
    const PolygonalMeshStructure<NDim>              meshStructure;                          // Contains the shape of the PolygonalMesh
    
    const int                                       numberSections;                         // Total number of sections
    const int                                       numberElements;                         // Total number of elements
    const int                                       dataVectorDim;                          // Total number of basis functions

    const Point<NDim>                               origin;                                 // origin of the Mesh
    const ArrayPoint<NDim,NDim>                     dGVec;                                  // Sides of the sections
    const ArrayPoint<NDim,NDim>                     invdGVec;                               // Inverse of Matrix formed by sides of the sections

    const Real                                      elemVolume;                             // Area of the element
    const GeometryCore::LinTransform<NDim>          massMatrix;                             // Mass matrix for linear functional representation
    const GeometryCore::LinTransform<NDim>          invMassMatrix;                          // Inverse of Mass matrix for linear functional representation
    

    // Methods
public:
    PolygonalMesh(const Mesh<NDim>& t_meshFull, const Eigen::Matrix<bool,1,Eigen::Dynamic>& t_whichSections);
    PolygonalMesh(const Mesh<NDim>* t_meshFull, const Eigen::Matrix<bool,1,Eigen::Dynamic>& t_whichSections);
    PolygonalMesh(const Function<NDim> &t_function, const std::function<bool(Real)>& boolFunctionToEvaluate, const RequirementType requirementType);
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
    PolygonalSectionIterator<NDim>      sectionIterator() const;
    PolygonalElementIterator<NDim>      elementIterator() const;
    VertexInElemIterator<NDim>          vertexInElemIterator() const;

    //=======================================================
    // Mesh geometry
    //===================
    // Mesh Points from section ID
    Point<NDim> operator()(const PolygonalSectionIterator<NDim>& t_sectionID) const;                    // Mesh point from section ID
    // Vertex from element ID and vertex in element ID
    Point<NDim> operator()(const PolygonalElementIterator<NDim>& t_elementID, const VertexInElemIterator<NDim>& t_nVertex) const;
    // Point from element ID and coordinates in reference element
    Point<NDim> operator()(const PolygonalElementIterator<NDim>& t_elementID, const Point<NDim>& refCoord) const;
    Point<NDim> operator()(const PolygonalElementIterator<NDim>& t_elementID, const Real refCoord) const requires (NDim == 1);
    // Baricentre of a given element
    Point<NDim> elemCentre(const PolygonalElementIterator<NDim>& t_elementID) const;

    //=======================================================
    // I/O
    //===================
    void plot(const std::string& t_title = "") const;
    void plotDetailed(const std::string& t_title = "") const;

    //=======================================================
    // Technicalities
    //=======================================================
public:
    // Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
    // See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // Removes the possibility of passing temporaries to the constructor
    PolygonalMesh(const Region<NDim>&& t_region, const Point<NDim>& t_relativeOrigin, const Point<NDim>& t_relativegVec, const CartIndex<NDim>& t_nSecSplits) = delete;

};

template<int NDim> std::ostream &operator<<(std::ostream &t_os, PolygonalMesh<NDim> const& t_mesh);
Plotter3D& operator << (Plotter3D& plotter, const PolygonalMesh<3>& t_mesh);
Plotter3D& operator << (Plotter3D& plotter, const PolygonalMesh<2>& t_mesh);
Plotter2D& operator << (Plotter2D& plotter, const PolygonalMesh<2>& t_mesh);
Plotter2D& operator << (Plotter2D& plotter, const PolygonalMesh<1>& t_mesh);

template<int NDim> bool operator==(const PolygonalMesh<NDim>& lhs, const PolygonalMesh<NDim>& rhs);

} // namespace Tortoise


#endif /* PolygonalMesh_hpp */
