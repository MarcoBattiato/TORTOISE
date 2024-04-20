////
////  Reticulation.hpp
////  TORTOISE
////
////  Created by Marco Battiato on 5/11/22.
////
//
//#ifndef Reticulation_hpp
//#define Reticulation_hpp
//
//#include <Geometry/GeometryCore/Geometry.hpp>
//#include <Geometry/Plotter/Plotter.hpp>
//#include <Geometry/Structured/DirectSpace/Region.hpp>
//
//namespace Tortoise {
//
//template <int NDim> class Reticulation {
//    
//    // Geometric Attributes
//public:
//    //*******************************
//    // Dimensionality constants
//    //*******************************
//    static constexpr int                            nVertElement = NDim+1;                  // Number of vertices in Element
//    static constexpr int                            nElementSection = factorial(NDim);      // Number of element per reference Section
//    static const VectorElement<NDim>                vertElemInRefSec;                       // Vertices of Ref Elements (They all share the first vertex at the origin)
//
//public: // DO NOT CHANCE THE ORDER!!! IT WILL BREAK THE CONSTRUCTOR!!!
//    //*******************************
//    // Reticulations features
//    //*******************************
//    const Region<NDim>*                             region;                                 // Region the mesh is part of
//
//    const Point<NDim>                               origin;                                 // origin of the Mesh
//    const ArrayPoint<NDim,NDim>                     gVec;                                   // Mesh sides vectors
//    const ArrayPoint<NDim,NDim>                     invgVec;                                // Inverse of Matrix formed by Mesh sides vectors
//
//    const CartIndex<NDim>                           nSecSplits;                             // Number of sectioning per dimensions
//
//    const int                                       numberSections;                         // Total number of sections
//    const int                                       numberElements;                         // Total number of elements
//
//    const ArrayPoint<NDim,NDim>                     dGVec;                                  // Sides of the sections
//    const ArrayPoint<NDim,NDim>                     invdGVec;                               // Inverse of Matrix formed by sides of the sections
//
//    const Real                                      elemVolume;                             // Area of the element
//    const LinTransform<NDim>                        massMatrix;                             // Mass matrix for linear functional representation
//    const LinTransform<NDim>                        invMassMatrix;                          // Inverse of Mass matrix for linear functional representation
//
//    const Point<NDim>                               relativeOrigin;                         // origin of the Mesh in fraction of Region
//    const Point<NDim>                               relativegVec;                           // Mesh sides vectors in fraction of Region
//    const Point<NDim>                               relativedGVec;                          // Mesh sides vectors in fraction of Region
//
//    const VectorLinearForm<NDim>                    basisFuncLF;                            // Linear forms corresponding to the basis functions
//    const std::vector<VectorLinearForm<NDim>>       basisFuncDerivLF;                       // Linear forms corresponding to the derivatives of the basis functions
//    const VecOfDataVector                           basisFuncDerivLF0thOrd;                 // Linear forms corresponding to the derivatives of the basis functions
//
//    // Methods
//public:
//    //=======================================================
//    // Constructors
//    // !!! Notice !!! gVec are the lattice vectors and NOT the end points of the mesh.
//    // The endpoints of the mesh are origin+gVec[1], origin+gVec[2], ...
//    // The origin has to be passed in relative position compared to the region's origin and gVec's
//    //===================
//    Reticulation(const Region<NDim>& t_region, const Point<NDim>& t_relativeOrigin, const Point<NDim>& t_relativegVec, const CartIndex<NDim>& t_nSecSplits);
//    Reticulation(const Region<NDim>* t_region, const Point<NDim>& t_relativeOrigin, const Point<NDim>& t_relativegVec, const CartIndex<NDim>& t_nSecSplits);
//    // In 1D a constructor accepting real numbers is provided
//    Reticulation(const Region<NDim>& t_region, const Real t_relativeOrigin, const Real t_relativegVec, const int t_nSecSplits) requires (NDim == 1);
//    Reticulation(const Region<NDim>* t_region, const Real t_relativeOrigin, const Real t_relativegVec, const int t_nSecSplits) requires (NDim == 1);
//    //===================
//
//    //********************************
//    //* Reticulation structure
//    //********************************
//    //=======================================================
//    // Reticulation Iterators and accessor
//    // Iterating through the mesh should be done as
//    // >>  for (auto elem = meshexmpl.elementIterator(); elem.unfinished; ++elem) {
//    // >>       something[elem]
//    // >>  }
//    // This allows for abstraction and also works for MeshSubsets that run only over some of the elements
//    // Notice that the iterator is always initialised to the first section, element, or vertex
//    //===================
//    // Indices conversions (more technical, the final user will not need these)
//    CartIndex<NDim> sectionID(const int t_elementID) const;
//    int             sectionID(const CartIndex<NDim>& SectionCartID) const;
//    int             elemInSectionID(const int t_elementID) const;
//    int             elemID(const CartIndex<NDim>& SectionCartID, const int elemInSect) const;
//    int             elemID(const std::array<int, NDim>& SectionCartID, const int elemInSect) const;
//
//    // Mesh indices
//    RetirculationSectionIndex<NDim>         sectionBeginIndex() const;
//    ReticulationElementIndex<NDim>          elementBeginIndex() const;
//    VertexInElemIndex<NDim>                 vertexInElemBeginIndex() const;
//    RetirculationSectionIndex<NDim>         sectionEndIndex() const;
//    ReticulationElementIndex<NDim>          elementEndIndex() const;
//    VertexInElemIndex<NDim>                 vertexInElemEndIndex() const;
//
//    // Mesh iterators
//    ReticulationSectionIterator<NDim>       sectionBeginIterator() const;
//    ReticulationElementIterator<NDim>       elementBeginIterator() const;
//    VertexInElemIterator<NDim>              vertexInElemBeginIterator() const;
//    ReticulationSectionIterator<NDim>       sectionEndIterator() const;
//    ReticulationElementIterator<NDim>       elementEndIterator() const;
//    VertexInElemIterator<NDim>              vertexInElemEndIterator() const;
//
//    
//    //=======================================================
//    // Reticulation geometry
//    //===================
//    // Reticulation Points from section ID
//    Point<NDim> operator()(const std::array<int, NDim>& SectionCartID) const;                                            // Mesh point in cartesian coordinates
//    Point<NDim> operator()(const SectionIterator<NDim>& t_sectionID) const;                                              // Mesh point from section ID
//    Point<NDim> operator()(const CartIndex<NDim>& t_sectionID) const;                                                    // Mesh point from section ID
//    // Vertex from element ID and vertex in element ID
//    Point<NDim> operator()(const MeshElementIterator<NDim>& t_elementID, const VertexInElemIterator<NDim>& t_nVertex) const;
//    Point<NDim> operator()(const int t_elementID, const int t_nVertex) const;
//    Point<NDim> operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const VertexInElemIterator<NDim>& t_nVertex) const;
//    // Point from element ID and coordinates in reference element
//    Point<NDim> operator()(const MeshElementIterator<NDim>& t_elementID, const Point<NDim>& refCoord) const;
//    Point<NDim> operator()(const MeshElementIterator<NDim>& t_elementID, const Real refCoord) const requires (NDim == 1);
//    Point<NDim> operator()(const int t_elementID, const Point<NDim>& refCoord) const;
//    Point<NDim> operator()(const int t_elementID, const Real refCoord) const requires (NDim == 1);
//    Point<NDim> operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const Point<NDim>& refCoord) const;
//    Point<NDim> operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const Real refCoord) const requires (NDim == 1);
//    // Baricentre of a given element
//    Point<NDim> elemCentre(const MeshElementIterator<NDim>& t_elementID) const;
//    Point<NDim> elemCentre(const MeshSubsetElementIterator<NDim>& t_elementID) const;
//    Point<NDim> elemCentre(const int t_elementID) const;
//
//    
//    
//    
//    
//    
//    //=======================================================
//    // Technicalities
//    //=======================================================
//    public:
//    // Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
//    // See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
//        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//    // Removes the possibility of passing temporaries to the constructor
//    Reticulation(const Region<NDim>&& t_region, const Point<NDim>& t_relativeOrigin, const Point<NDim>& t_relativegVec, const CartIndex<NDim>& t_nSecSplits) = delete;
//    Reticulation(const Region<NDim>&& t_region, const Real t_relativeOrigin, const Real t_relativegVec, const int t_nSecSplits) requires (NDim ==1) = delete;
//
//}
//
//template <int NDim> class RetirculationSectionIndex {
//    CartIndex<NDim>     cartID;
//    CartIndex<NDim>     maxCartID;
//public:
//    RetirculationSectionIndex(const CartIndex<NDim>& cartID, const CartIndex<NDim>& maxCartID);
//    SectionIterator& operator++();
//    SectionIterator& operator++(int);
//    
//    
//public:
//// Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
//// See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//}
//
//template <int NDim> class RetirculationElementIndex{
//    
//}
//
//template <int NDim> class RetirculationElemInSecIndex{
//    
//}
//
//template <int NDim> class RetirculationVertInElemIndex{
//    
//}
//
//
//
//
//} // namespace Tortoise
//
//
//#endif /* Reticulation_hpp */
