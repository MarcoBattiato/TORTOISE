//
//  PolygonalMeshIterators.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 20/9/22.
//

#ifndef PolygonalMeshIterators_hpp
#define PolygonalMeshIterators_hpp

#include <Geometry/Structured/DirectSpace/Polygonal/PolygonalMesh.hpp>

namespace Tortoise {

//*******************************
// Section Iterator
//*******************************
template <int NDim> class PolygonalSectionIterator;

template<> class PolygonalSectionIterator<2> {
public:
    const PolygonalMeshStructure<2>&                        polyMeshStruct;
    PolygonalMeshStructure<2>::const_iterator               currentStrip1stCoord;
    PolygonalMeshStructure<2>::mapped_type::const_iterator  currentStrip2ndCoord;
    int                                                     currentPosInStripe;
    CartIndex<2>                                            currentcoord;
    bool                                                    unfinished;
    
public:
    PolygonalSectionIterator(const PolygonalMesh<2>& t_polyMesh);
    PolygonalSectionIterator(const PolygonalMesh<2>* t_polyMesh);
    PolygonalSectionIterator(const PolygonalMeshStructure<2>& t_polyMeshStruct);
    //=======================================================
    // Increment and assignment
    //===================
    // Prefix increment
    PolygonalSectionIterator& operator++();
        
public:
    // Removes the possibility of passing temporaries to the constructor
    PolygonalSectionIterator(const PolygonalMesh<2>&& t_polyMesh) = delete;
    PolygonalSectionIterator(const PolygonalMeshStructure<2>&& t_polyMeshStruct) = delete;
    
    // See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};

//*******************************
// Element Iterators
//*******************************
template <int NDim> class PolygonalElementIterator;

template<> class PolygonalElementIterator<2> {
public:
    const PolygonalMeshStructure<2>&                        polyMeshStruct;
    PolygonalMeshStructure<2>::const_iterator               currentStrip1stCoord;
    PolygonalMeshStructure<2>::mapped_type::const_iterator  currentStrip2ndCoord;
    int                                                     currentPosInStripe;
    CartIndex<2>                                            currentSectCoord;
    int                                                     currentElemInSec;
    bool                                                    unfinished;
    
public:
    PolygonalElementIterator(const PolygonalMesh<2>& t_polyMesh);
    PolygonalElementIterator(const PolygonalMesh<2>* t_polyMesh);
    PolygonalElementIterator(const PolygonalMeshStructure<2>& t_polyMeshStruct);
    //=======================================================
    // Increment and assignment
    //===================
    // Prefix increment
    PolygonalElementIterator& operator++();
        
public:
    // Removes the possibility of passing temporaries to the constructor
    PolygonalElementIterator(const PolygonalMesh<2>&& t_polyMesh) = delete;
    PolygonalElementIterator(const PolygonalMeshStructure<2>&& t_polyMeshStruct) = delete;
    
    // See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};



} // namespace Tortoise

#endif /* PolygonalMeshIterators_hpp */
