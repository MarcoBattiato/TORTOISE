//
//  PolygonalMeshIterators.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 20/9/22.
//

#include <Geometry/Structured/DirectSpace/Polygonal/PolygonalMeshIterators.hpp>


namespace Tortoise {

CartIndex<2> buildCurrentcoord(const PolygonalMeshStructure<2>::const_iterator& currentStrip1stCoord, const PolygonalMeshStructure<2>::mapped_type::const_iterator&  currentStrip2ndCoord, int currentPosInStripe){
    CartIndex<2> toreturn;
    toreturn(0) = currentStrip1stCoord->first;
    toreturn(1) = currentStrip2ndCoord->first+currentPosInStripe;
    return toreturn;
};


PolygonalSectionIterator<2>::PolygonalSectionIterator(const PolygonalMeshStructure<2>& t_polyMeshStruct)
: polyMeshStruct(t_polyMeshStruct), currentStrip1stCoord(polyMeshStruct.begin()), currentStrip2ndCoord(currentStrip1stCoord->second.begin()), currentPosInStripe (0),
currentcoord(buildCurrentcoord(currentStrip1stCoord, currentStrip2ndCoord, currentPosInStripe)), unfinished (true) {};
PolygonalSectionIterator<2>::PolygonalSectionIterator(const PolygonalMesh<2>& t_polyMesh) : PolygonalSectionIterator(t_polyMesh.meshStructure) {}
PolygonalSectionIterator<2>::PolygonalSectionIterator(const PolygonalMesh<2>* t_polyMesh) : PolygonalSectionIterator(t_polyMesh->meshStructure) {}
PolygonalSectionIterator<2>& PolygonalSectionIterator<2>::operator++(){
        ++currentPosInStripe;
        if (currentPosInStripe == currentStrip2ndCoord->second[0]) {
            ++currentStrip2ndCoord;
            currentPosInStripe = 0;
            if (currentStrip2ndCoord == currentStrip1stCoord->second.end()) {
                ++currentStrip1stCoord;
                if ( currentStrip1stCoord != polyMeshStruct.end()) {
                    currentStrip2ndCoord = currentStrip1stCoord->second.begin();
                } else {
                    unfinished = false;
                    return *this;
                }
            }
        }
        currentcoord = buildCurrentcoord(currentStrip1stCoord, currentStrip2ndCoord, currentPosInStripe);
    return *this;
};


PolygonalElementIterator<2>::PolygonalElementIterator(const PolygonalMeshStructure<2>& t_polyMeshStruct): polyMeshStruct(t_polyMeshStruct), currentStrip1stCoord(polyMeshStruct.begin()), currentStrip2ndCoord(currentStrip1stCoord->second.begin()), currentPosInStripe (0),
currentSectCoord(buildCurrentcoord(currentStrip1stCoord, currentStrip2ndCoord, currentPosInStripe)), currentElemInSec(0), unfinished (true) {};
PolygonalElementIterator<2>::PolygonalElementIterator(const PolygonalMesh<2>& t_polyMesh) : PolygonalElementIterator(t_polyMesh.meshStructure) {}
PolygonalElementIterator<2>::PolygonalElementIterator(const PolygonalMesh<2>* t_polyMesh) : PolygonalElementIterator(t_polyMesh->meshStructure) {}
PolygonalElementIterator<2>& PolygonalElementIterator<2>::operator++(){
    ++currentElemInSec;
    if (currentElemInSec == Mesh<2>::nElementSection){
        currentElemInSec = 0;
        ++currentPosInStripe;
        if (currentPosInStripe == currentStrip2ndCoord->second[0]) {
            ++currentStrip2ndCoord;
            currentPosInStripe = 0;
            if (currentStrip2ndCoord == currentStrip1stCoord->second.end()) {
                ++currentStrip1stCoord;
                if ( currentStrip1stCoord != polyMeshStruct.end()) {
                    currentStrip2ndCoord = currentStrip1stCoord->second.begin();
                } else {
                    unfinished = false;
                    return *this;
                }
            }
        }
        currentSectCoord = buildCurrentcoord(currentStrip1stCoord, currentStrip2ndCoord, currentPosInStripe);
    }
    return *this;
};


} // namespace Tortoise
