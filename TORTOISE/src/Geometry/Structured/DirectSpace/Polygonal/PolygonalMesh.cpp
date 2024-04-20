//
//  PolygonalMesh.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 16/9/22.
//

#include <Geometry/Structured/DirectSpace/Polygonal/PolygonalMesh.hpp>
#include <Geometry/Structured/DirectSpace/Polygonal/PolygonalMeshIterators.hpp>

#include <vector>

using namespace Tortoise::GeometryCore;

namespace Tortoise {

void addNextElementVertexPoly(const std::vector<std::vector<Real>>& CurrentIncompleteVertices, std::vector<std::vector<std::vector<Real>>>& ListElementVertices, const int dim) {
    if (CurrentIncompleteVertices.size()==dim+1) {
        ListElementVertices.push_back(CurrentIncompleteVertices);
    }
    else {
        for (int i=0; i<dim; ++i){
            if((CurrentIncompleteVertices.back())[i]<0.001){
                std::vector<std::vector<Real>> newlistVertices(CurrentIncompleteVertices);
                newlistVertices.push_back(CurrentIncompleteVertices.back());
                (newlistVertices.back())[i]=1.0;
                addNextElementVertexPoly(newlistVertices,ListElementVertices,dim);
            }
        }
    }
};
template<int NDim>  VectorElement<NDim>           calculatevertElemInRefSecPoly(){
    VectorElement<NDim> toreturn;
    std::vector<std::vector<std::vector<Real>>> ListElementVertices;
    std::vector<Real> start(NDim,0.0);
    addNextElementVertexPoly({start}, ListElementVertices, NDim);
    for (int i=0; i<ListElementVertices.size(); ++i){
        Element<NDim> temporary;
        for (int col=0; col<ListElementVertices[i].size(); ++col){
            for (int row=0; row<ListElementVertices[i][col].size(); ++row){
                temporary(row,col)= ListElementVertices[i][col][row];
            }
        }
        toreturn.push_back(temporary);
    }
    return toreturn;
}
template<int NDim> const VectorElement<NDim>       PolygonalMesh<NDim>::vertElemInRefSec = calculatevertElemInRefSecPoly<NDim>();



template<int NDim> PolygonalMeshStructure<NDim> buildPolygonalMeshStructure(const Mesh<NDim>& t_meshFull, const Eigen::Matrix<bool,1,Eigen::Dynamic>& t_whichSections);

//template<> PolygonalMeshStructure<1> buildPolygonalMeshStructure(const Mesh<1>& t_meshFull, const Eigen::Matrix<bool,1,Eigen::Dynamic>& t_whichSections) {
//
//};
template<> PolygonalMeshStructure<2> buildPolygonalMeshStructure(const Mesh<2>& t_meshFull, const Eigen::Matrix<bool,1,Eigen::Dynamic>& t_whichSections) {
    PolygonalMeshStructure<2> meshStructure;
    auto vectorBools = t_whichSections.reshaped(t_meshFull.nSecSplits(0),t_meshFull.nSecSplits(1)).eval();
    int accumulated = 0;
    for(int i = 0; i < t_meshFull.nSecSplits(0); ++i){
        int j = 0;
        while (j < t_meshFull.nSecSplits(1)) {
            while ( (j < t_meshFull.nSecSplits(1)) && (vectorBools(i,j) != 1)) {
                ++j;
            }
            int length = 1;
            if ( j != t_meshFull.nSecSplits(1)) { // Then a 1 has been found
                while ( ((j+length) < t_meshFull.nSecSplits(1)) && (vectorBools(i,j+length) != 0)) { ++length;}
                meshStructure[i][j] = {length, accumulated };
                accumulated += length;
            }
            j += length;
        }
    }
    return meshStructure;
};

template<int NDim> int buildNumberSections(const PolygonalMeshStructure<NDim>& t_structure);
template<> int buildNumberSections<2>(const PolygonalMeshStructure<2>& t_structure) {
    return (--(--t_structure.end())->second.end())->second[1] + (--(--t_structure.end())->second.end())->second[0];
};

template<int NDim> PolygonalMesh<NDim>::PolygonalMesh(const Mesh<NDim>& t_meshFull, const Eigen::Matrix<bool,1,Eigen::Dynamic>& t_whichSections)
: region (t_meshFull.region),
  meshFull(&t_meshFull),
  meshStructure(buildPolygonalMeshStructure(t_meshFull,t_whichSections)),
  numberSections(buildNumberSections<NDim>(meshStructure)),
  numberElements(numberSections * nElementSection),
  dataVectorDim(numberElements * nVertElement),
  dGVec(t_meshFull.dGVec),
  invdGVec(t_meshFull.invdGVec),
  elemVolume(t_meshFull.elemVolume),
  massMatrix(t_meshFull.massMatrix),
  invMassMatrix(t_meshFull.invMassMatrix)
{
    assert( t_whichSections.sum() > 0 && "The polygonal Mesh is empty!");
};

template<int NDim> PolygonalMesh<NDim>::PolygonalMesh(const Mesh<NDim>* t_meshFull, const Eigen::Matrix<bool,1,Eigen::Dynamic>& t_whichSections) :
PolygonalMesh(* t_meshFull, t_whichSections){};
template<int NDim> PolygonalMesh<NDim>::PolygonalMesh(const Function<NDim> &t_function, const std::function<bool(Real)>& boolFunctionToEvaluate, const RequirementType requirementType):
PolygonalMesh(t_function.mesh, t_function.evaluateBoolCondOnSections(boolFunctionToEvaluate, requirementType)) {};

//===================

//********************************
//* Mesh structure
//********************************
// Mesh iterators
template<int NDim> PolygonalSectionIterator<NDim> PolygonalMesh<NDim>::sectionIterator() const {
    return PolygonalSectionIterator<NDim>(this);
};
template<int NDim> PolygonalElementIterator<NDim> PolygonalMesh<NDim>::elementIterator() const {
    return PolygonalElementIterator<NDim>(this);
};
template<int NDim> VertexInElemIterator<NDim> PolygonalMesh<NDim>::vertexInElemIterator() const {
    return VertexInElemIterator<NDim>();
};

//=======================================================
// Mesh geometry
//===================
// Mesh Points from section ID
template<int NDim> Point<NDim> PolygonalMesh<NDim>::operator()(const PolygonalSectionIterator<NDim>& t_sectionID) const{
    Point<NDim> toreturn(origin);
    for (int i=0; i<NDim; ++i){
        toreturn += dGVec.col(i)*static_cast<Real>(t_sectionID.currentcoord[i]);
    }
    return toreturn;
}
// Vertex from element ID and vertex in element ID
template<int NDim> Point<NDim> PolygonalMesh<NDim>::operator()(const PolygonalElementIterator<NDim>& t_elementID, const VertexInElemIterator<NDim>& t_nVertex) const {
    return origin + dGVec * (t_elementID.currentSectCoord.template cast<Real>() + vertElemInRefSec[t_elementID.currentElemInSec].col(t_nVertex.currentVertID));
};
// Point from element ID and coordinates in reference element
template<int NDim> Point<NDim> PolygonalMesh<NDim>::operator()(const PolygonalElementIterator<NDim>& t_elementID, const Point<NDim>& refCoord) const {
    return origin + dGVec * (t_elementID.currentSectCoord.template cast<Real>() + vertElemInRefSec[t_elementID.currentElemInSec].template block<NDim,NDim>(0,1) * refCoord);
};
template<int NDim> Point<NDim> PolygonalMesh<NDim>::operator()(const PolygonalElementIterator<NDim>& t_elementID, const Real refCoord) const requires (NDim == 1) {
    return origin + dGVec * (t_elementID.currentSectCoord.template cast<Real>() + vertElemInRefSec[t_elementID.currentElemInSec].template block<NDim,NDim>(0,1) * Point<NDim>::Constant(refCoord));
};
// Baricentre of a given element
template<> Point<2> PolygonalMesh<2>::elemCentre(const PolygonalElementIterator<2>& t_elementID) const {
    return (*this)(t_elementID,Point<2>({0.333333333,0.333333333}));
}


//=======================================================
// I/O
//===================

Plotter3D& operator << (Plotter3D& plotter, const PolygonalMesh<2>& t_mesh){
    plotter << OPENSTREAM << POLYGONMODE;
    for (auto elem = t_mesh.elementIterator(); elem.unfinished; ++elem){
        ArrayPoint<2,2+1> triang; int i=0;
        for (auto vertex = t_mesh.vertexInElemIterator(); vertex.unfinished ; ++vertex){
            triang.col(i++) = t_mesh(elem,vertex);
        }
        plotter << triang;
    }
    return plotter << CLOSESTREAM ;
};
Plotter2D& operator << (Plotter2D& plotter, const PolygonalMesh<2>& t_mesh){
    plotter << OPENSTREAM << POLYGONMODE;
    for (auto elem = t_mesh.elementIterator(); elem.unfinished; ++elem){
        ArrayPoint<2,2+1> triang; int i=0;
        for (auto vertex = t_mesh.vertexInElemIterator(); vertex.unfinished ; ++vertex){
            triang.col(i++) = t_mesh(elem,vertex);
        }
        plotter << triang;
    }
    return plotter << CLOSESTREAM ;
};

template<> void PolygonalMesh<2>::plot(const std::string& t_title) const {
    plotter3d << *region << NOLABEL << COLOR("blue") ;
    plotter3d << *this  << NOLABEL << COLOR("red");
    if (t_title!= "") plotter3d << PLOTTITLE(t_title);
    plotter3d << PLOT;
}


// Explicit template instantiations (For more info see: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl)

//template class PolygonalMesh<1>;
template class PolygonalMesh<2>;
//template class PolygonalMesh<3>;

//template  std::ostream &operator<<(std::ostream &os, const PolygonalMesh<1> & mesh);
//template  std::ostream &operator<<(std::ostream &os, const PolygonalMesh<2> & mesh);
//template  std::ostream &operator<<(std::ostream &os, const PolygonalMesh<3> & mesh);

//template bool operator==(const PolygonalMesh<1>& lhs, const PolygonalMesh<1>& rhs);
//template bool operator==(const PolygonalMesh<2>& lhs, const PolygonalMesh<2>& rhs);
//template bool operator==(const PolygonalMesh<3>& lhs, const PolygonalMesh<3>& rhs);

} // namespace Tortoise
