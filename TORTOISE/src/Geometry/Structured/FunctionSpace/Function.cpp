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
//  Function.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//

#include <Geometry/Structured/FunctionSpace/Function.hpp>

#include <cassert>

namespace Tortoise {

//=======================================================
// Constructors
//===================

template<int NDim> Function<NDim>::Function(const Mesh<NDim>& t_mesh): mesh(&t_mesh), vec(DataVector::Zero(t_mesh.dataVectorDim)){};
template<int NDim> Function<NDim>::Function(const Mesh<NDim>& t_mesh, const Real t_value): mesh(&t_mesh), vec(DataVector::Zero(t_mesh.dataVectorDim)){
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), mesh->nVertElement, mesh->numberElements);
    elementview.row(0).array() = t_value;
};
template<int NDim> Function<NDim>::Function(const Mesh<NDim>& t_mesh, const DataVector& t_nodalValues): mesh(&t_mesh), vec(t_nodalValues){
    assert(vec.size()==t_nodalValues.size());
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), mesh->nVertElement, mesh->numberElements);
    elementview.block(1,0,mesh->nVertElement-1,mesh->numberElements).rowwise() -= elementview.row(0);
};
template<int NDim> Function<NDim>::Function(const Mesh<NDim>& t_mesh, const std::function<Real(Point<NDim>)>& t_f): mesh(&t_mesh), vec(t_mesh.dataVectorDim){
    for (auto elem = mesh->elementIterator(); elem.unfinished ; ++elem) {
        for (auto vert = mesh->vertexInElemIterator(); vert.unfinished; ++vert){
            vec(vert.currentVertID + mesh->nVertElement * elem.currentElementID) = t_f(mesh->operator()(elem,vert));
        }
    }
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), mesh->nVertElement, mesh->numberElements);
    elementview.block(1,0,mesh->nVertElement-1,mesh->numberElements).rowwise() -= elementview.row(0);
}

template<int NDim> Function<NDim>::Function(const Mesh<NDim>* t_mesh): Function(*t_mesh) {}
template<int NDim> Function<NDim>::Function(const Mesh<NDim>* t_mesh, const Real t_value): Function(*t_mesh, t_value) {}
template<int NDim> Function<NDim>::Function(const Mesh<NDim>* t_mesh, const DataVector& t_nodalValues): Function(*t_mesh, t_nodalValues) {}
template<int NDim> Function<NDim>::Function(const Mesh<NDim>* t_mesh, const std::function<Real(Point<NDim>)>& t_f): Function(*t_mesh, t_f) {}

template<int NDim> Function<NDim>::Function(const std::function<Real(Real)>& t_f, const Function<NDim>& t_other): mesh(t_other.mesh), vec(t_other.vec){
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), mesh->nVertElement, mesh->numberElements);
    elementview.block(1,0,mesh->nVertElement-1,mesh->numberElements).rowwise() += elementview.row(0);
    for (auto elem = mesh->elementIterator(); elem.unfinished ; ++elem) {
        for (auto vert = mesh->vertexInElemIterator(); vert.unfinished; ++vert){
            vec(vert.currentVertID + mesh->nVertElement * elem.currentElementID) = t_f(vec(vert.currentVertID + mesh->nVertElement * elem.currentElementID));
        }
    }
    elementview.block(1,0,mesh->nVertElement-1,mesh->numberElements).rowwise() -= elementview.row(0);
};
template<int NDim> Function<NDim>::Function(const std::function<Real(std::vector<Real>)>& t_f, const std::vector<Function<NDim>>& t_others): mesh(t_others[0].mesh), vec(t_others[0].vec){
    for(auto it = ++t_others.begin(); it != t_others.end(); it++) { assert(t_others[0].mesh == it->mesh); }
    // The most convenience form is to transform the functions into a nodal representation
    // However that cannot be done since they have to be const
    // For that reason a local copy is made
    // This is not the most efficient way, but given that this constructor is going to be called very few times it does not matter
    std::vector<Function<NDim>> localFunctions;
    for(auto it = t_others.begin(); it != t_others.end(); it++) {
        localFunctions.emplace_back(*it);
        localFunctions.back().toNodalFromLinForm();
    }
    std::vector<Real> functionValues;
    functionValues.resize(t_others.size());
    for (auto elem = mesh->elementIterator(); elem.unfinished ; ++elem) {
        for (auto vert = mesh->vertexInElemIterator(); vert.unfinished; ++vert){
            for (int i=0; i<functionValues.size(); ++i){ functionValues[i] = localFunctions[i].vec(vert.currentVertID + mesh->nVertElement * elem.currentElementID);}
            vec(vert.currentVertID + mesh->nVertElement * elem.currentElementID) = t_f(functionValues);
        }
    }
    toLinFormFromNodal();
};

template<int NDim> Function<NDim>::Function(const Function<NDim>& other): mesh(other.mesh),vec(other.vec) {};
template<int NDim> Function<NDim>::Function(Function<NDim>&& other) : mesh(other.mesh), vec(std::move(other.vec)) {}
template<int NDim> void swap(Function<NDim>& first, Function<NDim>& second){
    assert(first.mesh == second.mesh);
    first.vec.swap(second.vec);
}


//********************************
//* Arithmetic
//********************************

template<int NDim> Function<NDim>& Function<NDim>::operator=(const Function<NDim>& other){
    assert(mesh==other.mesh);
    vec = other.vec;
    return *this;
}
template<int NDim> Function<NDim>& Function<NDim>::operator=(Function<NDim>&& other){  // Swap assignment
    assert(mesh==other.mesh);
//    vec.swap(other.vec);
    vec = std::move(other.vec);
    return *this;
}
template<int NDim> Function<NDim>& Function<NDim>::operator+=(const Function<NDim>& other){
    assert(mesh==other.mesh);
    vec+=other.vec;
    return *this;
};
template<int NDim> Function<NDim>& Function<NDim>::operator-=(const Function<NDim>& other){
    assert(mesh==other.mesh);
    vec-=other.vec;
    return *this;
};
template<int NDim> Function<NDim>& Function<NDim>::operator*=(const Function<NDim>& other){
    assert(mesh==other.mesh);
    toNodalFromLinForm();                                // Converts to nodal representation and stores in vec
    auto otherNodalRepr = other.toNodalFromLinForm();    // Converts other to nodal representation and stores in new matrix (since other is const)
    elementview().array() *= otherNodalRepr.array();// Executes the operation by treating the matrices as array (in Eigen sense)
    toLinFormFromNodal();                              // Converts back to lienar form representation and stores in vec
    return *this;
};
template<int NDim> Function<NDim>& Function<NDim>::operator/=(const Function<NDim>& other){
    assert(mesh==other.mesh);
    toNodalFromLinForm();                                // Converts to nodal representation and stores in vec
    auto otherNodalRepr = other.toNodalFromLinForm();    // Converts other to nodal representation and stores in new matrix (since other is const)
    elementview().array() /= otherNodalRepr.array();// Executes the operation by treating the matrices as array (in Eigen sense)
    toLinFormFromNodal();                              // Converts back to lienar form representation and stores in vec
    return *this;
};

template<int NDim> Function<NDim> Function<NDim>::operator-() const {
    Function<NDim> toreturn (*this);
    toreturn.vec *= -1.0;
    return toreturn;
}

template<int NDim> Function<NDim>& Function<NDim>::operator=(const Real scalar){
    vec = DataVector::Zero(mesh->dataVectorDim);
    elementview().col(0).array() = scalar;
    return *this;
}
template<int NDim> Function<NDim>& Function<NDim>::operator+=(const Real scalar){
    elementview().col(0).array() += scalar;
    return *this;
}
template<int NDim> Function<NDim>& Function<NDim>::operator-=(const Real scalar){
    elementview().col(0).array() -= scalar;
    return *this;
}
template<int NDim> Function<NDim>& Function<NDim>::operator*=(const Real scalar){
    vec *= scalar;
    return *this;
}
template<int NDim> Function<NDim>& Function<NDim>::operator/=(const Real scalar){
    vec /= scalar;
    return *this;
}
template<int NDim> Function<NDim> operator/(const Real lhs, Function<NDim> rhs) {
    rhs.toNodalFromLinForm();
    rhs.vec = rhs.vec.cwiseInverse();
    rhs.toLinFormFromNodal();
    return lhs*rhs;
}

template<int NDim> Function<NDim>& Function<NDim>::operator+=(const FunctionElement<NDim>& other){
    assert(mesh==other.mesh);
    vec.template segment<NDim+1>(other.elemID*(mesh->nVertElement)) += other.vec;
    return *this;
}
template<int NDim> Function<NDim>& Function<NDim>::operator-=(const FunctionElement<NDim>& other){
    assert(mesh==other.mesh);
    vec.template segment<NDim+1>(other.elemID*(mesh->nVertElement)) -= other.vec;
    return *this;
}

template<int NDim> Function<NDim>& Function<NDim>::operator=(const std::function<Real(Point<NDim>)>& other){
    for (auto elem = mesh->elementIterator(); elem.unfinished ; ++elem) {
        for (auto vert = mesh->vertexInElemIterator(); vert.unfinished; ++vert){
            vec(vert.currentVertID + mesh->nVertElement * elem.currentElementID) = other(mesh->operator()(elem,vert));
        }
    }
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), mesh->nVertElement, mesh->numberElements);
    elementview.block(1,0,mesh->nVertElement-1,mesh->numberElements).rowwise() -= elementview.row(0);
    return *this;
}
template<int NDim> Function<NDim>& Function<NDim>::operator+=(const std::function<Real(Point<NDim>)>& other){
    (*this) += Function<NDim>(mesh, other);
    return *this;
}
template<int NDim> Function<NDim>& Function<NDim>::operator-=(const std::function<Real(Point<NDim>)>& other){
    (*this) -= Function<NDim>(mesh, other);
    return *this;
}
template<int NDim> Function<NDim>& Function<NDim>::operator*=(const std::function<Real(Point<NDim>)>& other){
    (*this) *= Function<NDim>(mesh, other);
    return *this;
}
template<int NDim> Function<NDim>& Function<NDim>::operator/=(const std::function<Real(Point<NDim>)>& other){
    (*this) /= Function<NDim>(mesh, other);
    return *this;
}
template<int NDim> Function<NDim> operator/(const std::function<Real(Point<NDim>)>& lhs, Function<NDim> rhs){
    return Function<NDim>(rhs.mesh, lhs)/rhs;
}

//********************************
//* Mesh structure
//********************************
//=======================================================
// Mesh Iterators
//===================


template<int NDim> SectionIterator<NDim> Function<NDim>::sectionIterator() const {
    return mesh->sectionIterator();
}
template<int NDim> CartIndex<NDim> Function<NDim>::randomSection() const {
    return mesh->randomSection();
}
template<int NDim> MeshElementIterator<NDim> Function<NDim>::elementIterator() const {
    return mesh->elementIterator();
}
template<int NDim> int Function<NDim>::randomElement() const {
    return mesh->randomElement();
}
template<int NDim> VertexInElemIterator<NDim> Function<NDim>::vertexInElemIterator() const {
    return VertexInElemIterator<NDim>();
}

//********************************
//* Local Evaluation
//********************************
template<int NDim> FunctionElement<NDim> Function<NDim>::operator()(const int t_elementID) const{
    FunctionElement<NDim> toreturn (*mesh, t_elementID);
    toreturn.vec = vec.template segment<NDim+1>(t_elementID *(mesh->nVertElement));
    return toreturn;
}
template<int NDim> FunctionElement<NDim> Function<NDim>::operator()(const MeshElementIterator<NDim>& t_elementID) const{
    return (*this)(t_elementID.currentElementID);
}
template<int NDim> FunctionElement<NDim> Function<NDim>::operator()(const MeshSubsetElementIterator<NDim>& t_elementID) const{
    return (*this)(t_elementID.currentElementID);
}

template<int NDim> Real Function<NDim>::operator()(const int t_elementID, const Point<NDim>& t_point_ref) const{
    return vec.template segment<NDim>(t_elementID*(mesh->nVertElement)+1)*t_point_ref+vec(t_elementID*(mesh->nVertElement));
}
template<int NDim> Real Function<NDim>::operator()(const int t_elementID, const Real t_point_ref) const requires (NDim == 1){
    return vec.template segment<NDim>(t_elementID*(mesh->nVertElement)+1)*Point<NDim>::Constant(t_point_ref)+vec(t_elementID*(mesh->nVertElement));
}
template<int NDim> Real Function<NDim>::operator()(const int t_elementID, const int t_nVertex) const{
    return vec.template segment<NDim>(t_elementID*(mesh->nVertElement)+1)*refNodes<NDim>.col(t_nVertex)+vec(t_elementID*(mesh->nVertElement));
}
template<int NDim> Real Function<NDim>::operator()(const MeshElementIterator<NDim>& t_elementID, const Point<NDim>& t_point_ref) const {
    return apply(vec.template segment<NDim+1>(t_elementID.currentElementID *(mesh->nVertElement)),t_point_ref)(0,0);
}
template<int NDim> Real Function<NDim>::operator()(const MeshElementIterator<NDim>& t_elementID, const Real t_point_ref) const requires (NDim == 1){
    return apply(vec.template segment<NDim+1>(t_elementID.currentElementID *(mesh->nVertElement)),Point<NDim>::Constant(t_point_ref))(0,0);
}
template<int NDim> Real Function<NDim>::operator()(const MeshElementIterator<NDim>& t_elementID, const VertexInElemIterator<NDim>& t_nVertex) const {
    return apply(vec.template segment<NDim+1>(t_elementID.currentElementID *(mesh->nVertElement)), refNodes<NDim>.col(t_nVertex.currentVertID))(0,0);
}
template<int NDim> Real Function<NDim>::operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const Point<NDim>& t_point_ref) const {
    return apply(vec.template segment<NDim+1>(t_elementID.currentElementID *(mesh->nVertElement)),t_point_ref)(0,0);
}
template<int NDim> Real Function<NDim>::operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const Real t_point_ref) const requires (NDim == 1){
    return apply(vec.template segment<NDim+1>(t_elementID.currentElementID *(mesh->nVertElement)),Point<NDim>::Constant(t_point_ref))(0,0);
}
template<int NDim> Real Function<NDim>::operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const VertexInElemIterator<NDim>& t_nVertex) const {
    return apply(vec.template segment<NDim+1>(t_elementID.currentElementID *(mesh->nVertElement)),refNodes<NDim>.col(t_nVertex.currentVertID))(0,0);
}


//********************************
//* Other operations
//********************************
            
template<int NDim> Real Function<NDim>::elementMax(const MeshElementIterator<NDim>& t_elementID) const {
    FunctionElement<NDim> elemFun ((*this)(t_elementID));
    elemFun.vec.template segment<NDim>(1).array() += elemFun.vec(0);
    return elemFun.vec.maxCoeff();
}
template<int NDim> Real Function<NDim>::elementMin(const MeshElementIterator<NDim>& t_elementID) const{
    FunctionElement<NDim> elemFun ((*this)(t_elementID));
    elemFun.vec.template segment<NDim>(1).array() += elemFun.vec(0);
    return elemFun.vec.minCoeff();
}
template<int NDim> Real Function<NDim>::elementMax(const MeshSubsetElementIterator<NDim>& t_elementID) const {
    FunctionElement<NDim> elemFun ((*this)(t_elementID));
    elemFun.vec.template segment<NDim>(1).array() += elemFun.vec(0);
    return elemFun.vec.maxCoeff();
}
template<int NDim> Real Function<NDim>::elementMin(const MeshSubsetElementIterator<NDim>& t_elementID) const{
    FunctionElement<NDim> elemFun ((*this)(t_elementID));
    elemFun.vec.template segment<NDim>(1).array() += elemFun.vec(0);
    return elemFun.vec.minCoeff();
}
template<int NDim> Real Function<NDim>::elementMax(const int t_elementID) const {
    FunctionElement<NDim> elemFun ((*this)(t_elementID));
    elemFun.vec.template segment<NDim>(1).array() += elemFun.vec(0);
    return elemFun.vec.maxCoeff();
}
template<int NDim> Real Function<NDim>::elementMin(const int t_elementID) const{
    FunctionElement<NDim> elemFun ((*this)(t_elementID));
    elemFun.vec.template segment<NDim>(1).array() += elemFun.vec(0);
    return elemFun.vec.minCoeff();
}
template<int NDim> Real Function<NDim>::max() const{
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> nodalRepr = toNodalFromLinForm();  // Converts to nodal representation and stores in vec
    return nodalRepr.maxCoeff();
}
template<int NDim> Real Function<NDim>::min() const{
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> nodalRepr = toNodalFromLinForm();  // Converts to nodal representation and stores in vec
    return nodalRepr.minCoeff();
}
template<int NDim> std::pair<std::vector<int>,std::vector<Point<NDim>>> Function<NDim>::zeros() const requires (NDim == 1) {
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> nodalRepr = toNodalFromLinForm();
    std::vector<int> elems;
    std::vector<Point<NDim>>  zeroPos;
    for (int el = 0; el < mesh->numberElements; ++el){
        if(nodalRepr(el,0)*nodalRepr(el,1)<=0.){
            elems.emplace_back(el);
            if (nodalRepr(el,0)  - nodalRepr(el,1) != 0.) {
                zeroPos.emplace_back(Point<NDim>::Constant(nodalRepr(el,0)/(nodalRepr(el,0)  - nodalRepr(el,1))));
            } else {
                zeroPos.emplace_back(Point<NDim>::Constant(0.));
            }
        }
        if (el+1 < mesh->numberElements) {
            if(nodalRepr(el,1)*nodalRepr(el+1,0)<0.){
                elems.emplace_back(el);
                zeroPos.emplace_back(Point<NDim>::Constant(1.));
            }
        }
    }
    return std::make_pair(elems, zeroPos);
}
    
template<int NDim> Real Function<NDim>::integrate() const {
    auto elView = elementview();
    return mesh->elemVolume*(elView.sum()+NDim*elView.col(0).sum())/(NDim+1);
}
template<int NDim> Real Function<NDim>::integrate(const Function<NDim>& rhs) const{
    return mesh->elemVolume*((apply(elementview(),gaussPoints<NDim,2>).array()*apply(rhs.elementview(),gaussPoints<NDim,2>).array()).matrix()*(gaussWeigths<NDim,2>.asDiagonal())).sum();
}
template<int NDim> Real Function<NDim>::integrate(const Function<NDim>& rhs1, const Function<NDim>& rhs2) const{
    return mesh->elemVolume*((apply(elementview(),gaussPoints<NDim,3>).array()*apply(rhs1.elementview(),gaussPoints<NDim,3>).array()*apply(rhs2.elementview(),gaussPoints<NDim,3>).array()).matrix()*(gaussWeigths<NDim,3>.asDiagonal())).sum();
}

template<int NDim> Function<NDim> Function<NDim>::derivative(const int direction) const {
    assert(direction >= 0 && direction < NDim);
    Function<NDim> deriv(mesh);
    DataVector temp (vec.array() * mesh->basisFuncDerivLF0thOrd[direction].array());
    deriv.elementview().col(0) = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>::Map(temp.data(),Mesh<NDim>::nVertElement,mesh->numberElements).colwise().sum();
    return deriv;
}

template<int NDim> Function<NDim>& Function<NDim>::applyInverseMass() {
    elementview() = elementview() * mesh->invMassMatrix;
    return *this;
}

//********************************
//* I/O
//********************************

template<> void Function<2>::plot(const std::string& t_title) const {
    plotter3d  << *mesh << NOLABEL << COLOR("red") << *this << NOLABEL << COLOR("red") ;
    if (t_title!= "") plotter3d << PLOTTITLE(t_title);
    plotter3d << PLOT;
}
template<> void Function<1>::plot(const std::string& t_title) const {
    plotter2d  << *mesh << NOLABEL << COLOR("red") << *this << NOLABEL << COLOR("red") ;
    if (t_title!= "") plotter2d << PLOTTITLE(t_title);
    plotter2d << PLOT;
}

template<int NDim> std::pair<std::vector<Real>, std::vector<Real>> sideLineIntersection(const ArrayPoint<NDim, NDim>& tetra, const ArrayPoint<NDim, 2>& line);
template<> std::pair<std::vector<Real>, std::vector<Real>> sideLineIntersection(const ArrayPoint<2, 2>& tetra, const ArrayPoint<2, 2>& line){
    // Constructs the relative coordinates of the intersections
    auto T0 = tetra.col(0);
    auto T1 = tetra.col(1);
    auto L0 = line.col(0);
    auto L1 = line.col(1);
    Eigen::Matrix<Real, 2, 2>  intersectionRelativeCoords;
    Eigen::Matrix<Real,2,2> mat;
    mat.col(0) = T1 - T0; mat.col(1) = L0 - L1;
    intersectionRelativeCoords.col(0) = mat.lu().solve(L0-T0);
    // If the intersections are not NaN
    if (!isnan(intersectionRelativeCoords(0,0))) {
        // There is a Real intersection only if the two relative coordinates are between 0 and 1, otherwise the intersection is beyond the segments
        if (((intersectionRelativeCoords.col(0).array() >= 0.) * (intersectionRelativeCoords.col(0).array() <= 1.)).all())
            return std::pair<std::vector<Real>, std::vector<Real>>({intersectionRelativeCoords(0,0)}, {intersectionRelativeCoords(1,0)});
        else
            return std::pair<std::vector<Real>, std::vector<Real>>({}, {});
    }
    
    // If we are here we know that the two lines are parallel
    // We first check if the two lines are aligned, by seeing if tetra[1]-tetra[0] and line[1]-tetra[0] are linearly independent
    mat.col(1) = T1 - L0;
    if ( mat.determinant() != 0. ){
        return std::pair<std::vector<Real>, std::vector<Real>>({}, {});
    }
    
    // Find the 4 intersections
    int nfound = 0;
    // tetra[0]
    Real inters = (T0-L0).dot(L1-L0) / (L1-L0).dot(L1-L0);
    if (inters >=0. and inters <= 1.) { intersectionRelativeCoords.col(nfound++) = Point<2>{0, inters};  }
    // tetra[1]
    inters = (T1-L0).dot(L1-L0) / (L1-L0).dot(L1-L0);
    if (inters >=0. and inters <= 1.) { intersectionRelativeCoords.col(nfound++) = Point<2>{1, inters};  }
    // line[0]
    inters = (L0-T0).dot(T1-T0) / (T1-T0).dot(T1-T0);
    if (inters >0. and inters < 1.) {intersectionRelativeCoords.col(nfound++) = Point<2>{inters, 0};}
    // line[1]
    inters = (L1-T0).dot(T1-T0) / (T1-T0).dot(T1-T0);
    if (inters >0. and inters < 1.) { intersectionRelativeCoords.col(nfound++) = Point<2>{inters, 1.};  }

    switch (nfound) {
        case 0:
            return std::pair<std::vector<Real>, std::vector<Real>>({}, {});
            break;
        case 1:
            return std::pair<std::vector<Real>, std::vector<Real>>({intersectionRelativeCoords(0,0)}, {intersectionRelativeCoords(1,0)});
            break;
        case 2:
            return std::pair<std::vector<Real>, std::vector<Real>>({intersectionRelativeCoords(0,0), intersectionRelativeCoords(0,1)}, {intersectionRelativeCoords(1,0), intersectionRelativeCoords(1,1)});
            break;
        default:
            break;
    }
    return std::pair<std::vector<Real>, std::vector<Real>>({}, {});
}
template<int NDim> std::tuple<bool, ArrayPoint<NDim, 2>, Eigen::Vector<Real, 2>> elementLineIntersection(const ArrayPoint<NDim, NDim+1>& elemenVerts, const ArrayPoint<NDim, 2>& line);
template<> std::tuple<bool, ArrayPoint<2, 2>, Eigen::Vector<Real, 2>> elementLineIntersection(const ArrayPoint<2, 2+1>& elemenVerts, const ArrayPoint<2, 2>& line){
    int                     numfound = 0;
    ArrayPoint<2, 2>        relativeCoordsElem;
    Eigen::Vector<Real, 2>  relativeCoordsLine;

    {
        auto [ intersSide, intersLine ] = sideLineIntersection<2>(elemenVerts.block<2, 2>(0, 0), line);
        switch (intersSide.size()) {
            case 1:
                relativeCoordsElem.col(numfound) = Point<2>{intersSide[0], 0.};
                relativeCoordsLine[numfound++] = intersLine[0];
                break;
            case 2:
                relativeCoordsElem.col(0) = Point<2>{intersSide[0], 0.};
                relativeCoordsLine[0] = intersLine[0];
                relativeCoordsElem.col(1) = Point<2>{intersSide[1], 0.};
                relativeCoordsLine[1] = intersLine[1];
                return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Vector<Real, 2>>{true, relativeCoordsElem, relativeCoordsLine};
                break;
        }
    }
    {
        auto [ intersSide, intersLine ] = sideLineIntersection<2>(elemenVerts.block<2, 2>(0, 1), line);
        switch (intersSide.size()) {
            case 1:
                if (numfound == 0){
                    relativeCoordsElem.col(numfound) = Point<2>{1-intersSide[0], intersSide[0]};
                    relativeCoordsLine[numfound++] = intersLine[0];
                } else if (relativeCoordsLine[0] != intersLine[0]) {
                    relativeCoordsElem.col(1) = Point<2>{1-intersSide[0], intersSide[0]};
                    relativeCoordsLine[1] = intersLine[0];
                    return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Vector<Real, 2>>{true, relativeCoordsElem, relativeCoordsLine};
                }
                break;
            case 2:
                relativeCoordsElem.col(0) = Point<2>{1-intersSide[0], intersSide[0]};
                relativeCoordsLine[0] = intersLine[0];
                relativeCoordsElem.col(1) = Point<2>{1-intersSide[1], intersSide[1]};
                relativeCoordsLine[1] = intersLine[1];
                return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Vector<Real, 2>>{true, relativeCoordsElem, relativeCoordsLine};
                break;
        }
    }
    {
        ArrayPoint<2, 2> side; side.col(0) = elemenVerts.col(0); side.col(1) = elemenVerts.col(2);
        auto [ intersSide, intersLine ] = sideLineIntersection<2>(side, line);
        switch (intersSide.size()) {
            case 1:
                if (numfound == 0){
                    relativeCoordsElem.col(numfound) = Point<2>{0, intersSide[0]};
                    relativeCoordsLine[numfound++] = intersLine[0];
                } else if (relativeCoordsLine[0] != intersLine[0]) {
                    relativeCoordsElem.col(1) = Point<2>{0, intersSide[0]};
                    relativeCoordsLine[1] = intersLine[0];
                    return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Vector<Real, 2>>{true, relativeCoordsElem, relativeCoordsLine};
                }
                break;
            case 2:
                relativeCoordsElem.col(0) = Point<2>{0, intersSide[0]};
                relativeCoordsLine[0] = intersLine[0];
                relativeCoordsElem.col(1) = Point<2>{0, intersSide[1]};
                relativeCoordsLine[1] = intersLine[1];
                return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Vector<Real, 2>>{true, relativeCoordsElem, relativeCoordsLine};
                break;
        }
    }
    // Now we look for points inside
    ArrayPoint<2, 2> sideVects; sideVects.col(0) = elemenVerts.col(1)- elemenVerts.col(0); sideVects.col(1) = elemenVerts.col(2)- elemenVerts.col(0);
    Point<2> relCoordsVLine = sideVects.lu().solve(line.col(0)-elemenVerts.col(0));
    if (relCoordsVLine(0) > 0. and relCoordsVLine(1) > 0. and relCoordsVLine(0) +relCoordsVLine(1) < 1.) {
        relativeCoordsElem.col(numfound) = relCoordsVLine;
        relativeCoordsLine[numfound++] = 0.;
        if (numfound == 2) {
            return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Vector<Real, 2>>{true, relativeCoordsElem, relativeCoordsLine};
        }
    }
    relCoordsVLine = sideVects.lu().solve(line.col(1)-elemenVerts.col(0));
    if (relCoordsVLine(0) > 0. and relCoordsVLine(1) > 0. and relCoordsVLine(0) +relCoordsVLine(1) < 1.) {
        relativeCoordsElem.col(numfound) = relCoordsVLine;
        relativeCoordsLine[numfound++] = 1.;
        if (numfound == 2) {
            return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Vector<Real, 2>>{true, relativeCoordsElem, relativeCoordsLine};
        }
    }
    return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Vector<Real, 2>>{false, relativeCoordsElem, relativeCoordsLine};
}

template<> void Function<2>::plot(const VectorPoint<2>& fullLine, const std::vector<std::string>& names, const std::string& t_title) const{
    Real cumulativeStart = 0.;
    Real minF = min(), maxF = max();
    for (int linePiece = 0; linePiece < fullLine.cols(); ++linePiece ){
        ArrayPoint<2, 2> line;
        line.col(0) = Point<2>{cumulativeStart, minF};
        line.col(1) = Point<2>{cumulativeStart, maxF};
        plotter2d << line << NOLABEL << COLOR("black");
        plotter2d << TEXT(names[linePiece], Point<2>{cumulativeStart, minF-0.02*(maxF-minF)});
        if ( linePiece < fullLine.cols()-1)  cumulativeStart += (fullLine.col(linePiece+1) - fullLine.col(linePiece)).norm();
    }
    plotter2d << OPENSTREAM;
    cumulativeStart = 0.;
    Real currentLength;
    for (int linePiece = 0; linePiece < fullLine.cols()-1; ++linePiece ){
        currentLength = (fullLine.col(linePiece+1) - fullLine.col(linePiece)).norm();
        ArrayPoint<2, 2> line = fullLine.block<2, 2>(0, linePiece);
        for (int elem=0; elem < mesh->numberElements; ++ elem){
            ArrayPoint<2, 3> elemenVerts;
            elemenVerts.col(0) = (*mesh)(elem, Point<2>{0,0});
            elemenVerts.col(1) = (*mesh)(elem, Point<2>{1,0});
            elemenVerts.col(2) = (*mesh)(elem, Point<2>{0,1});
            auto [found, relCoord, lineCoord] = elementLineIntersection(elemenVerts, line);
            if (found) {
                ArrayPoint<2, 2> segment;
                segment(0,0) = cumulativeStart + currentLength * lineCoord(0);
                segment(1,0) = (*this)(elem, relCoord.col(0));
                segment(0,1) = cumulativeStart + currentLength * lineCoord(1);
                segment(1,1) = (*this)(elem, relCoord.col(1));
                plotter2d << segment;
            }
        }
        cumulativeStart += currentLength;
    }
    plotter2d << CLOSESTREAM << PLOT;
};

Plotter3D& operator << (Plotter3D& plotter, const Function<2>& function){
    bool wasLineMode = (plotter.polygonModeActive == false);
    plotter << OPENSTREAM << POLYGONMODE;
    for (auto elem = function.elementIterator(); elem.unfinished; ++elem){
        ArrayPoint<3,2+1> triang; int i=0;
        for (auto vertex = function.vertexInElemIterator(); vertex.unfinished ; ++vertex){
            triang.block<2,1>(0,i) = (*function.mesh)(elem,vertex);
            triang(2,i++) = function(elem,vertex);
        }
        plotter << triang;
    }
    if (wasLineMode) plotter << LINEMODE;           // returns to the plot mode that was active before this method was called
    return plotter << CLOSESTREAM ;
}
Plotter2D& operator << (Plotter2D& plotter, const Function<1>& function){
    bool wasPolygonMode = (plotter.polygonModeActive == true);
    plotter << OPENSTREAM << LINEMODE;
    VectorPoint<2> line;
    line.resize(2, 2*function.mesh->numberElements);
    int i=0;
    for (auto elem = function.elementIterator(); elem.unfinished; ++elem){
        for (auto vertex = function.vertexInElemIterator(); vertex.unfinished ; ++vertex){
            line.block<1,1>(0,i) = (*function.mesh)(elem,vertex);
            line(1,i++) = function(elem,vertex);
        }
    }
    plotter << line;
    if (wasPolygonMode) plotter << POLYGONMODE;           // returns to the plot mode that was active before this method was called
    return plotter << CLOSESTREAM ;
}
Plotter2D& operator << (Plotter2D& plotter, const Projected<2>& projFunction){
    Real cumulativeStart = 0.;
    Real minF = projFunction.fun.min(), maxF = projFunction.fun.max();
    for (int linePiece = 0; linePiece < projFunction.line.cols(); ++linePiece ){
        ArrayPoint<2, 2> line;
        line.col(0) = Point<2>{cumulativeStart, minF};
        line.col(1) = Point<2>{cumulativeStart, maxF};
        plotter << line << NOLABEL << COLOR("black");
        plotter << TEXT(projFunction.names[linePiece], Point<2>{cumulativeStart, minF-0.04*(maxF-minF)});
        if ( linePiece < projFunction.line.cols()-1)  cumulativeStart += (projFunction.line.col(linePiece+1) - projFunction.line.col(linePiece)).norm();
    }
    plotter << OPENSTREAM;
    cumulativeStart = 0.;
    Real currentLength;
    for (int linePiece = 0; linePiece < projFunction.line.cols()-1; ++linePiece ){
        currentLength = (projFunction.line.col(linePiece+1) - projFunction.line.col(linePiece)).norm();
        ArrayPoint<2, 2> line = projFunction.line.block<2, 2>(0, linePiece);
        for (int elem=0; elem < projFunction.fun.mesh->numberElements; ++ elem){
            ArrayPoint<2, 3> elemenVerts;
            elemenVerts.col(0) = (*projFunction.fun.mesh)(elem, Point<2>{0,0});
            elemenVerts.col(1) = (*projFunction.fun.mesh)(elem, Point<2>{1,0});
            elemenVerts.col(2) = (*projFunction.fun.mesh)(elem, Point<2>{0,1});
            auto [found, relCoord, lineCoord] = elementLineIntersection(elemenVerts, line);
            if (found) {
                ArrayPoint<2, 2> segment;
                segment(0,0) = cumulativeStart + currentLength * lineCoord(0);
                segment(1,0) = projFunction.fun(elem, relCoord.col(0));
                segment(0,1) = cumulativeStart + currentLength * lineCoord(1);
                segment(1,1) = projFunction.fun(elem, relCoord.col(1));
                plotter << segment;
            }
        }
        cumulativeStart += currentLength;
    }
    return plotter << CLOSESTREAM;
};


// ===========================================
// ===========================================
// ===========================================
// Explicit template instantiation (For more info see: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl)

template class Function<1>;
template class Function<2>;
template class Function<3>;

template void swap(Function<1>& first, Function<1>& second);
template void swap(Function<2>& first, Function<2>& second);
template void swap(Function<3>& first, Function<3>& second);

template Function<1> operator/(const Real lhs, Function<1> rhs);
template Function<2> operator/(const Real lhs, Function<2> rhs);
template Function<3> operator/(const Real lhs, Function<3> rhs);

} // namespace Tortoise
