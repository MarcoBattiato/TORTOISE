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
#include <fstream>

using namespace Tortoise::GeometryCore;

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
    toNodalFromLinForm();
    vec = vec.unaryExpr(t_f);
    toLinFormFromNodal();
    
//    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), mesh->nVertElement, mesh->numberElements);
//    elementview.block(1,0,mesh->nVertElement-1,mesh->numberElements).rowwise() += elementview.row(0);
//    for (auto elem = mesh->elementIterator(); elem.unfinished ; ++elem) {
//        for (auto vert = mesh->vertexInElemIterator(); vert.unfinished; ++vert){
//            vec(vert.currentVertID + mesh->nVertElement * elem.currentElementID) = t_f(vec(vert.currentVertID + mesh->nVertElement * elem.currentElementID));
//        }
//    }
//    elementview.block(1,0,mesh->nVertElement-1,mesh->numberElements).rowwise() -= elementview.row(0);
};
template<int NDim> Function<NDim>::Function(const std::function<Real(std::vector<Real>)>& t_f, const std::vector<Function<NDim>>& t_others): mesh(t_others[0].mesh), vec(t_others[0].vec){
    for(auto it = ++t_others.begin(); it != t_others.end(); it++) { assert(t_others[0].mesh == it->mesh); }
    // The most convenient form is to transform the functions into a nodal representation
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

template<> Real Function<1>::integrateDiracDelta(const Function<1>& functInDiracDelta) const {
// Authors: Ching Haoye Trevor, Marco Battiato
    assert(mesh == functInDiracDelta.mesh);
    Real totalIntegr = 0;
    Point<1> relPoint0(0.), relPoint1(1.);
    for (int i = 0; i < (mesh->numberElements); ++i) { // on the condition that the intersection point is the edge of 2 elements, something needs to be written to prevent repetition.
        const Real  gValue0 = functInDiracDelta(i, relPoint0),  gValue1 = functInDiracDelta(i, relPoint1);
        if (!((gValue0 > 0. && gValue1 > 0.) || (gValue0 < 0. && gValue1 < 0.))) {
            if (gValue0 == 0. && gValue1 == 0.) {
                return std::numeric_limits<double>::quiet_NaN();
            } else {
                Point<1> relPt(gValue0 / (gValue0 - gValue1));
                Real integrNum = (*this)(i, relPt);
                if (gValue0 == 0. || gValue1 == 0.) { integrNum *= 0.5; } // halves the numerator so that the edge repetition does not double the integral value
                // For the denominator, we find the gradient:
                Real integrDenom = std::fabs(functInDiracDelta.derivative(0)(i, relPt));
                // we can now sum each element integral to get the total integral:
                totalIntegr += (integrNum / integrDenom);
            }
        }
    }
    return totalIntegr;
}
template<> Real Function<2>::integrateDiracDelta(const Function<2>& functInDiracDelta) const {
// Authors: Ching Haoye Trevor, Marco Battiato
    assert(mesh == functInDiracDelta.mesh);
    Real totalIntegr = 0.;
    bool lineBool = false;
    for (int i = 0; i < (mesh->numberElements); ++i) {
        const Real gValue00 = functInDiracDelta(i, { 0,0 });
        const Real gValue10 = functInDiracDelta(i, { 1,0 });
        const Real gValue01 = functInDiracDelta(i, { 0,1 });
        if (!((gValue00 > 0. && gValue10 > 0. && gValue01 > 0.) || (gValue00 < 0. && gValue10 < 0. && gValue01 < 0.))) {
            int j = 0;
            Eigen::Matrix2d ptsMatrix = Eigen::Matrix2d::Zero(); // Matrix to hold the valid coords for A and B
            // case 1: all vertices are 0 (integral diverging)
            if (gValue00 == 0. && gValue10 == 0. && gValue01 == 0.) {
                return std::numeric_limits<double>::quiet_NaN();
            } else if ((gValue00 == 0 && gValue10 == 0)) {  // case 2: 2 vertices are 0
                ptsMatrix(0, 0) = 0.; ptsMatrix(1, 0) = 0.; ptsMatrix(0, 1) = 1.; ptsMatrix(1, 1) = 0.;
                j = 2; lineBool = true;
            } else if ((gValue00 == 0 && gValue01 == 0)) {
                ptsMatrix(0, 0) = 0.; ptsMatrix(1, 0) = 0.; ptsMatrix(0, 1) = 0.; ptsMatrix(1, 1) = 1.;
                j = 2; lineBool = true;
            } else if ((gValue10 == 0 && gValue01 == 0)) {
                ptsMatrix(0, 0) = 1.; ptsMatrix(1, 0) = 0.; ptsMatrix(0, 1) = 0.; ptsMatrix(1, 1) = 1.;
                j = 2; lineBool = true;
            } else { // case 3: at most 1 vertex is 0
                if ((gValue00 > 0. && gValue10 <= 0.) || (gValue00 < 0. && gValue10 >= 0.)) { // only allows one vertex to be an intersection per side
                    ptsMatrix(0, 0) = (gValue00 / (gValue00 - gValue10)); ptsMatrix(1, 0) = 0.;
                    ++j;
                }
                if ((gValue00 >= 0. && gValue01 < 0.) || (gValue00 <= 0. && gValue01 > 0.)) {
                    ptsMatrix(0, j) = 0.; ptsMatrix(1, j) = (gValue00 / (gValue00 - gValue01));
                    ++j;
                }
                if ((gValue10 > 0. && gValue01 <= 0.) || (gValue10 < 0. && gValue01 >= 0.)) {
                    Real relFrom10 = (gValue10 / (gValue10 - gValue01)); // relative distance from 1,0 along the 1,0 - 0,1 side
                    ptsMatrix(0, 1) = (1 - relFrom10); ptsMatrix(1, 1) = relFrom10;
                    ++j;
                }
            }

            if (j > 1) { // Excludes the case where one vertex is 0, but there is no intersection with the opposite side
                // now we have the relative coordinates for points A and B
                Point<2> relPtA = (ptsMatrix.col(0)), relPtB = (ptsMatrix.col(1));
                Point<2> ptA((*mesh)(i, relPtA)), ptB((*mesh)(i, relPtB));

                Real integrNum;
                if (lineBool == 1) { // Finally, we can find the numerator of the integral for each element.
                    integrNum = 0.5 * (ptB - ptA).norm() * (*this)(i, { ((relPtA + relPtB) / 2)(0), ((relPtA + relPtB) / 2)(1) });
                } else {
                    integrNum = (ptB - ptA).norm() * (*this)(i, { ((relPtA + relPtB) / 2)(0), ((relPtA + relPtB) / 2)(1) });
                }
                // For the denominator, we make use of the gradient vector
                Eigen::Vector2d gradVec = { functInDiracDelta.derivative(0)(i, relPtA) , functInDiracDelta.derivative(1)(i, relPtA) };
                // the denom is then the norm of the grad vec
                Real integrDenom = gradVec.norm();

                // we can now sum each element integral to get the total integral:
                totalIntegr += (integrNum / integrDenom);
            }
        }
    }
    return totalIntegr;
}
template<> Real Function<3>::integrateDiracDelta(const Function<3>& functInDiracDelta) const {
// Authors: Ching Haoye Trevor, Marco Battiato
    assert(mesh == functInDiracDelta.mesh);
    Real totalIntegr = 0;
    bool surfaceBool = false;
    for (int i = 0; i < (mesh->numberElements); ++i) {
        const Real gValue000 = functInDiracDelta(i, { 0,0,0 }), gValue100 = functInDiracDelta(i, { 1,0,0 });
        const Real gValue010 = functInDiracDelta(i, { 0,1,0 }), gValue001 = functInDiracDelta(i, { 0,0,1 });

        if (!((gValue000 > 0. && gValue100 > 0. && gValue010 > 0. && gValue001 > 0.) || (gValue000 < 0. && gValue100 < 0. && gValue010 < 0. && gValue001 < 0.))) {
            int numZeros = 0;  // Number of points where the functInDiracDelta is 0.
            Eigen::Matrix<Real, 3, 4> ptsMatrix = Eigen::Matrix<Real, 3, 4>::Zero(); // Matrix to hold the valid coords for A and B
            if (gValue000 == 0. && gValue100 == 0. && gValue010 == 0. && gValue001 == 0.) {  // case 1: all vertices are 0 (integral diverging)
                return std::numeric_limits<double>::quiet_NaN();
            } else if (gValue000 == 0. && gValue100 == 0. && gValue010 == 0.) { // case 2: 3 vertices are 0
                ptsMatrix(0, 0) = 0.; ptsMatrix(1, 0) = 0.; ptsMatrix(2, 0) = 0.;
                ptsMatrix(0, 1) = 1.; ptsMatrix(1, 1) = 0.; ptsMatrix(2, 1) = 0.;
                ptsMatrix(0, 2) = 0.; ptsMatrix(1, 2) = 1.; ptsMatrix(2, 2) = 0.;
                numZeros = 3; surfaceBool = true;
            } else if (gValue000 == 0. && gValue010 == 0. && gValue001 == 0.) {
                ptsMatrix(0, 0) = 0.; ptsMatrix(1, 0) = 0.; ptsMatrix(2, 0) = 0.;
                ptsMatrix(0, 1) = 0.; ptsMatrix(1, 1) = 1.; ptsMatrix(2, 1) = 0.;
                ptsMatrix(0, 2) = 0.; ptsMatrix(1, 2) = 0.; ptsMatrix(2, 2) = 1.;
                numZeros = 3; surfaceBool = true;
            } else if (gValue000 == 0. && gValue100 == 0. && gValue001 == 0.) {
                ptsMatrix(0, 0) = 0.;ptsMatrix(1, 0) = 0.; ptsMatrix(2, 0) = 0.;
                ptsMatrix(0, 1) = 1.; ptsMatrix(1, 1) = 0.; ptsMatrix(2, 1) = 0.;
                ptsMatrix(0, 2) = 0.; ptsMatrix(1, 2) = 0.; ptsMatrix(2, 2) = 1.;
                numZeros = 3; surfaceBool = true;
            } else if (gValue100 == 0. && gValue010 == 0. && gValue001 == 0.) {
                ptsMatrix(0, 0) = 1.; ptsMatrix(1, 0) = 0.; ptsMatrix(2, 0) = 0.;
                ptsMatrix(0, 1) = 0.; ptsMatrix(1, 1) = 1.; ptsMatrix(2, 1) = 0.;
                ptsMatrix(0, 2) = 0.; ptsMatrix(1, 2) = 0.; ptsMatrix(2, 2) = 1.;
                numZeros = 3; surfaceBool = true;
            } else {  // case 3: at most 2 vertices are 0
                if ((gValue000 > 0. && gValue100 <= 0.) || (gValue000 < 0. && gValue100 >= 0.)) { // only allows one vertex to be an intersection per side
                    ptsMatrix(0, 0) = (gValue000 / (gValue000 - gValue100)); ptsMatrix(1, 0) = 0.; ptsMatrix(2, 0) = 0.;
                    ++numZeros; // potentially 100 vertex intsct
                }
                if ((gValue000 >= 0. && gValue010 < 0.) || (gValue000 <= 0. && gValue010 > 0.)) {
                    ptsMatrix(0, numZeros) = 0.; ptsMatrix(1, numZeros) = (gValue000 / (gValue000 - gValue010)); ptsMatrix(2, numZeros) = 0.;
                    ++numZeros; // potentially 000 vertex intsct
                }
                if ((gValue000 > 0. && gValue001 <= 0.) || (gValue000 < 0. && gValue001 >= 0.)) {
                    ptsMatrix(0, numZeros) = 0.; ptsMatrix(1, numZeros) = 0.; ptsMatrix(2, numZeros) = (gValue000 / (gValue000 - gValue001));
                    ++numZeros; // potentially 001 vertex intsct
                }
                if ((gValue100 > 0. && gValue010 <= 0.) || (gValue100 < 0. && gValue010 >= 0.)) {
                    Real relFrom100 = (gValue100 / (gValue100 - gValue010)); // relative distance from 1,0,0 along the 1,0,0 - 0,1,0 side
                    ptsMatrix(0, numZeros) = (1 - relFrom100); ptsMatrix(1, numZeros) = relFrom100; ptsMatrix(2, numZeros) = 0.;
                    ++numZeros; // potentially 010 vertex intsct
                }
                if ((gValue001 > 0. && gValue010 < 0.) || (gValue001 < 0. && gValue010 > 0.)) {
                    Real relFrom010 = (gValue010 / (gValue010 - gValue001)); // relative distance from 0,1,0 along the 0,0,1 - 0,1,0 side
                    ptsMatrix(0, numZeros) = 0.; ptsMatrix(1, numZeros) = (1 - relFrom010); ptsMatrix(2, numZeros) = relFrom010;
                    ++numZeros;
                }
                if ((gValue001 > 0. && gValue100 < 0.) || (gValue001 < 0. && gValue100 > 0.)) {
                    Real relFrom001 = (gValue001 / (gValue001 - gValue100)); // relative distance from 0,0,1 along the 0,0,1 - 1,0,0 side
                    ptsMatrix(0, numZeros) = relFrom001; ptsMatrix(1, numZeros) = 0.; ptsMatrix(2, numZeros) = (1 - relFrom001);
                    ++numZeros;
                }

                if (gValue000 == 0. && gValue100 == 0.) // 000 is recorded
                {
                    ptsMatrix(0, numZeros) = 1.;
                    ptsMatrix(1, numZeros) = 0.;
                    ptsMatrix(2, numZeros) = 0.;
                    ++numZeros;
                }
                if (gValue000 == 0. && gValue010 == 0.) // 010 is recorded
                {
                    ptsMatrix(0, numZeros) = 0.;
                    ptsMatrix(1, numZeros) = 0.;
                    ptsMatrix(2, numZeros) = 0.;
                    ++numZeros;
                }
                if (gValue010 == 0. && gValue100 == 0.) // 100 is recorded
                {
                    ptsMatrix(0, numZeros) = 0.;
                    ptsMatrix(1, numZeros) = 1.;
                    ptsMatrix(2, numZeros) = 0.;
                    ++numZeros;
                }
                if (gValue000 == 0. && gValue001 == 0.) // 000 is recorded
                {
                    ptsMatrix(0, numZeros) = 0.;
                    ptsMatrix(1, numZeros) = 0.;
                    ptsMatrix(2, numZeros) = 1.;
                    ++numZeros;
                }
            }
            if (numZeros == 3) // Excludes the case where 2 vertices are 0, but there is no intersection with the opposite side
            {
                // assignment to pts and integration for 3 pt intsction
                Point<3> relPtA = (ptsMatrix.col(0));
                Point<3> relPtB = (ptsMatrix.col(1));
                Point<3> relPtC = (ptsMatrix.col(2));

                Point<3> ptA((*mesh)(i, relPtA));
                Point<3> ptB((*mesh)(i, relPtB));
                Point<3> ptC((*mesh)(i, relPtC));

                Eigen::Vector3d      gradVec = { functInDiracDelta.derivative(0)(i, relPtA) , functInDiracDelta.derivative(1)(i, relPtA) , functInDiracDelta.derivative(2)(i, relPtA) };
                Real integrDenom = gradVec.norm();
                Real integrNum;

                if (surfaceBool == 1)
                {
                    // the base area for the volume integral (numerator) is half the cross product norm of AB and AC
                    integrNum = 0.25 * ((ptB - ptA).cross((ptC - ptA))).norm() *
                                ((*this)(i, { ((relPtA + relPtB + relPtC) / 3)(0), ((relPtA + relPtB + relPtC) / 3)(1), ((relPtA + relPtB + relPtC) / 3)(2) }));
                }

                else
                {
                    integrNum = 0.5 * ((ptB - ptA).cross((ptC - ptA))).norm() *
                                ((*this)(i, { ((relPtA + relPtB + relPtC) / 3)(0), ((relPtA + relPtB + relPtC) / 3)(1), ((relPtA + relPtB + relPtC) / 3)(2) }));
                }

                totalIntegr += (integrNum / integrDenom);
            }

            if (numZeros == 4) // special case, 4 pt intsction
            {
                // assignment to pts and integration for 4 pt intsction
                Point<3> relPtA = (ptsMatrix.col(0));
                Point<3> relPtB = (ptsMatrix.col(1));
                Point<3> relPtC = (ptsMatrix.col(2));
                Point<3> relPtD = (ptsMatrix.col(3));

                Point<3> ptA((*mesh)(i, relPtA));
                Point<3> ptB((*mesh)(i, relPtB));
                Point<3> ptC((*mesh)(i, relPtC));
                Point<3> ptD((*mesh)(i, relPtD));

                // the base area for the volume integral (numerator) is the sum of 2 triangles (where the quadrilateral is split),
                // but the order of the 4 points can be random. easiest way is to find the midline splitting the quadrilateral, by probing the angles from a fixed point
                // fixed pt chosen to be A.
                Point<3> vecAB(ptB - ptA);
                Point<3> vecAC(ptC - ptA);
                Point<3> vecAD(ptD - ptA);

                Eigen::Vector3d  gradVec = { functInDiracDelta.derivative(0)(i, relPtA) , functInDiracDelta.derivative(1)(i, relPtA) , functInDiracDelta.derivative(2)(i, relPtA) };
                Real integrDenom = gradVec.norm();

                if (vecAB.norm() < 0.000001)
                {
                    Real integrNum = 0.5 * ((vecAC).cross((vecAD))).norm() *
                        ((*this)(i, { ((relPtA + relPtC + relPtD) / 3)(0), ((relPtA + relPtC + relPtD) / 3)(1), ((relPtA + relPtC + relPtD) / 3)(2) }));
                    totalIntegr += (integrNum / integrDenom);
                }
                else if ((vecAC.norm() < 0.000001) || ((ptC - ptB).norm() < 0.000001))
                {
                    Real integrNum = 0.5 * ((vecAB).cross((vecAD))).norm() *
                        ((*this)(i, { ((relPtA + relPtB + relPtD) / 3)(0), ((relPtA + relPtB + relPtD) / 3)(1), ((relPtA + relPtB + relPtD) / 3)(2) }));
                    totalIntegr += (integrNum / integrDenom);
                }
                else if ((vecAD.norm() < 0.000001) || ((ptD - ptB).norm() < 0.000001) || ((ptD - ptC).norm() < 0.000001))
                {
                    Real integrNum = 0.5 * ((vecAB).cross((vecAC))).norm() *
                        ((*this)(i, { ((relPtA + relPtB + relPtC) / 3)(0), ((relPtA + relPtB + relPtC) / 3)(1), ((relPtA + relPtB + relPtC) / 3)(2) }));
                    totalIntegr += (integrNum / integrDenom);
                }
                else
                {
                    Real angleABAC(std::acos((vecAB.dot(vecAC)) / (vecAB.norm() * vecAC.norm())));
                    Real angleACAD(std::acos((vecAC.dot(vecAD)) / (vecAC.norm() * vecAD.norm())));
                    Real angleABAD(std::acos((vecAB.dot(vecAD)) / (vecAB.norm() * vecAD.norm())));

                    if ((angleABAC >= angleACAD) && (angleABAC >= angleABAD)) // AD is the midline
                    {
                        Real integrNum = 0.5 * ((vecAB.cross(vecAD)).norm() + (vecAC.cross(vecAD)).norm()) *
                            ((*this)(i, { ((relPtA + relPtB + relPtC + relPtD) / 4)(0), ((relPtA + relPtB + relPtC + relPtD) / 4)(1), ((relPtA + relPtB + relPtC + relPtD) / 4)(2) }));
                        totalIntegr += (integrNum / integrDenom);
                    }
                    else if ((angleACAD >= angleABAC) && (angleACAD >= angleABAD)) // AB is the midline
                    {
                        Real integrNum = 0.5 * ((vecAB.cross(vecAC)).norm() + (vecAB.cross(vecAD)).norm()) *
                            ((*this)(i, { ((relPtA + relPtB + relPtC + relPtD) / 4)(0), ((relPtA + relPtB + relPtC + relPtD) / 4)(1), ((relPtA + relPtB + relPtC + relPtD) / 4)(2) }));
                        totalIntegr += (integrNum / integrDenom);
                    }
                    else // AC is the midline
                    {
                        Real integrNum = 0.5 * ((vecAB.cross(vecAC)).norm() + (vecAC.cross(vecAD)).norm()) *
                            ((*this)(i, { ((relPtA + relPtB + relPtC + relPtD) / 4)(0), ((relPtA + relPtB + relPtC + relPtD) / 4)(1), ((relPtA + relPtB + relPtC + relPtD) / 4)(2) }));
                        totalIntegr += (integrNum / integrDenom);
                    }
                }
            }

        }
    }

    return totalIntegr;
}

template<int NDim> Mesh<NDim> Function<NDim>::convolutionMesh() const {
    return mesh->convolutionMesh();
}

template<> Function<1> Function<1>::convolve(const Function<1>& functionG) const {
    const Mesh<1> convolFMesh(convolutionMesh());
    assert((*functionG.mesh) == convolFMesh);

    Point<1> relPoint0(0.);
    Point<1> relPoint1(1.);
    Eigen::Matrix<Real, Eigen::Dynamic, 2> integrStore;
    integrStore.resize((mesh->numberElements) + 1, 2);
    DataVector intStore;
    intStore.resize(1, 2 * (mesh->numberElements));

    for (int P = 0; P < (mesh->numberElements)+1; ++P)
    {
        // for mesh points ( of which there are numberElements+1)

        Real fMeshNum(mesh->numberElements);
        Point<1> fOrigin((mesh->relativeOrigin)); // gives starting boundary of functionF
        Point<1> fVecDir((mesh->region->gVec)*(mesh->relativegVec)); // gives the overall direction of the mesh vector for f

        //Point<1> ptPCoord(fOrigin + (P / fMeshNum) * (fGVec));
        //Point<1> newOrigin;
       // Point<1> newGVec;
        int functFStartElem;
        //int gVecPositive;

        if (fVecDir(0) > 0)
        {
            //newOrigin = (ptPCoord - fOrigin - fGVec);
            //newGVec = (fGVec);
            functFStartElem = fMeshNum - 1;
            //functFEndElem = 0;
            //gVecPositive = 1;
        }
        else
        {
            //newOrigin = (ptPCoord-fOrigin);
            //newGVec = (-fGVec);
            functFStartElem = 0;
            //functFEndElem = (functionF.mesh->numberElements)-1;
            //gVecPositive = 0;
        }

        Real elemWidth(std::fabs(mesh->dGVec(0)));
        int functGStartElem(P);
        
        // Point P should now be at 0, and the new functionF maps left to right
        // map and match the elements for functionG (functionG always maps left to right)
        
        // Now, to find the linear equation for each element
        Real pIntegr = 0;
        for (int i = 0; i < fMeshNum; ++i) // working left to right
        {

            // modified gauss quadrature
            Real elemIntegr = 0;
            for (int j = 0; j < ngausspoint(1, 2); ++j)
            {
                Point<1> gaussOrd2Pt(gaussPoints<1, 2>(0, j)), mirrorGaussOrd2Pt(1. - gaussPoints<1, 2>(0, j));

                if (fVecDir(0) > 0)
                {
                    elemIntegr += elemWidth * gaussWeigths<1, 2>(0, 1) * functionG(functGStartElem, gaussOrd2Pt) * (*this)(functFStartElem, mirrorGaussOrd2Pt);
                }
                
                else
                {
                    elemIntegr += elemWidth * gaussWeigths<1, 2>(0, 1) * functionG(functGStartElem, gaussOrd2Pt) * (*this)(functFStartElem, gaussOrd2Pt);
                }
                
            } // gives integral for each element, iterating thru gauss points
            // cout << "integral for elem " << i << " is: " << elemIntegr << endl;

            functGStartElem += 1;
            if (fVecDir(0) > 0)
            {
                functFStartElem -= 1;
            }
            else
            {
                functFStartElem += 1;
            }
            pIntegr += elemIntegr;
        }
        
        integrStore(P, 0) = P;
        integrStore(P, 1) = pIntegr;
        
        if (fVecDir(0) > 0) {
            if (P == 0) {
                intStore(0, 0) = pIntegr;
            } else if (P == fMeshNum) {
                intStore(0, 2 * (mesh->numberElements) - 1) = pIntegr;
            } else {
                intStore(0, 2 * P - 1) = pIntegr; intStore(0, 2 * P) = pIntegr;
            }
        }
        
        if (fVecDir(0) < 0) {
            if (P == 0) {
                intStore(0, 2 * (mesh->numberElements) - 1) = pIntegr;
            } else if (P == fMeshNum) {
                intStore(0, 0) = pIntegr;
            } else {
                intStore(0, 2*(mesh->numberElements) - (2 * P)) = pIntegr;
                intStore(0, 2*(mesh->numberElements) - (2 * P + 1)) = pIntegr;
            }
        }
        
        // translate left by P and reflect all points about P
        // or reflect all points about origin then shift right by P
    }

    return Function<1>(mesh, intStore);
}
template<> Function<2> Function<2>::convolve(const Function<2>& functionG) const{
    assert(false && "Function<2>::convolve not implemented");
    return functionG;
}
template<> Function<3> Function<3>::convolve(const Function<3>& functionG) const{
    assert(false && "Function<3>::convolve not implemented");
    return functionG;
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
template<int NDim> std::tuple<bool, ArrayPoint<NDim, 2>, Eigen::Matrix<Real, 2,1>> elementLineIntersection(const ArrayPoint<NDim, NDim+1>& elemenVerts, const ArrayPoint<NDim, 2>& line);
template<> std::tuple<bool, ArrayPoint<2, 2>, Eigen::Matrix<Real, 2,1>> elementLineIntersection(const ArrayPoint<2, 2+1>& elemenVerts, const ArrayPoint<2, 2>& line){
    int                     numfound = 0;
    ArrayPoint<2, 2>        relativeCoordsElem;
    Eigen::Matrix<Real, 2,1>  relativeCoordsLine;

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
                return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Matrix<Real, 2,1>>{true, relativeCoordsElem, relativeCoordsLine};
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
                    return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Matrix<Real, 2,1>>{true, relativeCoordsElem, relativeCoordsLine};
                }
                break;
            case 2:
                relativeCoordsElem.col(0) = Point<2>{1-intersSide[0], intersSide[0]};
                relativeCoordsLine[0] = intersLine[0];
                relativeCoordsElem.col(1) = Point<2>{1-intersSide[1], intersSide[1]};
                relativeCoordsLine[1] = intersLine[1];
                return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Matrix<Real, 2,1>>{true, relativeCoordsElem, relativeCoordsLine};
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
                    return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Matrix<Real, 2,1>>{true, relativeCoordsElem, relativeCoordsLine};
                }
                break;
            case 2:
                relativeCoordsElem.col(0) = Point<2>{0, intersSide[0]};
                relativeCoordsLine[0] = intersLine[0];
                relativeCoordsElem.col(1) = Point<2>{0, intersSide[1]};
                relativeCoordsLine[1] = intersLine[1];
                return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Matrix<Real, 2,1>>{true, relativeCoordsElem, relativeCoordsLine};
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
            return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Matrix<Real, 2,1>>{true, relativeCoordsElem, relativeCoordsLine};
        }
    }
    relCoordsVLine = sideVects.lu().solve(line.col(1)-elemenVerts.col(0));
    if (relCoordsVLine(0) > 0. and relCoordsVLine(1) > 0. and relCoordsVLine(0) +relCoordsVLine(1) < 1.) {
        relativeCoordsElem.col(numfound) = relCoordsVLine;
        relativeCoordsLine[numfound++] = 1.;
        if (numfound == 2) {
            return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Matrix<Real, 2,1>>{true, relativeCoordsElem, relativeCoordsLine};
        }
    }
    return std::tuple<bool, ArrayPoint<2, 2>, Eigen::Matrix<Real, 2,1>>{false, relativeCoordsElem, relativeCoordsLine};
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

template<> void Function<1>::writeToTxtFilePlot(const std::string &FileName) const {
    std::ofstream myfile (FileName);
    for (auto elem = elementIterator(); elem.unfinished; ++elem){
        for (auto vert = vertexInElemIterator(); vert.unfinished; ++vert){
            myfile << (*mesh)(elem,vert) << " " << (*this)(elem,vert) << "\n";
        }
    }
    myfile.close();
}

//Added Dec 2023 by Xavier for Function<2>
template<> void Function<2>::writeToTxtFilePlotPy(std::string &FileNameP)  const  {
    /*std::cout << "\nwriteToTxtFilePlotPy called\n";
    std::filesystem::path currentDir = std::filesystem::current_path();
    std::cout << "Current directory: " << currentDir << std::endl;
    std::cout << FileNameP;*/
    std::ofstream myfile(FileNameP);
    for (auto elem = elementIterator(); elem.unfinished; ++elem) {
        for (auto vert = vertexInElemIterator(); vert.unfinished; ++vert) {
            myfile << (*mesh)(elem, vert).transpose() << " " << (*this)(elem, vert) << "\n";
        }
    }
    myfile.close();
};

//Added Dec 2023 by Xavier for Function<2>
template<> void Function<2>::writeToTxtFilePlotAppend(std::string &FileName) const {
    //std::cout << "\nwriteToTxtFilePlotAppend called\n";
    //std::cout << FileName;
    std::ifstream myfile(FileName); // , std::ios::app
    std::ofstream tmpfile("./TORTOISE/TORTOISE/Plotter3D/Temp_Data_Files/Temp1.txt"); // Stores last column of new Function
    std::ofstream tmpfile2("./TORTOISE/TORTOISE/Plotter3D/Temp_Data_Files/Temp2.txt"); // Stores the combination of all Functions
    std::string line;
    std::string line2;
    for (auto elem = elementIterator(); elem.unfinished; ++elem) {
        for (auto vert = vertexInElemIterator(); vert.unfinished; ++vert) {
            tmpfile << (*this)(elem, vert) << std::endl;
        };
    };
    tmpfile.close();
    /*std::cout << "DEBUG HERE";
    std::filesystem::path currentDir = std::filesystem::current_path();
    std::cout << "Current directory: " << currentDir << std::endl;*/
    std::ifstream tmpfile3("./TORTOISE/TORTOISE/Plotter3D/Temp_Data_Files/Temp1.txt"); // Converts Temp1.txt to a readable file for std::getline


    while (std::getline(myfile, line) && std::getline(tmpfile3, line2)) {
        /*line += " " + (*this)(elem, vert);*/ // Append the last data to the line
        //line += " " + num;
        //std::cout << line << line2;


        tmpfile2 << line << " " << line2 << std::endl; // Write modified line to temporary file
    }

    myfile.close();
    tmpfile2.close();

    std::remove(FileName.c_str()); // Remove the existing file if it exists
    //std::cout << FileName.c_str();

    if (std::rename("./TORTOISE/TORTOISE/Plotter3D/Temp_Data_Files/Temp2.txt", FileName.c_str()) != 0) {
        std::cerr << "Error renaming file." << std::endl;
    }

    myfile.close();
    tmpfile.close();
    tmpfile2.close();
    tmpfile3.close();


};

//Added Dec 2023 by Xavier for Function<2>
template<> void Function<2>::writeToTxtFilePlotColour(std::string& FileNameP) const {
    /*std::cout << "Colour";
    std::filesystem::path currentDir = std::filesystem::current_path();
    std::cout << "Current directory: " << currentDir << std::endl;
    std::cout << FileNameP;*/
    std::ofstream myfile(FileNameP, std::ios::app);
    // Only appends 3rd var of Function Objects to a file
    for (auto elem = elementIterator(); elem.unfinished; ++elem) {
        for (auto vert = vertexInElemIterator(); vert.unfinished; ++vert) {
            myfile << (*this)(elem, vert) << "\n";
        }
    }
    myfile.close();
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
Plotter2D& operator << (Plotter2D& plotter, const ProjectedShifted<2>& projFunction){
    Real cumulativeStart = projFunction.shift;
    Real minF = projFunction.fun.min(), maxF = projFunction.fun.max();
    Real totalLength = 0;
    for (int linePiece = 0; linePiece < projFunction.line.cols() - 1; ++linePiece ){
        totalLength += (projFunction.line.col(linePiece+1) - projFunction.line.col(linePiece)).norm();
    }
        
    for (int linePiece = 0; linePiece < projFunction.line.cols(); ++linePiece ){
        ArrayPoint<2, 2> line;
        line.col(0) = Point<2>{cumulativeStart, minF};
        line.col(1) = Point<2>{cumulativeStart, maxF};
        plotter << line << NOLABEL << COLOR("black");
        if ( linePiece == 0 )
            plotter << TEXT(projFunction.names[linePiece], Point<2>{cumulativeStart + 0.01*totalLength, minF-0.04*(maxF-minF)});
        else if ( linePiece == projFunction.line.cols() -1 )
            plotter << TEXT(projFunction.names[linePiece], Point<2>{cumulativeStart - 0.02*totalLength, minF-0.04*(maxF-minF)});
        else
            plotter << TEXT(projFunction.names[linePiece], Point<2>{cumulativeStart - 0.05*totalLength, minF-0.04*(maxF-minF)});
        if ( linePiece < projFunction.line.cols()-1)  cumulativeStart += (projFunction.line.col(linePiece+1) - projFunction.line.col(linePiece)).norm();
    }
    plotter << OPENSTREAM;
    cumulativeStart = projFunction.shift;
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
