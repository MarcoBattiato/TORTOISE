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
//  FunctionElement.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//

#include <Geometry/Structured/FunctionSpace/FunctionElement.hpp>

#include <cassert>

namespace Tortoise {

// Output for EigenMatrix types
inline std::string to_string(const Eigen::MatrixXd& mat){
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

//=======================================================
// Constructors
//===================
//template<int NDim> FunctionElement<NDim>::FunctionElement(const Mesh<NDim>& t_mesh, const ElementID<NDim> elID) : mesh(&t_mesh), elemID(elID), vec(LinearForm<NDim>::Zero()) {}
//template<int NDim> FunctionElement<NDim>::FunctionElement(const Mesh<NDim>& t_mesh, const ElementID<NDim> elID, const Real t_value): mesh(&t_mesh), elemID(elID), vec(LinearForm<NDim>::Zero()) {vec(0)=t_value;}
//template<int NDim> template<typename Derived> FunctionElement<NDim>::FunctionElement(const Mesh<NDim>& t_mesh, const ElementID<NDim> elID, const Eigen::MatrixBase<Derived>& functRepr): mesh(&t_mesh), elemID(elID), vec(functRepr) {}
template<int NDim> FunctionElement<NDim>::FunctionElement(const FunctionElement<NDim>& other): mesh(other.mesh), elemID(other.elemID), vec(other.vec)  {}
template<int NDim> FunctionElement<NDim>::FunctionElement(FunctionElement<NDim>&& other): mesh(other.mesh), elemID(other.elemID), vec(std::move(other.vec))  {}

template<int NDim> FunctionElement<NDim>::FunctionElement(const Mesh<NDim>& t_mesh, const MeshElementIterator<NDim>& elID) : mesh(&t_mesh), elemID(elID.currentElementID), vec(LinearForm<NDim>::Zero()) {}
template<int NDim> FunctionElement<NDim>::FunctionElement(const Mesh<NDim>& t_mesh, const MeshElementIterator<NDim>& elID, const Real t_value): mesh(&t_mesh), elemID(elID.currentElementID), vec(LinearForm<NDim>::Zero()) {vec(0)=t_value;}
template<int NDim> FunctionElement<NDim>::FunctionElement(const Mesh<NDim>& t_mesh, const int elID) : mesh(&t_mesh), elemID(elID), vec(LinearForm<NDim>::Zero()) {}
template<int NDim> FunctionElement<NDim>::FunctionElement(const Mesh<NDim>& t_mesh, const int elID, const Real t_value): mesh(&t_mesh), elemID(elID), vec(LinearForm<NDim>::Zero()) {vec(0)=t_value;}
template<int NDim> FunctionElement<NDim>::FunctionElement(const Mesh<NDim>* t_mesh, const MeshElementIterator<NDim>& elID) : mesh(t_mesh), elemID(elID.currentElementID), vec(LinearForm<NDim>::Zero()) {}
template<int NDim> FunctionElement<NDim>::FunctionElement(const Mesh<NDim>* t_mesh, const MeshElementIterator<NDim>& elID, const Real t_value): mesh(t_mesh), elemID(elID.currentElementID), vec(LinearForm<NDim>::Zero()) {vec(0)=t_value;}
template<int NDim> FunctionElement<NDim>::FunctionElement(const Mesh<NDim>* t_mesh, const int elID) : mesh(t_mesh), elemID(elID), vec(LinearForm<NDim>::Zero()) {}
template<int NDim> FunctionElement<NDim>::FunctionElement(const Mesh<NDim>* t_mesh, const int elID, const Real t_value): mesh(t_mesh), elemID(elID), vec(LinearForm<NDim>::Zero()) {vec(0)=t_value;}



//********************************
//* Arithmetic
//********************************
// ATTENTION!!!!!!! Multiplications and divisions by anything except a scalar do not conserve precision!
// If you don't know what this means, do not use those operations! You will get unexpected results in most of the cases

template<int NDim> FunctionElement<NDim>& FunctionElement<NDim>::operator=(FunctionElement<NDim> other){
    assert(mesh == other.mesh && elemID == other.elemID);
    vec.swap(other.vec);
    return *this;
}
template<int NDim> FunctionElement<NDim>& FunctionElement<NDim>::operator=(FunctionElement<NDim>&& other){
    assert(mesh == other.mesh && elemID == other.elemID);
    vec.swap(other.vec);
    return *this;
}

template<int NDim> FunctionElement<NDim>& FunctionElement<NDim>::operator+=(const FunctionElement<NDim>& other){
    assert(mesh == other.mesh && elemID == other.elemID);
    vec += other.vec;
    return *this;
}
template<int NDim> FunctionElement<NDim>& FunctionElement<NDim>::operator-=(const FunctionElement<NDim>& other) {
    assert(mesh == other.mesh && elemID == other.elemID);
    vec -= other.vec;
    return *this;
}
template<int NDim> FunctionElement<NDim>& FunctionElement<NDim>::operator*=(const FunctionElement<NDim>& other){
    assert(mesh == other.mesh && elemID == other.elemID);
    // The product is calculated to make sure that the final values at the nodes are the product of the two functions
    // function 1 : values at nodes {f1.vec(0), f1.vec(0)+f1.vec(1), ... }
    // function 2 : values at nodes {f2.vec(0), f2.vec(0)+f2.vec(1), ... }
    // product    : values at nodes {f1.vec(0)*f2.vec(0), (f1.vec(0)+f1.vec(1))*(f2.vec(0)+f2.vec(1)), ... }
    // product in linear form       {f1.vec(0)*f2.vec(0), (f1.vec(0)+f1.vec(1))*(f2.vec(0)+f2.vec(1)) -  f1.vec(0)*f2.vec(0), ... }
    vec.template block<1,NDim>(0,1).array() += vec(0);
    vec(0) *= other.vec(0);
    vec.template block<1,NDim>(0,1) = vec.template block<1,NDim>(0,1).array()*( other.vec.template block<1,NDim>(0,1).array() + other.vec(0)) - vec(0);
    return *this;
}
template<int NDim> FunctionElement<NDim>& FunctionElement<NDim>::operator/=(const FunctionElement<NDim>& other) {
    assert(mesh == other.mesh && elemID == other.elemID);
    vec.template block<1,NDim>(0,1).array() += vec(0);
    vec(0) /= other.vec(0);
    vec.template block<1,NDim>(0,1) = vec.template block<1,NDim>(0,1).array() / ( other.vec.template block<1,NDim>(0,1).array() + other.vec(0)) - vec(0);
    return *this;
}

template<int NDim> FunctionElement<NDim> FunctionElement<NDim>::operator-() const {
    FunctionElement<NDim> toreturn (*this);
    toreturn.vec = -toreturn.vec;
    return toreturn;
}

template<int NDim> FunctionElement<NDim>& FunctionElement<NDim>::operator=(const Real& other){
    vec(0) = other;
    vec.template block<1,NDim>(0,1) = Eigen::Matrix<Real, 1, NDim>::Zero();
    return *this;
}

template<int NDim> FunctionElement<NDim>& FunctionElement<NDim>::operator+=(const Real& other){
    vec(0) += other;
    return *this;
}
template<int NDim> FunctionElement<NDim>& FunctionElement<NDim>::operator-=(const Real& other){
    vec(0) -= other;
    return *this;
}
template<int NDim> FunctionElement<NDim>& FunctionElement<NDim>::operator*=(const Real& other) {
    vec *= other;
    return *this;
}
template<int NDim> FunctionElement<NDim>& FunctionElement<NDim>::operator/=(const Real& other) {
    vec /= other;
    return *this;
}
template<int NDim> FunctionElement<NDim> operator/(const Real& lhs, FunctionElement<NDim> rhs) {
    rhs.vec.template block<1,NDim>(0,1).array() += rhs.vec(0);
    rhs.vec(0) = lhs/rhs.vec(0);
    rhs.vec.template block<1,NDim>(0,1).array() = lhs/rhs.vec.template block<1,NDim>(0,1).array() -rhs.vec(0);
    return rhs;
}

//********************************
//* Local Evaluation
//********************************
template<int NDim> Real FunctionElement<NDim>::operator()(const Point<NDim>& t_point_ref) const{
    return apply(vec,t_point_ref).sum();
}
template<int NDim> template <typename Derived> auto FunctionElement<NDim>::operator() (const Eigen::MatrixBase<Derived>& points_ref) const {
    return apply(vec,points_ref);
}


//********************************
//* Other operations
//********************************

template<int NDim>  Real FunctionElement<NDim>::max() const {
    static const Eigen::Matrix<Real, NDim+1, NDim+1> representconversion = [] {
        Eigen::Matrix<Real, NDim+1, NDim+1> tmp;
        tmp(0,0) = 1.0;
        for (int i=0; i<NDim; ++i) {tmp(0,i+1) = 1.0;tmp(i+1,i+1) = 1.0;}
        return tmp; }();
    return (vec*representconversion).maxCoeff();
}
template<int NDim>  Real FunctionElement<NDim>::min() const {
    static const Eigen::Matrix<Real, NDim+1, NDim+1> representconversion = [] {
        Eigen::Matrix<Real, NDim+1, NDim+1> tmp;
        tmp(0,0) = 1.0;
        for (int i=0; i<NDim; ++i) {tmp(0,i+1) = 1.0;tmp(i+1,i+1) = 1.0;}
        return tmp; }();
    return (vec*representconversion).minCoeff();
}
template<int NDim> Real FunctionElement<NDim>::average() const {
    return (vec.template block<1,NDim>(0,1).sum()+vec(0)*(static_cast<Real>(NDim+1)))/static_cast<Real>(NDim+1);
}
template<int NDim> Real FunctionElement<NDim>::integrate() const {
    return (apply(vec, gaussPoints<NDim,1>).array()*gaussWeigths<NDim,1>.array()).sum()*mesh->elemVolume;
}
template<int NDim> Real FunctionElement<NDim>::integrate(const FunctionElement<NDim>& rhs) const {
    return (apply(vec, gaussPoints<NDim,2>).array()*apply(rhs.vec, gaussPoints<NDim,2>).array()*gaussWeigths<NDim,2>.array()).sum()*mesh->elemVolume;
}
template<int NDim> Real FunctionElement<NDim>::integrate(const FunctionElement<NDim>& rhs1, const FunctionElement<NDim>& rhs2) const {
    return (apply(vec, gaussPoints<NDim,3>).array()*apply(rhs1.vec, gaussPoints<NDim,3>).array()*apply(rhs2.vec, gaussPoints<NDim,3>).array()*gaussWeigths<NDim,3>.array()).sum()*mesh->elemVolume;
}


//********************************
//* I/O
//********************************

template<int NDim> std::ostream &operator<<(std::ostream &t_os, const FunctionElement<NDim>& funct){
    std::string output;
    static const Eigen::Matrix<Real, NDim+1, NDim+1> representconversion = [] {
    Eigen::Matrix<Real, NDim+1, NDim+1> tmp;
    tmp(0,0) = 1.0;
    for (int i=0; i<NDim; ++i) {tmp(0,i+1) = 1.0;tmp(i+1,i+1) = 1.0;}
    return tmp; }();
    output+= "********************************************\n";
    output+= "elemID  : " + std::to_string(funct.elemID) + "\n";
    output+= "linFun  : " + to_string(funct.vec) + "\n";
    output+= "nodeVal : " + to_string(funct.vec*representconversion) + "\n";
    output+= "********************************************\n";
    return t_os << output;
}

// Explicit template instantiation (For more info see: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl)

template class FunctionElement<1>;
template class FunctionElement<2>;
template class FunctionElement<3>;

template  std::ostream &operator<<(std::ostream &os, const FunctionElement<1> & funct);
template  std::ostream &operator<<(std::ostream &os, const FunctionElement<2> & funct);
template  std::ostream &operator<<(std::ostream &os, const FunctionElement<3> & funct);


} // namespace Tortoise
