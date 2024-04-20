//
//  PolygonalFunction.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 4/11/22.
//

#include <Geometry/Structured/FunctionSpace/Polygonal/PolygonalFunction.hpp>

#include <cassert>

namespace Tortoise {

//=======================================================
// Constructors
//===================
template<int NDim> PolygonalFunction<NDim>::PolygonalFunction(const PolygonalMesh<NDim>& t_mesh): mesh(&t_mesh), vec(DataVector::Zero(t_mesh.dataVectorDim)){};

template<int NDim> PolygonalFunction<NDim>::PolygonalFunction(const PolygonalMesh<NDim>& t_mesh, const Real t_value): mesh(&t_mesh), vec(DataVector::Zero(t_mesh.dataVectorDim)){
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), mesh->nVertElement, mesh->numberElements);
    elementview.row(0).array() = t_value;
};
template<int NDim> PolygonalFunction<NDim>::PolygonalFunction(const PolygonalMesh<NDim>& t_mesh, const DataVector& t_nodalValues): mesh(&t_mesh), vec(t_nodalValues){
    assert(vec.size()==t_nodalValues.size());
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), mesh->nVertElement, mesh->numberElements);
    elementview.block(1,0,mesh->nVertElement-1,mesh->numberElements).rowwise() -= elementview.row(0);
};
template<int NDim> PolygonalFunction<NDim>::PolygonalFunction(const PolygonalMesh<NDim>& t_mesh, const std::function<Real(Point<NDim>)>& t_f): mesh(&t_mesh), vec(t_mesh.dataVectorDim){
    for (auto elem = mesh->elementIterator(); elem.unfinished ; ++elem) {
        for (auto vert = mesh->vertexInElemIterator(); vert.unfinished; ++vert){
            vec(vert.currentVertID + mesh->nVertElement * elem.currentElementID) = t_f(mesh->operator()(elem,vert));
        }
    }
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), mesh->nVertElement, mesh->numberElements);
    elementview.block(1,0,mesh->nVertElement-1,mesh->numberElements).rowwise() -= elementview.row(0);
}

template<int NDim> PolygonalFunction<NDim>::PolygonalFunction(const PolygonalMesh<NDim>* t_mesh): PolygonalFunction(*t_mesh) {}
template<int NDim> PolygonalFunction<NDim>::PolygonalFunction(const PolygonalMesh<NDim>* t_mesh, const Real t_value): PolygonalFunction(*t_mesh, t_value) {}
template<int NDim> PolygonalFunction<NDim>::PolygonalFunction(const PolygonalMesh<NDim>* t_mesh, const DataVector& t_nodalValues): PolygonalFunction(*t_mesh, t_nodalValues) {}
template<int NDim> PolygonalFunction<NDim>::PolygonalFunction(const PolygonalMesh<NDim>* t_mesh, const std::function<Real(Point<NDim>)>& t_f): PolygonalFunction(*t_mesh, t_f) {}

template<int NDim> PolygonalFunction<NDim>::PolygonalFunction(const std::function<Real(Real)>& t_f, const PolygonalFunction<NDim>& t_other): mesh(t_other.mesh), vec(t_other.vec){
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), mesh->nVertElement, mesh->numberElements);
    elementview.block(1,0,mesh->nVertElement-1,mesh->numberElements).rowwise() += elementview.row(0);
    for (auto elem = mesh->elementIterator(); elem.unfinished ; ++elem) {
        for (auto vert = mesh->vertexInElemIterator(); vert.unfinished; ++vert){
            vec(vert.currentVertID + mesh->nVertElement * elem.currentElementID) = t_f(vec(vert.currentVertID + mesh->nVertElement * elem.currentElementID));
        }
    }
    elementview.block(1,0,mesh->nVertElement-1,mesh->numberElements).rowwise() -= elementview.row(0);
};
template<int NDim> PolygonalFunction<NDim>::PolygonalFunction(const std::function<Real(std::vector<Real>)>& t_f, const std::vector<PolygonalFunction<NDim>>& t_others): mesh(t_others[0].mesh), vec(t_others[0].vec){
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

//********************************
//* Arithmetic
//********************************

template<int NDim> PolygonalFunction<NDim>& PolygonalFunction<NDim>::operator=(const PolygonalFunction<NDim>& other){
    assert(mesh==other.mesh);
    vec = other.vec;
    return *this;
}
template<int NDim> PolygonalFunction<NDim>& PolygonalFunction<NDim>::operator=(PolygonalFunction<NDim>&& other){  // Swap assignment
    assert(mesh==other.mesh);
//    vec.swap(other.vec);
    vec = std::move(other.vec);
    return *this;
}
template<int NDim> PolygonalFunction<NDim>& PolygonalFunction<NDim>::operator+=(const PolygonalFunction<NDim>& other){
    assert(mesh==other.mesh);
    vec+=other.vec;
    return *this;
};
template<int NDim> PolygonalFunction<NDim>& PolygonalFunction<NDim>::operator-=(const PolygonalFunction<NDim>& other){
    assert(mesh==other.mesh);
    vec-=other.vec;
    return *this;
};
template<int NDim> PolygonalFunction<NDim>& PolygonalFunction<NDim>::operator*=(const PolygonalFunction<NDim>& other){
    assert(mesh==other.mesh);
    toNodalFromLinForm();                                // Converts to nodal representation and stores in vec
    auto otherNodalRepr = other.toNodalFromLinForm();    // Converts other to nodal representation and stores in new matrix (since other is const)
    elementview().array() *= otherNodalRepr.array();// Executes the operation by treating the matrices as array (in Eigen sense)
    toLinFormFromNodal();                              // Converts back to lienar form representation and stores in vec
    return *this;
};
template<int NDim> PolygonalFunction<NDim>& PolygonalFunction<NDim>::operator/=(const PolygonalFunction<NDim>& other){
    assert(mesh==other.mesh);
    toNodalFromLinForm();                                // Converts to nodal representation and stores in vec
    auto otherNodalRepr = other.toNodalFromLinForm();    // Converts other to nodal representation and stores in new matrix (since other is const)
    elementview().array() /= otherNodalRepr.array();// Executes the operation by treating the matrices as array (in Eigen sense)
    toLinFormFromNodal();                              // Converts back to lienar form representation and stores in vec
    return *this;
};

template<int NDim> PolygonalFunction<NDim> PolygonalFunction<NDim>::operator-() const {
    Function<NDim> toreturn (*this);
    toreturn.vec *= -1.0;
    return toreturn;
}

template<int NDim> PolygonalFunction<NDim>& PolygonalFunction<NDim>::operator=(const Real scalar){
    vec = DataVector::Zero(mesh->dataVectorDim);
    elementview().col(0).array() = scalar;
    return *this;
}
template<int NDim> PolygonalFunction<NDim>& PolygonalFunction<NDim>::operator+=(const Real scalar){
    elementview().col(0).array() += scalar;
    return *this;
}
template<int NDim> PolygonalFunction<NDim>& PolygonalFunction<NDim>::operator-=(const Real scalar){
    elementview().col(0).array() -= scalar;
    return *this;
}
template<int NDim> PolygonalFunction<NDim>& PolygonalFunction<NDim>::operator*=(const Real scalar){
    vec *= scalar;
    return *this;
}
template<int NDim> PolygonalFunction<NDim>& PolygonalFunction<NDim>::operator/=(const Real scalar){
    vec /= scalar;
    return *this;
}
template<int NDim> PolygonalFunction<NDim> operator/(const Real lhs, PolygonalFunction<NDim> rhs) {
    rhs.toNodalFromLinForm();
    rhs.vec = rhs.vec.cwiseInverse();
    rhs.toLinFormFromNodal();
    return lhs*rhs;
}

template<int NDim> PolygonalFunction<NDim>& PolygonalFunction<NDim>::operator+=(const FunctionElement<NDim>& other){
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




} // namespace Tortoise
