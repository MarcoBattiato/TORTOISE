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
//  MaterialStatus.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//

#include <Physics/Status/MaterialStatus.hpp>
#include <Geometry/Plotter/WeightedFunction.hpp>

#include <cassert>

using Tortoise::Containers::StringMap;

namespace Tortoise {

//=======================================================
// Constructors
//===================
template<int NDim> MaterialStatus<NDim>::MaterialStatus(Material<NDim> &t_material): material(&t_material) {
    for (int i=0; i < t_material.size(); i++) {
        this->emplace_back(t_material.id(i),Function<NDim>(*(t_material[i].mesh)));
    }
};
template<int NDim> MaterialStatus<NDim>::MaterialStatus(Material<NDim> *t_material): material(t_material) {
    for (int i=0; i < t_material->size(); i++) {
        this->emplace_back(t_material->id(i),Function<NDim>(*((*t_material)[i].mesh)));
    }
};
template<int NDim> MaterialStatus<NDim>::MaterialStatus(const MaterialStatus<NDim> &other) : material(other.material), StringMap<Function<NDim>>(other) {};

// The following move contructor has been built following
// https://stackoverflow.com/questions/15351341/move-constructors-and-inheritance
template<int NDim> MaterialStatus<NDim>::MaterialStatus(MaterialStatus<NDim>&& other) : material(other.material), StringMap<Function<NDim>>(std::move(other)) {};


//template<int NDim> void swap(MaterialStatus<NDim>& first, MaterialStatus<NDim>& second){
//    assert(first.material == second.material);
//    swap(first.populations, second.populations);
//};

//=======================================================
// Algebraic operations
//===================
template<int NDim> MaterialStatus<NDim>& MaterialStatus<NDim>::operator=(const MaterialStatus<NDim>& other){
    assert(material == other.material);
    return static_cast<MaterialStatus<NDim>&>(StringMap<Function<NDim>>::operator=(other));
}
// The following move assignment has been built following  https://stackoverflow.com/questions/15351341/move-constructors-and-inheritance
template<int NDim> MaterialStatus<NDim>& MaterialStatus<NDim>::operator=(MaterialStatus<NDim>&& other) {
    assert(material == other.material);
    return static_cast<MaterialStatus<NDim>&>(StringMap<Function<NDim>>::operator=(std::move(other)));
}

template<int NDim> MaterialStatus<NDim>& MaterialStatus<NDim>::operator+=(const MaterialStatus<NDim>& other){
    assert(material == other.material);
    for (int i=0; i < this->size(); i++) { (*this)[i]+=other[i];}
    return *this;
};
template<int NDim> MaterialStatus<NDim>& MaterialStatus<NDim>::operator-=(const MaterialStatus<NDim>& other){
    assert(material == other.material);
    for (int i=0; i < this->size(); i++) { (*this)[i]-=other[i];}
    return *this;
};
template<int NDim> MaterialStatus<NDim>& MaterialStatus<NDim>::operator*=(const MaterialStatus<NDim>& other){
    assert(material == other.material);
    for (int i=0; i < this->size(); i++) { (*this)[i]*=other[i];}
    return *this;
};
template<int NDim> MaterialStatus<NDim>& MaterialStatus<NDim>::operator/=(const MaterialStatus<NDim>& other){
    assert(material == other.material);
    for (int i=0; i < this->size(); i++) { (*this)[i]/=other[i];}
    return *this;
};
template<int NDim> MaterialStatus<NDim> MaterialStatus<NDim>::operator-() const {
    MaterialStatus<NDim> toreturn (*this);
    for (int i=0; i < this->size(); i++) { toreturn[i] *= -1.0;}
    return toreturn;
}

template<int NDim> MaterialStatus<NDim>& MaterialStatus<NDim>::operator=(const Real scalar){
    for (int i=0; i < this->size(); i++) { (*this)[i]=scalar;}
    return *this;
};
template<int NDim> MaterialStatus<NDim>& MaterialStatus<NDim>::operator+=(const Real scalar){
    for (int i=0; i < this->size(); i++) { (*this)[i]+=scalar;}
    return *this;
};
template<int NDim> MaterialStatus<NDim>& MaterialStatus<NDim>::operator-=(const Real scalar){
    for (int i=0; i < this->size(); i++) { (*this)[i]-=scalar;}
    return *this;
};
template<int NDim> MaterialStatus<NDim>& MaterialStatus<NDim>::operator*=(const Real scalar){
    for (int i=0; i < this->size(); i++) { (*this)[i]*=scalar;}
    return *this;
};
template<int NDim> MaterialStatus<NDim>& MaterialStatus<NDim>::operator/=(const Real scalar){
    for (int i=0; i < this->size(); i++) { (*this)[i]/=scalar;}
    return *this;
};
template<int NDim> MaterialStatus<NDim> operator/(const Real scalar, MaterialStatus<NDim> rhs){
    for (int i=0; i < rhs.size(); i++) { rhs[i]=scalar/rhs[i];}
    return rhs;
};

template<int NDim> MaterialStatus<NDim>& MaterialStatus<NDim>::applyInverseMass(){
    for (int i=0; i < this->size(); i++) { (*this)[i].applyInverseMass();}
    return *this;
};


//********************************
// Propagations
//********************************


template<int NDim> MaterialStatus<NDim> MaterialStatus<NDim>::propagate() const{
    MaterialStatus<NDim> toreturn(material);
    for (int i=0; i<material->scatteringChannelsList2Legs.size(); ++i){
        if (material->scatteringChannelsList2Legs[i].amplitude != 0.)
            toreturn += material->scatteringChannelsList2Legs[i].propagate(*this);
    }
    for (int i=0; i<material->scatteringChannelsList3Legs.size(); ++i){
        if (material->scatteringChannelsList3Legs[i].amplitude != 0.)
            toreturn += material->scatteringChannelsList3Legs[i].propagate(*this);
    }
    for (int i=0; i<material->scatteringChannelsList4Legs.size(); ++i){
        if (material->scatteringChannelsList4Legs[i].amplitude != 0.)
            toreturn += material->scatteringChannelsList4Legs[i].propagate(*this);
    }
    return toreturn;
}


template<int NDim> MaterialStatus<NDim> MaterialStatus<NDim>::scatteringRates() const{
    MaterialStatus<NDim> toreturn(material);
    assert(false);
    return toreturn;
};

template<> MaterialStatus<2> MaterialStatus<2>::scatteringRates() const{
    MaterialStatus<2> toreturn(material);
    for (auto band : material->matchListId("*")) {
        for (int i=0; i<material->scatteringChannelsList2Legs.size(); ++i){
            if (material->scatteringChannelsList2Legs[i].amplitude != 0.)
                toreturn[band] += material->scatteringChannelsList2Legs[i].scatteringRates(*this, band);
        }
        for (int i=0; i<material->scatteringChannelsList3Legs.size(); ++i){
            if (material->scatteringChannelsList3Legs[i].amplitude != 0.)
                toreturn[band] += material->scatteringChannelsList3Legs[i].scatteringRates(*this, band);
        }
        for (int i=0; i<material->scatteringChannelsList4Legs.size(); ++i){
            if (material->scatteringChannelsList4Legs[i].amplitude != 0.)
                toreturn[band] += material->scatteringChannelsList4Legs[i].scatteringRates(*this, band);
        }
    }
    return toreturn;
};

template<int NDim> Real MaterialStatus<NDim>::relativeError(const MaterialStatus<NDim>& population) const {
    Real  error = 0.;
    for (auto band: this->matchListId("*")){
        Real newError = (*material)[band].relativeError((*this)[band], population[band]);
        if (newError >= error) error = newError;
    }
    return error;
}

template<int NDim> MaterialStatus<NDim> MaterialStatus<NDim>::applyConstraints(Real time) const{
    MaterialStatus<NDim> toReturn(*this);
    for (auto band: this->matchListId("*")){
        if ((*material)[band].constrained) {
            toReturn[band] = (*material)[band].constrainedPopulation(time);
        }
    }
    return toReturn;
}

//********************************
//* Utilities
//********************************

template<int NDim> void MaterialStatus<NDim>::setToEquilibrium(Real temperature, Real chemicalPotentialFermions, Real chemicalPotentialBosons){
    for (auto band : material->matchListIndex("*")) {
        if ((*material)[band].statistics == fermion) {
            (*this)[band] = 1./(exp( ((*material)[band].dispersion - chemicalPotentialFermions)/temperature ) + 1.);
        } else {
            (*this)[band] = 1./(exp( ((*material)[band].dispersion - chemicalPotentialBosons)/temperature ) - 1.);
        }
    }
}

//********************************
//* I/O
//********************************


Plotter3D& operator << (Plotter3D& plotter, const MaterialStatus<2>& status) {
    // RED, BLUE, LIME, FUCHSIA, AQUA, YELLOW, MAROON, NAVY, GREEN
    static std::vector<std::string> colors  {"FF0000", "0000FF", "00FF00", "FF00FF", "00FFFF", "FFFF00", "800000", "000080", "008000"};
    plotter3d << *(*(status.material)).region  << LABEL("BZ") << COLOR("black") ;
    for (int band = 0; band < status.material->size(); ++band ) {
        plotter3d << *(*(status.material))[band].mesh << NOLABEL << COLOR("#AA" + colors[band%colors.size()]) ;
        plotter3d << WeightedFunction<2>((*(status.material))[band].dispersion,status[band])  << LABEL(status.material->id(band)) ;
    }
    return plotter;
}

template<> void MaterialStatus<2>::plot() const {plotter3d << *this << PLOT;};

template<> void MaterialStatus<2>::plot(const std::vector<std::string> &bandstoplot) const {
    // RED, BLUE, LIME, FUCHSIA, AQUA, YELLOW, MAROON, NAVY, GREEN
    static std::vector<std::string> colors  {"FF0000", "0000FF", "00FF00", "FF00FF", "00FFFF", "FFFF00", "800000", "000080", "008000"};
    plotter3d << *(*((*this).material)).region  << LABEL("BZ") << COLOR("black") ;
    for (int band = 0; band < bandstoplot.size(); ++band ) {
        plotter3d << *(*((*this).material))[band].mesh << NOLABEL << COLOR("#AA" + colors[band%colors.size()]) ;
        plotter3d << WeightedFunction<2>((*((*this).material))[band].dispersion,(*this)[band])  << LABEL((*this).material->id(band)) ;
    }
    plotter3d << PLOT;
};


// Explicit template instantiation (For more info see: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl)

template class MaterialStatus<1>;
template class MaterialStatus<2>;
template class MaterialStatus<3>;
//template void swap(MaterialStatus<1>& first, MaterialStatus<1>& second);
//template void swap(MaterialStatus<2>& first, MaterialStatus<2>& second);
//template void swap(MaterialStatus<3>& first, MaterialStatus<3>& second);
//template MaterialStatus<1> propagateDeterministic(Real, const MaterialStatus<1>& status);
//template MaterialStatus<2> propagateDeterministic(Real, const MaterialStatus<2>& status);
//template MaterialStatus<3> propagateDeterministic(Real, const MaterialStatus<3>& status);

} // namespace Tortoise
