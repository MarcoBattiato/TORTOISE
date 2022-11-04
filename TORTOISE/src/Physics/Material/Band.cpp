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
//  Band.cpp
//  Tortoise
//
//  Created by Marco Battiato on 23/11/19.
//

#include <Physics/Material/Band.hpp>

#include <cassert>

namespace Tortoise {

template<int NDim> std::array<Function<NDim>,NDim> costructCrystalMomentum(const Mesh<NDim> &t_mesh);

template<> std::array<Function<1>,1> costructCrystalMomentum<1>(const Mesh<1> &t_mesh){
    std::array<Function<1>,1> toreturn = { Function<1>(t_mesh,[](Point<1> P) {return P(0);})};
    return toreturn;
}
template<> std::array<Function<2>,2> costructCrystalMomentum<2>(const Mesh<2> &t_mesh){
    std::array<Function<2>,2> toreturn = { Function<2>(t_mesh,[](Point<2> P) {return P(0);}), Function<2>(t_mesh,[](Point<2> P) {return P(1);})};
    return toreturn;
}
template<> std::array<Function<3>,3> costructCrystalMomentum<3>(const Mesh<3> &t_mesh){
    std::array<Function<3>,3> toreturn = { Function<3>(t_mesh,[](Point<3> P) {return P(0);}), Function<3>(t_mesh,[](Point<3> P) {return P(1);}), Function<3>(t_mesh,[](Point<3> P) {return P(2);})};
    return toreturn;
}

template<int NDim> Band<NDim>::Band(const Real t_charge, const Real t_statistics, const Function<NDim>& t_dispersion)
: charge(t_charge), statistics(t_statistics), mesh(t_dispersion.mesh), dispersion(t_dispersion), densityOfStates(Function<NDim>(*t_dispersion.mesh,1.0)), crystalMomentum(costructCrystalMomentum<NDim>(*t_dispersion.mesh)) {
    minElemEnergies.resize(1, mesh->numberElements); maxElemEnergies.resize(1, mesh->numberElements);
    for (auto iter = dispersion.elementIterator(); iter.unfinished; ++iter){
        minElemEnergies(iter.currentElementID) = dispersion.elementMin(iter);
        maxElemEnergies(iter.currentElementID) = dispersion.elementMax(iter);
    }
};
template<int NDim> Band<NDim>::Band(const Real t_charge, const Real t_statistics, const Mesh<NDim> &t_mesh, const std::function<Real(Point<NDim>)>& t_dispersion)
: charge(t_charge), statistics(t_statistics), mesh(&t_mesh), dispersion(t_mesh,t_dispersion), densityOfStates(Function<NDim>(t_mesh,1.0)), crystalMomentum(costructCrystalMomentum<NDim>(t_mesh)) {
    minElemEnergies.resize(1, mesh->numberElements); maxElemEnergies.resize(1, mesh->numberElements);
    for (auto iter = dispersion.elementIterator(); iter.unfinished; ++iter){
        minElemEnergies(iter.currentElementID) = dispersion.elementMin(iter);
        maxElemEnergies(iter.currentElementID) = dispersion.elementMax(iter);
    }
};
template<int NDim> Band<NDim>::Band(const Real t_charge, const Real t_statistics, const Mesh<NDim> &t_mesh, const Real t_dispersion)
: charge(t_charge), statistics(t_statistics), mesh(&t_mesh), dispersion(t_mesh,t_dispersion), densityOfStates(Function<NDim>(t_mesh,1.0)), crystalMomentum(costructCrystalMomentum<NDim>(t_mesh)) {
    minElemEnergies.resize(1, mesh->numberElements); maxElemEnergies.resize(1, mesh->numberElements);
    for (auto iter = dispersion.elementIterator(); iter.unfinished; ++iter){
        minElemEnergies(iter.currentElementID) = dispersion.elementMin(iter);
        maxElemEnergies(iter.currentElementID) = dispersion.elementMax(iter);
    }
};

template<int NDim> bool operator==(const Band<NDim>& lhs, const Band<NDim>& rhs){
    return lhs.ID == rhs.ID;
}       // Compares only the ID

template<int NDim> Real Band<NDim>::relativeError(const Function<NDim>& error, const Function<NDim>& population) const {
    assert(error.mesh == population.mesh);
    
    if (statistics == -1) {
        Function<NDim> scale (relativeErrorTolerance * min(fabs(population), fabs(1-population)) + absoluteErrorTolerance);
        Function<NDim> relErrorFunct (error/scale);
        Function<NDim> penalizUnphysical ( (-min(population, 0.) + max(population-1,0.))/scale );
        return std::sqrt((relErrorFunct.integrate(relErrorFunct) + penalisationUnphysicalPopulation * penalizUnphysical.integrate(penalizUnphysical))/(population.mesh->numberElements*population.mesh->elemVolume));
    } else {
        Function<NDim> scale (relativeErrorTolerance * population + absoluteErrorTolerance);
        Function<NDim> relErrorFunct (error/scale);
        Function<NDim> penalizUnphysical ( -min(population, 0.)/scale );
        return std::sqrt((relErrorFunct.integrate(relErrorFunct) + penalisationUnphysicalPopulation * penalizUnphysical.integrate(penalizUnphysical))/(population.mesh->numberElements*population.mesh->elemVolume));
    }
}


// Explicit template instantiation (For more info see: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl)

template class Band<1>;
template class Band<2>;
template class Band<3>;

} // namespace Tortoise
