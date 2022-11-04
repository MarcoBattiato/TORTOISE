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
//  MaterialTimeStatus.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//

#include <Physics/Status/MaterialTimeStatus.hpp>
#include <Generics/Algorithms/DormandPrince54.hpp>


namespace Tortoise {
// Explicit template instantiation (For more info see: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl)

//=======================================================
// Constructors
//===================
template<int NDim> MaterialTimeStatus<NDim>::MaterialTimeStatus(const Material<NDim> &Materialpass):
UnivariateRelationContainer<std::deque<Real>, std::deque<MaterialStatus<NDim>>>(std::deque<Real>{}, std::deque<MaterialStatus<NDim>>{}), material(&Materialpass){}
template<int NDim> MaterialTimeStatus<NDim>::MaterialTimeStatus(Real time0, const MaterialStatus<NDim>& matStatus0):
UnivariateRelationContainer<std::deque<Real>, std::deque<MaterialStatus<NDim>>>(time0, matStatus0), material(matStatus0.material){}
template<int NDim> MaterialTimeStatus<NDim>::MaterialTimeStatus(const MaterialStatus<NDim>& matStatus0):
MaterialTimeStatus(0.0, matStatus0){}


//********************************
// Propagations
//********************************

template<int NDim> void MaterialTimeStatus<NDim>::propagateDeterministic(double finalTime, double denseStep){
    dormandPrince54(
                    // input/output
                    *this,
                    // end time of propagation
                    finalTime,
                    // step for dense grid output
                    denseStep,
                    // operator that gives the derivative
                    // It must be >   timeOp(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> ContainerType::ValType
                    [](Real t, const MaterialStatus<NDim>& status){
                        return status.propagateDeterministic();
                    },
                    // Function that modifies the value of the solution to enforce certain constraints
                    // it must be >   enforceValue(ContainerType::ArgType currentTime, const ContainerType::ValType& currentY) -> ContainerType::ValType
                    [](Real t, const MaterialStatus<NDim>& status)-> MaterialStatus<NDim> {
                        return status.applyConstraints(t);
                    },
                    //  normalised error
                    // It must be >    normErrorForm(ContainerType::ValType localErrY, ContainerType::ValType currentY) -> ContainerType::ArgType
                    [](const MaterialStatus<NDim>& error, const MaterialStatus<NDim>& currentY){
                        return error.relativeError(currentY);
                    },// Defines some output operations
                    // It must be in the form >     outputFunct(ContainerType::ArgType currentTime, ContainerType::ValType currentY)
                    [](Real t, const MaterialStatus<NDim>& currentY){}
                    );
};




template class MaterialTimeStatus<1>;
template class MaterialTimeStatus<2>;
template class MaterialTimeStatus<3>;

} // namespace Tortoise
