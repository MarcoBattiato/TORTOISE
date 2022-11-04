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
//  MaterialTimeStatus.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//

#ifndef MaterialTimeStatus_hpp
#define MaterialTimeStatus_hpp

#include <Physics/Status/MaterialStatus.hpp>
#include <Generics/Containers/UnivariateRelationContainer.hpp>


#include <deque>

namespace Tortoise {

template<int NDim> class MaterialTimeStatus :

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Inheritance
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public
UnivariateRelationContainer<
    std::deque<Real>,              // Times
    std::deque<MaterialStatus<NDim>>        // MaterialStatus<NDim> containing the populations at different times
>

{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main properties
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    const Material<NDim> *          material;       // Pointer to the material the populations refer to
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main methods
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    //=======================================================
    // Constructors
    //===================
    explicit MaterialTimeStatus(const Material<NDim> &Materialpass);
    explicit MaterialTimeStatus(const MaterialStatus<NDim>& matStatus0); // Assumes that initial time is 0.
    MaterialTimeStatus(Real time0, const MaterialStatus<NDim>& matStatus0);

    
    //********************************
    // Propagations
    //********************************
    
    void propagateDeterministic(double finalTime, double denseStep);    // Propagates with adaptive DP5(4)
    
    MaterialTimeStatus(const Material<NDim>&& Materialpass) = delete;
};
    
    
} // namespace Tortoise

#endif /* MaterialTimeStatus_hpp */
