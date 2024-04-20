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
//  Band.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//

#ifndef Band_hpp
#define Band_hpp

#include <Geometry/Structured/FunctionSpace/Function.hpp>
#include <Generics/Containers/StringMap.hpp>

#include <vector>
#include <array>
#include <functional>

namespace Tortoise {

template<int NDim> class Band {                                 // NDim = Number of Dimensions
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main properties
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    const Mesh<NDim> *                      mesh;               // Pointer to the mesh object (if the Mesh object gets deallocated, the behaviour is undefined)

// Scalar attributes
    Real                                    charge;             // Stores the value of the charge and the multiplier for the statistic
    Real                                    statistics;         // Stores the value of the multiplier for the statistic
    Containers::StringMap<Real>             attribute;          // Used to store user-defined scalar attributes.
    
// Function attributes
    Function<NDim>                          dispersion;         // Stores the dispersion
    Function<NDim>                          densityOfStates;    // Stores the density of states (at this stage a simple constant function)
    std::array<Function<NDim>,NDim>         crystalMomentum;    // Stores the crystal momentum (at this stage a simple f(k)=kx )
    Containers::StringMap<Function<NDim>>   functionAttribute;  // Used to store attributes that are functions (for instance spin polarisation)
    
// Time propagation properties
    bool                                    constrained = false;// If constrained, it means that the population in the band should not be allowed to propagate, but should have a value fixed by the user
    // This is useful, for instance, to
    // 1) contrain the population within the present band to be fixed (for instance in bands that should be used as reservoirs)
    // 2) contrain the population to have a specific value at specific times (useful to define the temporal and spectral profile of a laser)
    std::function<Function<NDim>(Real)>     constrainedPopulation; // Function of time that gives the population at the specified time
    
// Parameters for error calculation in adaptive time step time propagation
    Real                                    absoluteErrorTolerance = 1.e-3;         // Used in adaptive time step: acceptable absolute error in population for this band
    Real                                    relativeErrorTolerance = 1.e-3;         // Used in adaptive time step: acceptable relative error in population as fraction of population for this band
    // It is discouraged to set too high tolerances. Given the high order of convergence of the used time stepping methods, the savings
    // when setting high tolerances are usually marginal, yet high tolerances can lead to accepted errors that can be even higher than
    // the population in the bands.
    Real                                    penalisationUnphysicalPopulation = 0.;  // Amplification factor of contribution to error due to appearence of unphysical population for this band
                                                                                    // A population is unphysical if <0.0 and >1.0 for fermions, or if <0.0 for bosons
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main methods
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    //=======================================================
    // Constructors
    //===================
    Band(const Real t_charge, const Real t_statistics, const Function<NDim>& t_dispersion);
    Band(const Real t_charge, const Real t_statistics, const Mesh<NDim> &t_mesh, const std::function<Real(Point<NDim>)>& t_dispersion);
    Band(const Real t_charge, const Real t_statistics, const Mesh<NDim> &t_mesh, const Real t_dispersion);
    
    //=======================================================
    // Stepper errors
    //===================
    Real relativeError(const Function<NDim>& error, const Function<NDim>& population) const;
    
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Secondary properties (implementation details)
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public: // Implementation details
    // Store max and min of dispersion within each element for quick lookup
    DataVector maxElemEnergies, minElemEnergies;
     
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Friend methods
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
template<int NDim> bool operator==(const Band<NDim>& lhs, const Band<NDim>& rhs);       // Compares only the ID!!!!


} // namespace Tortoise
    
#endif /* Band_hpp */
