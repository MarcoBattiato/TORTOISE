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
//  MaterialStatus.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//

#ifndef MaterialStatus_hpp
#define MaterialStatus_hpp

#include <Physics/Material/Material.hpp>

namespace Tortoise {

template<int NDim> class Material;

template<int NDim> class MaterialStatus :

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Inheritance
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
// Material status is a container of populations (Function<NDim>)
// Populations of each single band can be accessed by materialstatus["bandname"]
// Populations can be iterated over using wildcards
// for (auto i : materialstatus.matchListIndex("electr*")) { materialstatus[i];}
// or
// for (auto i : materialstatus.matchList("electr*")) { *i;}
public Containers::StringMap<Function<NDim>>,   // It contains the populations

// Math field operations are defined (+,-,*,/).
// The same operations are defined with scalars.
// This inheritance constructs all the possible combinations of operators
public Tortoise::Features::MathFieldSpace<MaterialStatus<NDim>> ,   // Inherits all the operations of a field from MathFieldSpace and AsymmetricMathFieldSpace
public Tortoise::Features::AsymmetricMathFieldSpaceByValue<MaterialStatus<NDim>,Real>
{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main properties
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    Material<NDim> * const          material;       // Pointer to the material the populations refer to

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main methods
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    //=======================================================
    // Constructors
    //===================
    explicit MaterialStatus(Material<NDim> &t_material);          // Initialises all the populations to zero
    explicit MaterialStatus(Material<NDim>* const t_material);          // Initialises all the populations to zero
    MaterialStatus(const MaterialStatus<NDim> &other);                  // Copy constructor
    MaterialStatus(MaterialStatus<NDim>&& other);                       // Move constructor
    
    //********************************
    //* Arithmetic
    //********************************
    // ATTENTION!!!!!!! Multiplications and divisions by anything except a scalar do not conserve precision!
    // If you don't know what this means, do not use those operations. You will get unexpected results in many cases
    MaterialStatus& operator=(const MaterialStatus<NDim>& other);
    MaterialStatus& operator=(MaterialStatus<NDim>&& other);
    MaterialStatus& operator+=(const MaterialStatus<NDim>& other);
    MaterialStatus& operator-=(const MaterialStatus<NDim>& other);
    MaterialStatus& operator*=(const MaterialStatus<NDim>& other);                          // !!!! Loses precision. DO NOT USE if you do not know what this means!!!!
    MaterialStatus& operator/=(const MaterialStatus<NDim>& other);                          // !!!! Loses precision. DO NOT USE if you do not know what this means!!!!

    MaterialStatus operator-() const;
        
    MaterialStatus& operator=(const Real scalar);
    MaterialStatus& operator+=(const Real scalar);
    MaterialStatus& operator-=(const Real scalar);
    MaterialStatus& operator*=(const Real scalar);
    MaterialStatus& operator/=(const Real scalar);
    friend MaterialStatus operator/(const Real lhs, MaterialStatus<NDim> rhs);
        
    MaterialStatus& applyInverseMass();                            // Multiply all populations by appropriate inverse mass
    
    //********************************
    // Propagations
    //********************************
    MaterialStatus propagate() const;
    MaterialStatus scatteringRates() const;
    
    Real           relativeError(const MaterialStatus<NDim>& population) const;         // Assumes that the object is the error and the population at the next time step is passed
    
    MaterialStatus applyConstraints(Real time) const;
    
    //********************************
    // Utilities
    //********************************
    void setToEquilibrium(Real temperature, Real chemicalPotentialFermions, Real chemicalPotentialBosons); 
                                                                     
    
    //********************************
    //* I/O
    //********************************
    void save(const std::string &t_fileName) const;                                 // Save status to file
    void load(const std::string &t_fileName) ;                                      // Load status from file

    void plot() const;
    void plot(const std::vector<std::string> &bandstoplot) const;
//    void plot(const std::vector<int> &bandstoplot) const;
    
    MaterialStatus(const Material<NDim>&& t_material) = delete;

};
Plotter3D& operator << (Plotter3D& plotter, const MaterialStatus<2>& material);


//=======================================================
// Constructors: memory managers
//===================
template<int NDim> void swap(MaterialStatus<NDim>& first, MaterialStatus<NDim>& second); // Swap function

template<int NDim> auto sin(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=sin(other[i]);} return other; }
template<int NDim> auto cos(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=cos(other[i]);} return other; }
template<int NDim> auto tan(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=tan(other[i]);} return other; }
template<int NDim> auto asin(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=asin(other[i]);} return other; }
template<int NDim> auto acos(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=acos(other[i]);} return other; }
template<int NDim> auto atan(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=atan(other[i]);} return other; }
template<int NDim> auto cosh(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=cosh(other[i]);} return other; }
template<int NDim> auto sinh(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=sinh(other[i]);} return other; }
template<int NDim> auto tanh(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=tanh(other[i]);} return other; }
template<int NDim> auto pow(MaterialStatus<NDim> other, Real exponent) { for (int i=0; i < other.size(); i++) { other[i]=pow(other[i], exponent);} return other; }
template<int NDim> auto sqrt(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=sqrt(other[i]);} return other; }
template<int NDim> auto log(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=log(other[i]);} return other; }
template<int NDim> auto floor(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=floor(other[i]);} return other; }
template<int NDim> auto ceil(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=ceil(other[i]);} return other; }
template<int NDim> auto abs(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=abs(other[i]);} return other; }
template<int NDim> auto fabs(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i]=fabs(other[i]);} return other; }
template<int NDim> auto exp(MaterialStatus<NDim> other) { for (int i=0; i < other.size(); i++) { other[i] = exp(other[i]);} return other; }

} // namespace Tortoise


#endif /* MaterialStatus_hpp */
