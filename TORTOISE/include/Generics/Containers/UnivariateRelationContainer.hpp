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
//  UnivariateRelationContainer.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 3/11/20.
//
// This class is meant to allow for the treatment of functions of one variable y = f(x), discretised by providing y at discrete (non necessarily uniformly distributed)
// values of x. x is called the argument of the function and y is called the value. This naming convention is used throughout the class.
//
// This class is meant to provide a series of standardised methods to access the elements, to allow for the application of standardised ODE algorithms.
// The template parameters are the types of the containers to be used for the argument (x) and the value (y).
//  containerArg sould be of types like std::vector<something>, std::deque<something>
//  containerVal sould be of types like std::vector<somethingelse>, std::deque<somethingelse>
//  Other types of containers can be used but they must provide resizing capabilities


#ifndef UnivariateRelationContainer_hpp
#define UnivariateRelationContainer_hpp

namespace Tortoise {

template <class containerArg , class containerVal> class UnivariateRelationContainer {

    //=======================================================
    // Properties
    //===================
    containerArg    argVec;
    containerVal    valVec;

public:
    //=======================================================
    // Types
    //===================
    typedef     typename containerArg::value_type       ArgType;
    typedef     typename containerVal::value_type       ValType;
    
    //=======================================================
    // Constructors
    //===================
    // Note that it does not implement a default constructor, to allow for the use of non default constructible types
    UnivariateRelationContainer(const containerArg& argVec_, const containerVal& valVec_): argVec(argVec_), valVec(valVec_){assert(argVec_.size() == valVec_.size());};
    UnivariateRelationContainer(const ArgType& arg_, const ValType& val_): argVec({arg_}), valVec({val_}){};                    // Single entry constructor
    UnivariateRelationContainer(const UnivariateRelationContainer& other): argVec(other.argVec), valVec(other.valVec){};                    // Copy constructor
    UnivariateRelationContainer(UnivariateRelationContainer&& other): argVec(std::move(other.argVec)), valVec(std::move(other.valVec)){};   // Move constructor
    // https://stackoverflow.com/questions/3413470/what-is-stdmove-and-when-should-it-be-used

    //=======================================================
    // Assignment
    //===================
    // Move assignment https://en.cppreference.com/w/cpp/language/move_assignment
    UnivariateRelationContainer& operator=(const UnivariateRelationContainer& other){
        argVec = other.argVec; valVec = other.valVec; return *this;}                        // Copy assignment
    UnivariateRelationContainer& operator=(UnivariateRelationContainer&& other){
        argVec = std::move(other.argVec); valVec = std::move(other.valVec); return *this;}  // Move assignment

    //=======================================================
    // Status
    //===================
    bool                        empty() const {return argVec.empty();};
    int                         size() const {return argVec.size();};
    
    //=======================================================
    // Array subscripting
    //===================
    // Returns the argument of the function (i.e. the x variable)
    ArgType&                    arg(const unsigned long i) {return argVec[i];};                 // Array Subscripting x value of the numtime element in the Container
    const ArgType&              arg(const unsigned long i) const {return argVec[i];};           // Array Subscripting x value of the numtime element in the Container
    ArgType&                    arg(const int i) {return argVec[i];};                           // Array Subscripting x value of the numtime element in the Container
    const ArgType&              arg(const int i) const {return argVec[i];};                     // Array Subscripting x value of the numtime element in the Container
    ArgType&                    argfront() {return argVec.front();};                            // x of the first element
    const ArgType&              argfront() const {return argVec.front();};                      // x of the first element
    ArgType&                    argback() {return argVec.back();};                              // x of the last element
    const ArgType&              argback() const {return argVec.back();};                        // x of the last element
    
    // Returns the value of the function (i.e. the y variable)
    ValType&                    val(const unsigned long i) {return valVec[i];};                 // Array Subscripting y value of the numtime element in the Container
    const ValType&              val(const unsigned long i) const {return valVec[i];};           // Array Subscripting y value of the numtime element in the Container
    ValType&                    val(const int i) {return valVec[i];};                           // Array Subscripting y value of the numtime element in the Container
    const ValType&              val(const int i) const {return valVec[i];};                     // Array Subscripting y value of the numtime element in the Container
    ValType&                    valfront() {return valVec.front();};                            // y of the first element
    const ValType&              valfront() const {return valVec.front();};                      // y of the first element
    ValType&                    valback() {return valVec.back();};                              // y of the last element
    const ValType&              valback() const {return valVec.back();};                        // y of the last element

    // Same as above, simply shorter notation
    ValType&                    operator[](const unsigned long i) {return valVec[i];};          // Array Subscripting y value of the numtime element in the Container
    const ValType&              operator[](const unsigned long i) const {return valVec[i];};    // Array Subscripting y value of the numtime element in the Container
    ValType&                    operator[](const int i) {return valVec[i];};                    // Array Subscripting y value of the numtime element in the Container
    const ValType&              operator[](const int i) const {return valVec[i];};              // Array Subscripting y value of the numtime element in the Container
    ValType&                    front() {return valVec.front();};                               // y of the first element
    const ValType&              front() const {return valVec.front();};                         // y of the first element
    ValType&                    back() {return valVec.back();};                                 // y of the last element
    const ValType&              back() const {return valVec.back();};                           // y of the last element

    //=======================================================
    // Filling of new elements
    //===================
    UnivariateRelationContainer&     emplace_back(const ArgType& x, const ValType& y) {argVec.emplace_back(x); valVec.emplace_back(y); return *this;};
    
};

} // namespace Tortoise

#endif /* UnivariateRelationContainer_hpp */
