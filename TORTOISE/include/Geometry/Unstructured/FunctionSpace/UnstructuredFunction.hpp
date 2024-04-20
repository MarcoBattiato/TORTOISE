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
//  UnstructuredFunction.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 7/8/21.
//
//  INCOMPLETE

#ifndef UnstructuredFunction_hpp
#define UnstructuredFunction_hpp


#include <Generics/Features/VectorSpace.hpp>
#include <Generics/Features/MathFieldSpace.hpp>
#include <Geometry/Unstructured/DirectSpace/UnstructuredMesh.hpp>
#include <Geometry/Plotter/Plotter.hpp>

#include <functional>
#include <math.h>

namespace Tortoise {

class UnstructuredFunction1D :
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Inheritance
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
// Inherits all the operations of a field from MathFieldSpace and AsymmetricMathFieldSpace
// Math field operations are defined (+,-,*,/).
// The same operations are defined with scalars.
// Only + and - are defined on FunctionElement. This is because it would be unclear what the user might expect to happen on the remaining
// elements. The algebraic operations are done assuming that the FunctionElement is non zero only on the specified element. However if the same
// interpretation were to be used for * and /, problems arise (especially in the case of /). If the user wants to do those operations he/she should first
// construct a Function from the FunctionElement and explicitly set the value everywhere else to avoid confusion.
// This inheritance constructs all the possible combinations of operators
public Features::MathFieldSpace<UnstructuredFunction1D>,
public Features::AsymmetricMathFieldSpaceByValue<UnstructuredFunction1D,Real> {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main properties
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    const UnstructuredMesh1D*   mesh;
    DataVector                  vec;
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main methods
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    //=======================================================
    // Constructors
    //===================
    explicit UnstructuredFunction1D(const UnstructuredMesh1D& mesh);                                                  // Constructs a constant 0 function
    UnstructuredFunction1D(const UnstructuredMesh1D& mesh, const Real value);                                // Constructs a constant function
    UnstructuredFunction1D(const UnstructuredMesh1D& mesh, const DataVector& nodalValues);                            // Constructs using the values at nodes
    UnstructuredFunction1D(const UnstructuredMesh1D& mesh, const std::function<Real(Real)>& f);     // Constructs using the function passed
    explicit UnstructuredFunction1D(const UnstructuredMesh1D* mesh);                                                  // Constructs a constant 0 function
    UnstructuredFunction1D(const UnstructuredMesh1D* mesh, const Real value);                                // Constructs a constant function
    UnstructuredFunction1D(const UnstructuredMesh1D* mesh, const DataVector& nodalValues);                            // Constructs using the values at nodes
    UnstructuredFunction1D(const UnstructuredMesh1D* mesh, const std::function<Real(Real)>& f);       // Constructs using the function passed
    UnstructuredFunction1D(const std::function<Real(Real)>& f, const UnstructuredFunction1D& other);// Constructs a function of a function f(other)
    UnstructuredFunction1D(const std::function<Real(std::vector<Real>)>& f, const std::vector<UnstructuredFunction1D>& t_others);// Constructs a function of a function f(other[0], other[1], ...)
    UnstructuredFunction1D(const UnstructuredFunction1D& t_other);                                                        // Copy constructor
    UnstructuredFunction1D(UnstructuredFunction1D&& t_other);                                                             // Move constructor
    
        
        //=======================================================
        // Implementation Details
        //===================
public:
    auto elementview();         // Returns a map. Can be used to modify directly vec, but it is accessed in element view
    auto elementview() const ;  // Returns a map. Cannot be used to modify directly vec, since it is constant

};


//=======================================================
// Implementation Details
//===================

auto UnstructuredFunction1D::elementview() {
    return Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vec.data(),  mesh->numberElements, 2);
}
auto UnstructuredFunction1D::elementview() const {
    return Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vec.data(), mesh->numberElements, 2);
}

} // namespace Tortoise

#endif /* UnstructuredFunction_hpp */
