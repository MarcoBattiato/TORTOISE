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
//  main.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 16/12/21.
//
// TUTORIAL 4: Functions
//
// Some parts of the tutorial are marked with >>> DEVELOPERS <<<
// That means that the final user most probably will not need to know about those functionalities


// #define NDEBUG            // Uncomment to deactivate debug mode

#include <TORTOISE>
#include <iostream>
#include <functional>

using namespace Tortoise;
using std::cout;

int main(int argc, const char * argv[]) {
    
    //**********************
    // Function (Function.hpp)
    //**********************
    // Functions are more complex objects and provide a large amount of functionalities
    // Almost all mathematical operations relevant to the Boltzmann equation are already implemented.
    // If you find some mathematical operation not implemented yet that you need to perform, let me know.
    
    //======================
    // Constructors
    //======================
    // Function objects can be initialised in several ways
    Function<2> funct0(myMesh);         // Creates a function object on the Mesh mesh, and initialises it to constant 0.0
    Function<2> funct1(myMesh, 3.4);    // Creates a function object on the Mesh mesh, and initialises it to constant 3.4
    
    // We can also initialise Function<n> objects by providing an analytic expression to discretise.
    // Notice that the analytic expression is used only for the construction of the discretised object and then discarded.
    // We can construct using lambda expressions
    Function<2> funct2(myMesh, [](Point<2> P){return P(0)+P(1)*P(1);});  // Creates a function object on the Mesh mesh, and initialises it to the discretised version of x+y^2
    // But also TORTOISE provides functionalities to pass nicer looking expressions using x, y, and z,
    // or kx, ky, kz and k (with the latter intended as the vector of coordinates)
    Function<2> funct3(myMesh, x + pow(y,2));
    // Please notice that it is not advisable to use the objects x, y, z, kx, ky, kz and k beyond the usage suggested below, unless the user is not fully aware of
    // what these objects are. They are actually functors built with expression templates, to make their evaluation as efficient as possible. In particular expression
    // templates require some understanding, to be used properly beyond what is suggested here.
    
    // Notice that, while regions and meshes cannot be modified after being created, functions can be modifed
    
    // There are more advanced constructor
    
    // This construct a new discretised function as the result of some expression of another discretised function
    Function<2> funct4([](realValueType z){return 1./z;}, funct2);
    Function<2> funct4_1(1./x, funct2);
    // Notice nhow it is not necessary to pass a mesh, since this constructor creates a function object on the same
    // mesh of funct2, and initialises it to 1/funct2(x,y)
    
    // As a function of more variables
    Function<2> funct5([](std::vector<realValueType> z){return 1./z[0] + 3.*z[0]*z[1];}, {funct2, funct3});
    // Creates a function object on the same mesh of funct2, and initialises it to 1/funct2(x,y) + funct2(x,y) * funct3(x,y)
    // Notice that if funct2 and funct3 are defined on different meshes, the constructor will return an error
    Function<2> funct5_1(1./x + 3.* x * y, {funct2, funct3});
    
    // Of course the usual copy constructor is provided
    Function<2> funct6(funct3);    // Makes a new function, copying funct4
    
    //======================
    // Evaluation
    //======================
    // A Function can be evaluated at a point. However, given that the function can be discontinuous at the elements' edges,
    // over there the result depends on which element one is looking at.
    // Therefore there is no method to evaluate a function at a point in global coordinates
    // The evaluation of a function is done by specifying the element and then the local coordinate of the point
    for (auto elem = myMesh.elementIterator(); elem.unfinished; ++elem) {
        cout << funct4(elem, Point<2>({0.333,0.333})) << " ";     // Value of the function at the baricentre of every element
    }
    cout << "\n";
    
    for (auto elem = myMesh.elementIterator(); elem.unfinished; ++elem) {
        for (auto vert = myMesh.vertexInElemIterator(); vert.unfinished; ++vert) {
            cout << funct4(elem, vert) << " ";     // Value of the function at each vertex of every element
        }
        cout << "\n";
    }
    cout << "\n";
    
    //======================
    // Mathematical Operations
    //======================
    // Function can be operated upon using overloaded operators.
    // All operators have been overloaded to provide the most natural platform for arithmetic operations and assignments
    // All these operations can operate on functions defined on the same mesh. If the meshes are not the same a compilation error is produced if debug mode is active.
    funct0 = funct1 + 3.5 * funct2 - 1.3;      // Calculates the linear combination and assign to funct0.
    funct0 += 2.3;                             // For convenience also the sum of a Function and a scalar has been overloaded
    
    // Notice that all the methods can be applied to temporary Function
    (funct0 + 2.4 * funct3).plot();
    Function<2> newfunc(funct3 + funct0 / 0.3);

    // Multiplication and division between functions are also implemented.
    // However notice that these operations cause a LOSS OF PRECISION!!! If you don't know what this means, then use with caution!
    funct0 = funct1 * 3.5 * funct2 / funct3;
    
    // It is also possible to use common mathematical expressions. However notice that these operations cause a LOSS OF PRECISION!!!
    // For a full list, look at the Function.hpp file
    funct0 = funct1 * sin(funct2) / ( 1. + exp(funct0));
    
    // It is possible to find the global maximum and minimum
    cout << funct0.max() << " " << funct0.min() << "\n\n";
    
    // Integration has a special importance.
    // TORTOISE provides methods to perform integration without loss of information
    cout << funct0.integrate() << "\n\n";
    // This integrates with respect of all the variable of the function over the mesh
    
    // In case it is necessary to perform the integral of the product of two functions (we will see that it is frequently required) TORTOISE
    // provides a method to perform this integral with no loss of information
    cout << funct0.integrate(funct1) << "\n\n";
    // Notice that if you perform the integration in the way below
    cout << (funct0 * funct1).integrate() << "\n\n";    // Do not do this!!!
    // the result is different
    
    // Same holds for the integral of the product of three functions
    cout << funct0.integrate(funct1,funct2) << "\n\n";

    
    
    return 0;
}
