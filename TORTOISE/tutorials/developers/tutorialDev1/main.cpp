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
// TUTORIAL 1: Points and lines in TORTOISE
//
// Some parts of the tutorial are marked with >>> DEVELOPERS <<<
// That means that the final user most probably will not need to know about those functionalities



// ******************************
// How to run
// ******************************

// TORTOISE requires -O3 optimisation (as a large part of the code is built to trigger certain optimisations)
// It also requires to be compiled with c++20, since TORTOISE uses several features introduced in c++20.


// ******************************
// Debug mode
// ******************************

// TORTOISE can run in two modes: debug and release.
// In debug mode several errors not caught at compile time can be caught at runtime through the use of assert (https://www.cplusplus.com/reference/cassert/assert/)
// Debug mode is activated by default. From the tests I did it seems that all the included asserts are very cheap.
// Yet it is possible to deactivate debug mode by uncommenting the line below

// #define NDEBUG            // Uncomment to deactivate debug mode

// Be careful as several tests will not be performed when debug mode is deactivated and bugs may remain hidden.
// It is strongly suggested to let debug mode on while writing code and performing tests, and deactivate it only for the final high resolution calculation.


// ******************************
// Includes
// ******************************

// TORTOISE defines a large number of classes over several files.
// For these definitions to be available to the main, those files need to be included (same as you need to include iostream to be able to use cout).
// Given the large number of files and the fact that there are some dependencies, three files have been provided to simplify the inclusion of correct
// headers
// Includes/DiscretizationCore.hpp
// Includes/PhysicsCore.hpp
// Includes/MiscUtilities.hpp

// DiscretizationCore.hpp containes includes to the files that collectively provide geometric, plotting and function discretisation functionalities,
// but are not necesarily related to the Boltzmann equation
// PhysicsCore.hpp containes includes to the files that collectively provide functionalities directly related to the Boltzmann equation. In particular,
// how to construct materials, define scattering channels and do time propagation.
// MiscUtilities.hpp is a collaction of utilities that I programmed and turned out useful for a number of uses of TORTOISE. We will see later what is provided.


#include "DiscretizationCore.hpp"
#include <iostream>

// All the classes in TORTOISE are within the namespace Totroise. Therefore all classes and objects within TORTOISE need to be called with Tortoise::name_of_class
// That will be very cumbersome.
// It is usually a bad idea to use the keyword "using" with an entire namespace, unless that namespace
// is the key one. Therfore in this case we will do
using namespace Tortoise;
// There is little risk involved in this since the names within TORTOISE tend to be rather specific

// Notice how intead we do not do "using namespace std". That really is a very bad practice. But intead we do that only for the commands we know we are going to use
using std::cout;


int main(int argc, const char * argv[]) {
    
    // First we learn few fundamental functionalities offered by "Geometry.hpp"
    
    // ************************************************************
    // The type alias: Real (Geometry.hpp)
    // ************************************************************

    // The whole library is built in a way that it is easy to change the precision of the used real numbers.
    // Instead of explicitly using float or double, etc, the user should use the type alias: Real
    // In the version this tutorial was written for, Real = double
    // However in a later version the precision might be decreased or increased.

    // Do not use
    // double myValue = 1.35;
    // Use instead
    Real myValue = 1.35;

    // This will ensure that all the numbers within TORTOISE have the same precision
    
    // If for whatever reason you want to change the type aliased by Real
    // modify file Geometry/GeometryCore/Geometry.hpp  line 111
    // Notice that this change will automatically propagate everywhere in TORTOISE.
    // This is a very convenient way of changing precision everywhere by changing a single line.
    
    
    // ************************************************************
    // Points (Geometry.hpp)
    // ************************************************************

    // The library offers definition of types used to define points in real space
    // All the types are templated, meaning that they need to have some template values specified.
    // Usually the first template parameter is the number of dimensions of the space
    // A point can be initialised easily     Point<2> exmplPoint = {1.2,-2.5};
    
    cout << "Example 1.1 \n";
    Point<2> exmplPoint = {1.2,-2.5};
    cout << exmplPoint << "\n\n";

    // Note that standard operations can be conducted on the points
    cout << "Example 1.2 \n";
    Point<2> exmplPoint2 = {1.5, -0.5};
    cout << exmplPoint + (0.5 + myValue) * exmplPoint2 << "\n\n";

    // Notice that you cannot sum to a list, you need to construct a Point first
    cout << "Example 1.3 \n";
    cout << exmplPoint + 0.5 * Point<2>({1.0, 1.0}) << "\n\n";

    // It is possible to extract the cartesian components of the Point as
    cout << "Example 1.4 \n";
    cout << exmplPoint(0)<< " " << exmplPoint(1) << "\n\n";
    // The x component is the compoentnt 0, y is the component 1, and z the component 2 (of course if they exist, meaning if the
    // specified dimension is sufficient)
    
    
    // ************************************************************
    // Arrays of Points (Geometry.hpp) >>> DEVELOPERS <<<
    // ************************************************************

    // In case one needs to store several points, one can use two types of containers:
    // ArrayPoint: for a small number of Points
    // VectorPoint: for a large number of Points

    // ArrayPoint requires two template arguments: dimension of the space and number of points
    // There are two ways to initialise an array of points
    cout << "Example 2.1 \n";
    ArrayPoint<2,4> exmplPointArray ({{1.2,  1.5,  1.,  1.5},     // x components of all the points
                                     {-2.5, -0.5, 0.5, -0.5}});   // y components of all the points
    cout << exmplPointArray << "\n\n";
    // The one above stores 4 points in 2 dimensions. Notice the structure of the initialisation that requires to pass the x coordinates
    // for all the points and then the y coordinates. In 3D it would also require the z coordinates.
    
    // The other possibility is to use the sintax below
    // Notice that the order is: all the x coordinates for all the points, then the y coordinates for all the points and so on.
    cout << "Example 2.2 \n";
    ArrayPoint<2,4> exmplPointArray2nd ;
    exmplPointArray2nd <<   1.2,  1.5,  1.,  1.5,
                           -2.5, -0.5, 0.5, -0.5;
    cout << exmplPointArray2nd << "\n\n";

    // To access the i-th point one needs to use the method col(point_number)
    cout << "Example 2.3 \n";
    cout << exmplPointArray.col(2) << "\n\n";

    // It is possible to assign a given point within an array of points
    // Notice how one needs to first construct a Point to assign
    cout << "Example 2.4 \n";
    exmplPointArray.col(2) = Point<2>({1.3, 4.5});
    cout << exmplPointArray.col(2) << "\n\n";

    // To access a y component (component 1) of the 3rd point
    cout << "Example 2.5 \n";
    cout << exmplPointArray.col(2)(1) << "\n";
    cout << exmplPointArray(1,2) << "\n\n";
    // Please pay attention to the order of the indices
    
    
    // ************************************************************
    // Vectors of Points (Geometry.hpp) >>> DEVELOPERS <<<
    // ************************************************************
    // In case the number of Points that one needs to store is large, it is better to use VectorPoint
    // VectorPoint requires only template arguments: dimension of the space. The number of points can be changed at runtime
    // To assign the dimension one needs to use resize(Eigen::NoChange, >>number of points<< );
    cout << "Example 3.1 \n";
    VectorPoint<2> exmplPointVector ;
    exmplPointVector.resize(Eigen::NoChange, 4);
    exmplPointVector <<   1.2,  1.5,  1.,  1.5,
                        -2.5, -0.5, 0.5, 0.75;
    cout << exmplPointVector << "\n\n";

    // It is possible to change the number of points, yet it is a destructive operation: previous data is lost
    cout << "Example 3.2 \n";
    VectorPoint<2> exmplPointVector2 = exmplPointVector;
    exmplPointVector2.resize(Eigen::NoChange, 5);
    cout << exmplPointVector2 << "\n\n";

    // To preserve the data one needs to do (it is obviously slower)
    cout << "Example 3.3 \n";
    exmplPointVector.conservativeResize(Eigen::NoChange, 5);
    cout << exmplPointVector << "\n\n";
    exmplPointVector(0,4) = 0.4; exmplPointVector(1,4) = 0.4;

    // Acessing points and coordinates is done as in the case of an array of points
    
    return 0;
}
