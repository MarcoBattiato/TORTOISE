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
// TUTORIAL 1: Introduction to TORTOISE



// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

// ******************************
// How to run
// ******************************

// TORTOISE requires -O3 optimisation (as a large part of the code is built to trigger certain optimisations)
// It also requires to be compiled with c++20, since TORTOISE uses several features introduced in c++20.



// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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

// We include TORTOISE with

#include <TORTOISE>

// Please notice that you are expected to have added the directory of TORTOISE header files
// to the list of header search paths

// All the classes in TORTOISE are within the namespace Totroise. Therefore all classes and objects within TORTOISE need to be called with Tortoise::name_of_class
// That will be very cumbersome.
// It is usually a bad idea to use the keyword "using" with an entire namespace, unless that namespace
// is the key one. Therfore in this case we will do

using namespace Tortoise;

// There is little risk involved in this since the names within TORTOISE tend to be rather specific

#include <iostream>
using std::cout;

int main(int argc, const char * argv[]) {
    
    // First we learn few fundamental functionalities offered by "Geometry.hpp"
    
    // ************************************************************
    // The type alias: real (Geometry.hpp)
    // ************************************************************

    // The whole library is built in a way that it is easy to change the precision of the used real numbers.
    // Instead of explicitly using float or double, etc, the user should use the type alias: real
    // In the version this tutorial was written for, real = double
    // However in a later version the precision might be decreased or increased.

    // Do not use
    // double myValue = 1.35;
    // Use instead
    real myValue = 1.35;

    // This will ensure that all the numbers within TORTOISE have the same precision
    
    // If for whatever reason you want to change the type aliased by real
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
    // Plotter3d (Plotter.hpp)
    // ************************************************************

    // Let us now have a look at how graphics is produced within TORTOISE
    // NOTE: The plotting part has been implemented rushly so it is not performant
    //       Feel free to suggest improvements
    
    // By including Plotter.hpp an external object plotter3d is defined.
    // It is built to have an interface and usage as similar as possible to std::cout.
    
    // The plotter3D object is nothing else than a frontend to GNUPlot, built to simplify the transfer of data from
    // TORTOISE objects to GNUPlot. The plotter3D object takes care of producing the appropriate input for GNUPlot
    
    // Before doing that there is one important this to do.
    // Since plotter3D needs to launch GNUPlot, it must interact with the operative system.
    // Programming in c++ in Linux/Mac and Windows have some important differences, especially when the c++ program has to
    // interact with the operative system.
    
    // Therefore you need to change some line in TORTOISE depending on your opertive system.
    // Before continuing with this tutorial,
    // go to Geometry/Plotter/Plotter.hpp   Line 90 and follow the instructions.
    // Then save that file.
    
    
    
    // ************************************************************
    // Plotter2d and Plotter3d (Plotter.hpp)
    // ************************************************************
    
    // Before using the objects the user should make sure that the object knows the path to GNUPlot
    // The user can assign it directly by assigning
    
    // plotter3d.gnuplotcommand = "/usr/local/bin/gnuplot > /dev/null 2>&1";
    // plotter2d.gnuplotcommand = "/usr/local/bin/gnuplot > /dev/null 2>&1";
    
    // where you need to use the appropriate path to GNUPlot.
    // On Windows it might look something like
    // plotter3d.gnuplotcommand = "C:\\Users\\admin\\gnuplot\\bin\\gnuplot.exe";
    
    // The remaining aprt of the command gives instruction to GNUPlot about where the send the output.
    
    // There is however a more convenient way f doing that.
    // It is  possible to directly modify this parameter in Includes/Configuration.hpp
    
    // Once the plotters are configured, to plot something, it is necessary to send
    // graphica objects to the plotter3d object with the << operator

    // We will see below how to use them
    
    
    // ************************************************************
    // Discretised functions
    // ************************************************************
    
    // Below we will learn how to discretise and use discretised functions within TORTOISE.
    
    // It is not possible in TORTOISE to discretise function over an infinite domain. They can only be discretised over
    // finite domains. For reasons that will be addressed in a later tutorial, the only domains that can be used in the
    // current version of TORTOISE are lines in 1D, parallelograms in 2D and parallelepypeds in 3D.
    // The domain of the function is then meshed with a regular mesh: segments in 1D, triangles in 2D and tetrahedra in 3D
    
    // Keep in mind that we will use discretised functions, among other cases, to describe bands. Each band is a function
    // defined over the whole brillouin zone. That however does not mean that we are forced to discretise the band across the
    // whole Brillouin zone. Some parts of the band wil be at energies that are too high to be participating to the dynamics we
    // want to study. For that reason TORTOISE gives flexibility in choosing the domain within the brillouin zone.
    // Therefore, we can have several functions, defining the dispersion of bands, that are discretised over different domains.
    
    // Let us start by learning how to build meshes. However before constructing a mesh, we need to build a region.
    
    // ******************************
    // Regions (Geometry/Structured/DirectSpace/Region.hpp)
    // ******************************
    // Regions are the simplest gemetric shapes that are handled by the library.
    // They are simply NDim-dimensional parallelograms.
    // From a less abstract point of view, regions are used to store the Brillouin zone of a material.
    // We will build a region that represent the Brillouin zone and then extract domains from the region
    // over which we will define the dispersion of bands (but we will see this later).
    
    // Let us see how we build a region

    // Region needs one template argument: the number of dimensions
    // They also need several arguments for their constructor:
    // -> origin
    // -> the sides
    
    // Let us construct what we will send to the constructor
    Point<2>  origin({0.0,-0.1});
    Point<2>  sideRegion0({1.0,0.1}), sideRegion1({.1,1.1});
    Region<2> region(origin, {sideRegion0, sideRegion1} );
    // You can also define a region simply as
    // Region<2> region({0.0,-0.1}, {{1.0,0.1}, {.1,1.1}} );
    // But I personally prefer the first way as it makes it clear what each vector means.


    // It is possible to plot the region
    // The method below uses plotter3d
    region.plot("Example 1.1: Region plot");
    // The function above simply calls the plotter3D and finalises the plot: you cannot add anything else. Go and see how it is implemented

    // If one needs to add further elements to the plot, it is necessary to use
    plotter3d << region << COLOR("red") << LABEL("region");
    plotter3d << PLOTTITLE("Example 1.2: Region Plot with other lines") << PLOT;
    // Have a look at how the method plot and the << operator are implemented in Region.cpp


    // ******************************
    // Mesh (Geometry/Structured/DirectSpace/Mesh.hpp)
    // ******************************
    // We now start building the geometrical objects that are required to discretise the functions
    // The mesh object can be constructed on a part of a region.
    // For reasons that are not relevant here, the part of the region to be meshed must have the sides parallel
    // to the sides of the region.
    // The mesh requires a template parameter: the number of dimensions
    // The arguments of the constructor are
    // -> region over which to build the mesh
    // -> origin of the mesh as fraction of the region
    // -> sides of the mesh as fraction of the region's sides
    // -> splits per side
    Mesh<2> mesh(region, Point<2>({0.1,0.2}), Point<2>({0.6,0.5}),CartIndex<2>({4,3}));
    mesh.plot("Example 2.1: Mesh plot");
    // Similarly to regions, meshes can be directly sent to plotter3d for a higher personalisation of the plot
    // Have a look at how the method plot and the << operator are implemented in Mesh.cpp

    // The mesh has several characteristics
    // Let us access the simplest ones
    // The other ones can be looked up on the header file Geometry/Structured/DirectSpace/Mesh.hpp
    cout << "Example 2.2 \n";
    cout << mesh.origin << "\n";

    // We can print a variety of information about the mesh
    cout << "Example 2.3 \n";
    cout << mesh;


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
    Function<2> funct0(mesh);         // Creates a function object on the Mesh mesh, and initialises it to constant 0.0
    Function<2> funct1(mesh, 3.4);    // Creates a function object on the Mesh mesh, and initialises it to constant 3.4
    
    // We can also initialise Function<n> objects by providing an analytic expression to discretise.
    // Notice that the analytic expression is used only for the construction of the discretised object and then discarded.
    // We can construct using lambda expressions
    Function<2> funct2(mesh, [](Point<2> P){return P(0)+P(1)*P(1);});  // Creates a function object on the Mesh mesh, and initialises it to the discretised version of x+y^2
    // But also TORTOISE provides functionalities to pass nicer looking expressions using x, y, and z,
    // or kx, ky, kz and k (with the latter intended as the vector of coordinates)
    Function<2> funct3(mesh, x + pow(y,2));
    // Please notice that it is not advisable to use the objects x, y, z, kx, ky, kz and k beyond the usage suggested below, unless the user is not fully aware of
    // what these objects are. They are actually functors built with expression templates, to make their evaluation as efficient as possible. In particular expression
    // templates require some understanding, to be used properly beyond what is suggested here.
    
    // Notice that, while regions and meshes cannot be modified after being created, functions can be modifed
    
    // There are more advanced constructor
    
    // This construct a new discretised function as the result of some expression of another discretised function
    Function<2> funct4([](real z){return 1./z;}, funct2);
    // Notice how it is not necessary to pass a mesh, since this constructor creates a function object on the same
    // mesh of funct2, and initialises it to 1/funct2(x,y)
    
    // As a function of more variables
    Function<2> funct5([](std::vector<real> z){return 1./z[0] + 3.*z[0]*z[1];}, {funct2, funct3});
    // Creates a function object on the same mesh of funct2, and initialises it to 1/funct2(x,y) + funct2(x,y) * funct3(x,y)
    // Notice that if funct2 and funct3 are defined on different meshes, the constructor will return a runtime error if debug
    // mode is active
    
    // Of course the usual copy constructor is provided
    Function<2> funct6(funct3);    // Makes a new function, copying funct4
    
    //======================
    // Evaluation
    //======================
    // A Function can be evaluated at a point. However, given that the function can be discontinuous at the elements' edges,
    // over there the result depends on which element one is looking at.
    // Therefore there is no method to evaluate a function at a point in global coordinates
    // The evaluation of a function is done by specifying the element and then the local coordinate of the point
    for (auto elem = mesh.elementIterator(); elem.unfinished; ++elem) {
        cout << funct4(elem, Point<2>({0.333,0.333})) << " ";     // Value of the function at the baricentre of every element
    }
    cout << "\n";
    
    for (auto elem = mesh.elementIterator(); elem.unfinished; ++elem) {
        for (auto vert = mesh.vertexInElemIterator(); vert.unfinished; ++vert) {
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
