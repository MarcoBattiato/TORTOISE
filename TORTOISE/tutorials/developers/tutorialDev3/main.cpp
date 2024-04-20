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
// TUTORIAL 3: Regions and Meshes
//
// Some parts of the tutorial are marked with >>> DEVELOPERS <<<
// That means that the final user most probably will not need to know about those functionalities


// #define NDEBUG            // Uncomment to deactivate debug mode

#include "DiscretizationCore.hpp"
#include <iostream>
#include <functional>

using namespace Tortoise;
using std::cout;

int main(int argc, const char * argv[]) {
    
    // Below we will learn how to discretise and use discretised functions within TORTOISE.
    
    // It is not possible in TORTOISE to discretise function over an infinite domain. They can only be discretised over
    // finite domains. For reasons that will be addressed in a later tutorial, the only domains that can be used in
    // TORTOISE are lines in 1D, parallelograms in 2D and parallelepypeds in 3D.
    // The domain of the function is then meshed with a regular mesh: segments in 1D, triangles in 2D and tetrahedra in 3D
    
    // Keep in mind that we will use discretised functions, among other cases, to describe bands. Each band is a function
    // defined over the whole brillouin zone. That however does not mean that we are forced to discretise the band across the
    // whole Brillouin zone. Some parts of the band wil be a energies that are too high to be participating to the dynamics we
    // want to study. For that reason TORTOISE gives flexibility in choosing the domain within the brillouin zone,.
    // Therefore, we can have several functions, defining the dispersion of bands, that are discretised over different domains.
    // However notice that logically they are functions of the same variable
    
    // Let us start by learning how to build meshes
    
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

    // ******************************
    // Mesh >>> DEVELOPERS <<<
    // ******************************

    // Sometimes we might want to extract more information from a mesh, or perhaps create plots which highlight
    // something more than what is already implemented. In this case it is better to understand a bit more how the
    // mesh is structured.

    // The mesh is composed of elements (for instance in 2D elements are triangles). TORTOISE has assigned an order to
    // the elements, however it is not really relevant to try and keep in mind what is the numbering convention, as in
    // all realistic applications I can think of it is not really relevant. Therefore I am not going to explain it.
    // However there can be cases where we want to do something to each element, so we need to iterate through the
    // elements (again the order will often not matter).

    // As a first exercise let us plot the index that TORTOISE assigned to each element for the 2D mesh that we have
    // already constructed
    plotter3d << mesh;   // We send the mesh to the plotter
    for (int i=0; i <mesh.numberElements; ++i){
        plotter3d << TEXT( i, mesh.elemCentre(i));  // We send a text to plot
        // The first argument is a string or a number (which will be converted in a string), while the second the position
        // where to plot the text (sorry font size cannot be controlled in the current version)
    }
    plotter3d << PLOTTITLE("Example 3.1: Element indexing") << PLOT;

    // Let us now suppose that we want to shift the text closer to one of the vertices of each element.
    
    // To do that we need to understand how points within an element are addressed.
    // Points within each element are referenced using relative coordinates. Meaning that a one-to-one correspondence is made
    // with the reference simplex: in 1D the [0,1] segment, in 2D the triangle with vertices [0,0], [1,0], [0,1] and in 3D the
    // tetrahedron with vertices [0,0,0], [1,0,0],  [0,1,0], [0,0,1]

    // Let us make an example in 2D, first
    // We plot the reference triangle and a generic triangle
    ArrayPoint<2, 3>  referTriangle ({{0, 1, 0},{0, 0, 1}}), triangle ({{1.3, .6, 2.7},{0.9, 1.7, 2.2}});

    std::function<void(Plotter2D&)> prepareTwoTrianglesPlot = [&referTriangle, &triangle](Plotter2D&){
        ArrayPoint<2, 4>  figBox ({{-0.5, 3, 3, -0.5},{-0.5, -0.5, 3, 3}});
        ArrayPoint<2, 2>  axisX ({{-0.5, 3},{0, 0}}), axisY ({{0,0},{-0.5, 3}});
        plotter2d << POLYGONMODE << figBox << NOLABEL << COLOR("white") ;
        plotter2d << referTriangle << LABEL("Reference Triangle") << COLOR("green") ;
        plotter2d << triangle << LABEL("Triangle")<< COLOR("blue");
        for (int i=0; i<3; ++i){
            plotter2d << TEXT("VR"+std::to_string(i), (referTriangle.col(i)+Point<2>({0.03,0.05})).eval());
            plotter2d << TEXT("VA"+std::to_string(i), (triangle.col(i)+Point<2>({0.03,0.05})).eval());
        }
        plotter2d << ARROW(axisX, 1.5) << ARROW(axisY, 1.5);
    };
    prepareTwoTrianglesPlot(plotter2d);
    plotter2d << PLOTTITLE("Example 4.1: Relative Coordinates") << PLOT;
    
    // We now highlight the corresponding vertices
    prepareTwoTrianglesPlot(plotter2d);
    for (int i=0; i<3; ++i){
        ArrayPoint<2, 2> arrow;
        arrow.col(0) = referTriangle.col(i); arrow.col(1) = triangle.col(i);
        plotter2d << ARROW(arrow, 1.5);
    }
    plotter2d << PLOTTITLE("Example 4.2: Corresponding vertices") << PLOT;
    
    // With a similar logic we can highlight other corresponding points
    prepareTwoTrianglesPlot(plotter2d);
    ArrayPoint<2, 2> arrow;
    arrow.col(0) = 0.5*(referTriangle.col(0)+referTriangle.col(1));
    arrow.col(1) = 0.5*(triangle.col(0)+triangle.col(1));
    plotter2d << ARROW(arrow, 1.5);
    
    arrow.col(0) = 0.333*(referTriangle.col(0)+referTriangle.col(1)+referTriangle.col(2));
    arrow.col(1) = 0.333*(triangle.col(0)+triangle.col(1)+triangle.col(2));
    plotter2d << ARROW(arrow, 1.5);
    plotter2d << PLOTTITLE("Example 4.3: Corresponding Points") << PLOT;
    
    // For more info see the manual
    
    
    // Summarising, a point within the mesh is identified by the element it belongs to and by its
    // relative coordinates.
    // To obtain the absolute coordinates of a point within a specific element at specific relative
    // coordinates we use
    
    cout << "Example 5.1 \n";
    cout << mesh(5, {0.3,0.5}); // Coordinates of the point in element 5 with relative coordinates {0.3,0.5}
    cout << "\n\n";

    // Notice that meaningful relative coordinates have to be within the reference simplex. For instance in 2D,
    // the relative coordinates {2.3,1.5} or {0.6,0.7} do not make sense as they are outside the reference triangle.
    // TORTOISE does not check whether you send meaningful relative coordinates. It is your responsability.
    
    
    // We now can use what we learned to print the relative coordinates of each vertex of each element
    Mesh<2> mesh2(region, Point<2>({0.1,0.2}), Point<2>({0.6,0.5}),CartIndex<2>({3,2}));
    plotter3d << mesh2 << LABEL("mesh");   // We send the mesh to the plotter
    for (int i=0; i <mesh2.numberElements; ++i){
        plotter3d << TEXT("[0,0]", mesh2(i, {0.1,0.1}));
        plotter3d << TEXT("[1,0]", mesh2(i, {0.8,0.1}));
        plotter3d << TEXT("[0,1]", mesh2(i, {0.1,.8}));
    }
    plotter3d << PLOTTITLE("Example 6.1: Element indexing") << PLOT;
    
    
    
    return 0;
}
