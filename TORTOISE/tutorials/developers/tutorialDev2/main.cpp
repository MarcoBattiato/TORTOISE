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
// TUTORIAL 2: Plotting in TORTOISE
//
// Some parts of the tutorial are marked with >>> DEVELOPERS <<<
// That means that the final user most probably will not need to know about those functionalities


// #define NDEBUG            // Uncomment to deactivate debug mode

#include "DiscretizationCore.hpp"
#include <iostream>

using namespace Tortoise;
using std::cout;


int main(int argc, const char * argv[]) {
    
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
    // Plotter3d (Plotter.hpp)
    // ************************************************************
    
    // Before using the object the user should make sure that the object knows the path to GNUPlot
    // The user can assign it directly by assigning
    // plotter3d.gnuplotcommand = "/usr/local/bin/gnuplot > /dev/null 2>&1";
    // The above command works on my Mac.
    // First of all, you need to select the folder wehre your GNUPlot is
    // On Windows it might look something like
    // plotter3d.gnuplotcommand = "C:\\Users\\admin\\gnuplot\\bin\\gnuplot.exe";

    plotter3d.gnuplotcommand = "/usr/local/bin/gnuplot > /dev/null 2>&1";
    
    // It is however possible to directly modify this parameter in Includes/Configuration.hpp
    
    // To plot something, it is necessary to send it to the plotter3d object with the << operator
    
    // ************************************************************
    // Plotting Lines >>> DEVELOPERS <<<
    // ************************************************************

    // To plot a line let us create a vector of Points (an array of points would work in the same way).
    
    VectorPoint<3> exmplPointVector3D { {1.2,  1.5,  1.,  1.5},
                                        {-2.5, -0.5, 0.5, -0.5},
                                        {0., 1.,  2., 3.}};
    // When the vector of points is sent to plotter3d via the << operator, plotter3d reads the cartesian components and
    // starts producing the input for GNUPlot
    plotter3d << exmplPointVector3D ;
    // Notice however that that does not trigger a plot. plotter3d will still be waiting for more instructions.
    // We will see later how to tell plotter3d to collect everything, send the data to GNUPlot and produce a plot.
    
    // ************************************************************
    // Labels >>> DEVELOPERS <<<
    // ************************************************************

    // Each graphical object that is send to plotter3d will automatically create an entry in the legend. If we do not explicitly
    // assign a label, the legend will show "-". We can assign a label to the graphical pbject that has been entered last by
    plotter3d << LABEL("line 1");
    
    // What can we do if we do not want the last graphical object to appear in the legend (for whatever reason)?
    // We can use the line
    // plotter3d << NOLABEL;
    
    // ************************************************************
    // Color >>> DEVELOPERS <<<
    // ************************************************************

    // If we do not specify the color of an entry, GNUPlot will decide it. We can assign a color with
    plotter3d << COLOR("red") ;
    // Notice again that also this formatting flag refers to the last entered object. Or in other terms, we can pass formatting commands
    // after the related object.
    
    // ************************************************************
    // In a single line >>> DEVELOPERS <<<
    // ************************************************************

    // We can add another line.
    VectorPoint<3> exmplPointVector3D_other {{0.2,  0.5,  0.,  0.5, 0.3},{-2.5, -0.5, 0.5, -0.5, 1.2},{0., 1.,  2., 3., 3.}};
    plotter3d << exmplPointVector3D_other << LABEL("line 2") << COLOR("blue");
    
    // Like this we are adding a second line to the plot.
    // plotter3d, noticing that a new graphical onbject is being passed, will finilise the GNUPlot input for the previous graphical object
    // and get ready for the new one.
    // Notice that the formatting flags do not apply to following objects. Here we need to specify them again.
    
    // Also notice how all the commands can be concatenated in a single line (similarly to std::cout)
   
    // Notice plotter3d is still not producing the plot, but waiting for further input.
    
    // ************************************************************
    // 2D points in 3D plots >>> DEVELOPERS <<<
    // ************************************************************

    // If only 2D points are sent, they are plotted at z=0
    VectorPoint<2> exmplPointVector ({{1.2,  1.5,  1.,  1.5},{-2.5, -0.5, 0.5, 0.75}});
    plotter3d << exmplPointVector << NOLABEL << COLOR("green") ;
    
    // Notice that with NOLABEL, the corresponding line will not be added to the legend
    
    // ************************************************************
    // Plot title and producing a plot >>> DEVELOPERS <<<
    // ************************************************************

    // Let us suppose that we are done and we want to wrap everything a produce the plot.
    // Yet, before doing that we give a title to the plot
    
    plotter3d << PLOTTITLE("Example 4.1: First plot") << PLOT;
    
    // When the command PLOT is passed, plotter3d now finalises all the input data required for GNUPlot,
    // launches through the operative system GNUPlot, sends all the data, and then forgets everything to be ready for a new plot
    
    
    // ************************************************************
    // More graphical primitives >>> DEVELOPERS <<<
    // ************************************************************

    // We have seen so far only lines produced by linking a series of points with straight lines.
    // However we can plot more types of primitives
    
    VectorPoint<3> segment ({{1.2,  1.5},{-2.5, -0.5},{-2.5, -0.5}});
    plotter3d << segment << NOLABEL << COLOR("green") ;
    
    // Let us add some text
    // The first argument is a string representing the text to display, and the second the point at which to draw the text
    plotter3d << TEXT("MY TEXT 1", Point<3>({1.2,  1.5, 0.3})) << TEXT("MY TEXT 2", segment.col(1).eval());
    
    
    // We also add two arrows
    // The first argument is a vector or an array of 2 Points 3D and gives the tail and head position of the arrow
    // The second argument is a real number giving the thickness of the arrow
    VectorPoint<3> myArrow ({{1.3,  1.3},{-2.2, -0.2},{-2.2, -0.2}});
    plotter3d << ARROW(myArrow, 1.5);

    // Notice that neither the text or the arrows produce entries in the legend


    // We can also add a series of points. This is done by sending vector of points,
    // but specifying that they must be plotted as points
    plotter3d << ASPOINTS(exmplPointVector3D_other) << LABEL("Some points") << COLOR("red");
    plotter3d << PLOTTITLE("Example 4.2: Other graphical primitives") << PLOT;
    
    // ************************************************************
    // LINEMODE vs POLYGONMODE >>> DEVELOPERS <<<
    // ************************************************************

    // plotter3d has some modes to address some commonly used cases.
    // LINEMODE vs POLYGONMODE
    // The default mode is LINEMODE. It means that the line starts at the first point and ends at the last
    // Howevr often we would like to close the line (like a polygon)
    // POLYGONMODE simply plots a further line that goes from the last point to the first.
    // The mode must be selected before the data is sent in
    VectorPoint<3> exmplPointVector3D_other2 ({{0.2,  0.7,  0.3,  0.6, 1.3},{-1.5, -0.7, 0.1, 0.3, 1.2},{0.3,   1.2,  1., 1., 1.4}});
    plotter3d << exmplPointVector3D <<  LABEL("LINEMODE active by default") << COLOR("red");
    plotter3d << POLYGONMODE << exmplPointVector3D_other <<  LABEL("POLYGONMODE line") << COLOR("blue");
    plotter3d << LINEMODE << exmplPointVector3D_other2 <<  LABEL("Back to LINEMODE") << COLOR("green");
    plotter3d << PLOTTITLE("Example 4.3: LINEMODE vs POLYGONMODE") << PLOT;

    
    // ************************************************************
    // Stream mode >>> DEVELOPERS <<<
    // ************************************************************

    // Another very useful mode is stream mode
    // It allows to send several lines or polygons and have them plot under the same line name
    // For instance it is very useful to draw several triangles without having a legend line for each one of them
    // The default input mode is the non-stream mode. It means that every dataset sent creates a different line
    plotter3d << OPENSTREAM << POLYGONMODE << exmplPointVector3D << LINEMODE << exmplPointVector3D_other << exmplPointVector3D_other2;
    plotter3d << LABEL("Single legend line") << COLOR("red") << CLOSESTREAM;
    plotter3d << PLOTTITLE("Example 4.4: Stream mode") << PLOT;
    // Notice how one has to open and close the stream. Also notice how one can change from POLYGONMODE to LINEMODE during a stream
    
    
    // ************************************************************
    // Plotter2d (Plotter.hpp) >>> DEVELOPERS <<<
    // ************************************************************
    
    // There is also a plotter for 2D plots
    // It works in the same way, but before using it we need to remember to specify the command to starg GNUPlot

    plotter2d.gnuplotcommand = "/usr/local/bin/gnuplot > /dev/null 2>&1";

    
    
    return 0;
}
