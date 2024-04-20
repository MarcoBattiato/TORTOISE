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
//  Plotter.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//

#include <Geometry/Plotter/Plotter.hpp>

#include <iostream>

namespace Tortoise {

//======================================
// PLOTTER 3D
//======================================

Plotter3D::Plotter3D(const std::string& t_Gnuplotcommand, bool t_interactive): gnuplotcommand(t_Gnuplotcommand), interactive(t_interactive) {
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    gnuplotPipeList.emplace_back(_popen(gnuplotcommand.c_str(), "w"));  // Open a pipe to gnuplot for the first plot
#else
    gnuplotPipeList.emplace_back(popen(gnuplotcommand.c_str(), "w"));  // Open a pipe to gnuplot for the first plot
#endif
    plotCommand = "splot ";
    plotFormat = "set xlabel 'x'\n set ylabel 'y'\n set zlabel 'z'\n set xyplane at 0\n set hidden3d nooffset\n";
    plotData = "";
};
Plotter3D::~Plotter3D(){
    if (interactive){
        std::cout << "press return/enter to close all 3D windows...";
        getchar();
    }
    for (int i=0; i<gnuplotPipeList.size(); ++i){
        fprintf(gnuplotPipeList[i],"exit \n");   // exit gnuplot
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
        _pclose(gnuplotPipeList[i]);
#else
        pclose(gnuplotPipeList[i]);
#endif
    }
}

Plotter3D& Plotter3D::operator << (PlotCommands command) {
    switch (command) {
        case PLOT:
            if (gnuplotPipeList[currentGnuplotPipe]) {
                plotFormat += "\n";
                plotCommand += "\n";
                fprintf(gnuplotPipeList[currentGnuplotPipe],"%s", plotFormat.c_str());
                fprintf(gnuplotPipeList[currentGnuplotPipe],"%s", plotCommand.c_str());
                fprintf(gnuplotPipeList[currentGnuplotPipe],"%s", plotData.c_str());
                fflush(gnuplotPipeList[currentGnuplotPipe]);
                ++currentGnuplotPipe;
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
                gnuplotPipeList.emplace_back(_popen(gnuplotcommand.c_str(), "w"));  // Open a pipe to gnuplot
#else
                gnuplotPipeList.emplace_back(popen(gnuplotcommand.c_str(), "w"));  // Open a pipe to gnuplot
#endif
                plotCommand = "splot ";
                plotFormat = "set xlabel 'x'\nset ylabel 'y'\nset zlabel 'z'\nset xyplane at 0\n set hidden3d nooffset\n";
                plotData = "";
                numberplots= 0;
                numberobjects =0;
            }
            break;
        case NOLABEL:
            if (numberplots>0 && !streamModeActive) plotCommand += " not ";
            break;
        case LINEMODE:
            polygonModeActive = false;
            break;
        case POLYGONMODE:
            polygonModeActive = true;
            break;
        case OPENSTREAM:
            if(numberplots++>0) plotCommand +=" , ";
//            plotCommand += " for [i=0:*] '-' index i with lines ";
            plotCommand += " '-' with lines ";
            streamModeActive = true;
            break;
        case CLOSESTREAM:
            plotData += "e\n";
            streamModeActive = false; // firstStreamed = false;
            break;
        default:
            break;
    }
    return *this;
}
Plotter3D& Plotter3D::operator << (COLOR color) {if (numberplots>0) plotCommand += " lt rgb '" + color.color + "' "; return *this;};
Plotter3D& Plotter3D::operator << (LABEL label) {
    if (numberplots>0) {
//        if (streamModeActive) {
//            plotFormat += "myTitle(i) = i==0 ? 'only one title' : ''\n";
//            plotCommand += " title myTitle(i)";
//        } else {
            plotCommand += " title '" + label.label + "' ";
//        }
    }
    return *this;};
Plotter3D& Plotter3D::operator << (PLOTTITLE title) {plotFormat += "set title '" + title.title + "'\n"; return *this;};
Plotter3D& Plotter3D::operator << (GNUPLOTCOMMAND command) {plotFormat += command.command; return *this;}
Plotter3D& Plotter3D::operator << (TEXT text) {
    plotFormat += " set label '"+ text.label + "' at ";
    for (int i=0; i<text.pos.size()-1; ++i) {plotFormat += std::to_string(text.pos[i]) + ",";}
    plotFormat += std::to_string(text.pos[text.pos.size()-1]) + "\n";
    return *this;}
Plotter3D& Plotter3D::operator << (ASPOINTS points) {
    if( points.vertices.rows() == 2){
        if (streamModeActive) {
            auto siz = points.vertices.cols();
            for (int i=0; i<siz; i++){
                plotData += std::to_string(points.vertices(0,i)) + " " + std::to_string(points.vertices(1,i)) + " 0. \n";
            }
            if (polygonModeActive){
                plotData += std::to_string(points.vertices(0,0)) + " " + std::to_string(points.vertices(1,0)) + " 0. \n";
            }
            plotData += "\n\n";
//                if (firstStreamed){ plotData += "\n\n";} else {plotData += "e\n";}
//                firstStreamed = true;
        } else {
            if(numberplots++>0) plotCommand +=" , ";
            plotCommand += " '-' using 1:2:(0) with points ";
            auto siz = points.vertices.cols();
            for (int i=0; i<siz; i++){
                plotData += std::to_string(points.vertices(0,i)) + " " + std::to_string(points.vertices(1,i)) + "\n";
            }
            if (polygonModeActive){
                plotData += std::to_string(points.vertices(0,0)) + " " + std::to_string(points.vertices(1,0)) + "\n";
            }
            plotData += "e\n";
        }
    } else if (points.vertices.rows() == 3){
        if (streamModeActive) {
            auto siz = points.vertices.cols();
            for (int i=0; i<siz; i++){
                plotData += std::to_string(points.vertices(0,i)) + " " + std::to_string(points.vertices(1,i)) + " " + std::to_string(points.vertices(2,i)) + "\n";
            }
            if (polygonModeActive){
                plotData += std::to_string(points.vertices(0,0)) + " " + std::to_string(points.vertices(1,0)) + " " + std::to_string(points.vertices(2,0))+ "\n";
            }
            plotData += "\n\n";
//                if (firstStreamed){ plotData += "\n\n";} else {plotData += "e\n";}
//                firstStreamed = true;
        } else {
            if(numberplots++>0) plotCommand +=" , ";
            plotCommand += " '-' with points ";
            auto siz = points.vertices.cols();
            for (int i=0; i<siz; i++){
                plotData += std::to_string(points.vertices(0,i)) + " " + std::to_string(points.vertices(1,i)) + " " + std::to_string(points.vertices(2,i)) + "\n";
            }
            if (polygonModeActive){
                plotData += std::to_string(points.vertices(0,0)) + " " + std::to_string(points.vertices(1,0)) + " " + std::to_string(points.vertices(2,0))+ "\n";
            }
            plotData += "e\n";
        }
    }
    return *this;
}
Plotter3D& Plotter3D::operator << (ARROW arrow) {
    plotFormat += " set arrow from "+ std::to_string(arrow.arrow(0,0)) + ", "+ std::to_string(arrow.arrow(1,0))+ ", "+ std::to_string(arrow.arrow(2,0)) + " to " + std::to_string(arrow.arrow(0,1)) + ", " + std::to_string(arrow.arrow(1,1))+ ", " + std::to_string(arrow.arrow(2,1)) + " lw " + std::to_string(arrow.thickness) + "\n";
    return *this;
}

Plotter3D plotter3d = Plotter3D();




//======================================
// PLOTTER 2D
//======================================


Plotter2D::Plotter2D(const std::string& t_Gnuplotcommand, bool t_interactive): gnuplotcommand(t_Gnuplotcommand), interactive(t_interactive) {
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    gnuplotPipeList.emplace_back(_popen(gnuplotcommand.c_str(), "w"));  // Open a pipe to gnuplot for the first plot
#else
    gnuplotPipeList.emplace_back(popen(gnuplotcommand.c_str(), "w"));  // Open a pipe to gnuplot for the first plot
#endif
    plotCommand = "plot ";
    plotFormat = "set xlabel 'x'\n set ylabel 'y'\n\n";
    plotData = "";
};
Plotter2D::~Plotter2D(){
    if (interactive){
        std::cout << "press return/enter to close all 2D windows...";
        getchar();
    }
    for (int i=0; i<gnuplotPipeList.size(); ++i){
        fprintf(gnuplotPipeList[i],"exit \n");   // exit gnuplot
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
        _pclose(gnuplotPipeList[i]);
#else
        pclose(gnuplotPipeList[i]);
#endif
    }
}

Plotter2D& Plotter2D::operator << (PlotCommands command) {
    switch (command) {
        case PLOT:
            if (gnuplotPipeList[currentGnuplotPipe]) {
                plotFormat += "\n";
                plotCommand += "\n";
                fprintf(gnuplotPipeList[currentGnuplotPipe],"%s", plotFormat.c_str());
                fprintf(gnuplotPipeList[currentGnuplotPipe],"%s", plotCommand.c_str());
                fprintf(gnuplotPipeList[currentGnuplotPipe],"%s", plotData.c_str());
                fflush(gnuplotPipeList[currentGnuplotPipe]);
                ++currentGnuplotPipe;
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
                gnuplotPipeList.emplace_back(_popen(gnuplotcommand.c_str(), "w"));  // Open a pipe to gnuplot
#else
                gnuplotPipeList.emplace_back(popen(gnuplotcommand.c_str(), "w"));  // Open a pipe to gnuplot
#endif
                plotCommand = "plot ";
                plotFormat = "set xlabel 'x'\nset ylabel 'y'\n";
                plotData = "";
                numberplots= 0;
                numberobjects =0;
            }
            break;
        case NOLABEL:
            if (numberplots>0 && !streamModeActive) plotCommand += " not ";
            break;
        case LINEMODE:
            polygonModeActive = false;
            break;
        case POLYGONMODE:
            polygonModeActive = true;
            break;
        case OPENSTREAM:
            if(numberplots++>0) plotCommand +=" , ";
//            plotCommand += " for [i=0:*] '-' index i with lines ";
            plotCommand += " '-' with lines ";
            streamModeActive = true;
            break;
        case CLOSESTREAM:
            plotData += "e\n";
            streamModeActive = false; // firstStreamed = false;
            break;
        default:
            break;
    }
    return *this;
}
Plotter2D& Plotter2D::operator << (COLOR color) {if (numberplots>0) plotCommand += " lt rgb '" + color.color + "' "; return *this;};
Plotter2D& Plotter2D::operator << (LABEL label) {
    if (numberplots>0) {
//        if (streamModeActive) {
//            plotFormat += "myTitle(i) = i==0 ? 'only one title' : ''\n";
//            plotCommand += " title myTitle(i)";
//        } else {
            plotCommand += " title '" + label.label + "' ";
//        }
    }
    return *this;};
Plotter2D& Plotter2D::operator << (PLOTTITLE title) {plotFormat += "set title '" + title.title + "'\n"; return *this;};
Plotter2D& Plotter2D::operator << (GNUPLOTCOMMAND command) {plotFormat += command.command; return *this;}
Plotter2D& Plotter2D::operator << (TEXT text) {
    plotFormat += " set label '"+ text.label + "' at ";
    for (int i=0; i<text.pos.size()-1; ++i) {plotFormat += std::to_string(text.pos[i]) + ",";}
    plotFormat += std::to_string(text.pos[text.pos.size()-1]) + "\n";
    return *this;}
Plotter2D& Plotter2D::operator << (ASPOINTS points) {
    if( points.vertices.rows() == 1){
        if (streamModeActive) {
            auto siz = points.vertices.cols();
            for (int i=0; i<siz; i++){
                plotData += std::to_string(points.vertices(0,i))  + " 0. \n";
            }
            if (polygonModeActive){
                plotData += std::to_string(points.vertices(0,0))  + " 0. \n";
            }
            plotData += "\n\n";
//                if (firstStreamed){ plotData += "\n\n";} else {plotData += "e\n";}
//                firstStreamed = true;
        } else {
            if(numberplots++>0) plotCommand +=" , ";
            plotCommand += " '-' using 1:(0) with points ";
            auto siz = points.vertices.cols();
            for (int i=0; i<siz; i++){
                plotData += std::to_string(points.vertices(0,i))  + "\n";
            }
            if (polygonModeActive){
                plotData += std::to_string(points.vertices(0,0))  + "\n";
            }
            plotData += "e\n";
        }
    } else if (points.vertices.rows() == 2){
        if (streamModeActive) {
            auto siz = points.vertices.cols();
            for (int i=0; i<siz; i++){
                plotData += std::to_string(points.vertices(0,i)) + " " + std::to_string(points.vertices(1,i)) + "\n";
            }
            if (polygonModeActive){
                plotData += std::to_string(points.vertices(0,0)) + " " + std::to_string(points.vertices(1,0)) + "\n";
            }
            plotData += "\n\n";
//                if (firstStreamed){ plotData += "\n\n";} else {plotData += "e\n";}
//                firstStreamed = true;
        } else {
            if(numberplots++>0) plotCommand +=" , ";
            plotCommand += " '-' with points ";
            auto siz = points.vertices.cols();
            for (int i=0; i<siz; i++){
                plotData += std::to_string(points.vertices(0,i)) + " " + std::to_string(points.vertices(1,i)) + "\n";
            }
            if (polygonModeActive){
                plotData += std::to_string(points.vertices(0,0)) + " " + std::to_string(points.vertices(1,0)) + "\n";
            }
            plotData += "e\n";
        }
    }
    return *this;
}
Plotter2D& Plotter2D::operator << (ARROW arrow) {
    plotFormat += " set arrow from "+ std::to_string(arrow.arrow(0,0)) + ", "+ std::to_string(arrow.arrow(1,0)) + " to " + std::to_string(arrow.arrow(0,1)) + ", " + std::to_string(arrow.arrow(1,1)) + " lw " + std::to_string(arrow.thickness) + "\n";
    return *this;
}
Plotter2D& Plotter2D::operator << (RANGE plotRange){
// Author: Indrajit Wadgaonkar
    plotFormat += "set yrange ["+ std::to_string(plotRange.plotLimits(0))+":" + std::to_string(plotRange.plotLimits(1))+"]\n";
    return *this;
}

Plotter2D plotter2d = Plotter2D();

} // namespace Tortoise
