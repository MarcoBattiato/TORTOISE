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
//  Plotter.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//
// By including the header of this class an external object plotter3d is defined.
// It is built to have an interface and usage as similar as possible to std::cout.
//
// Before using the object the user should make sure that the object knows the path to gnuplot
// The user can assign it directly by assigning
// plotter3d.gnuplotcommand = "/usr/local/bin/gnuplot > /dev/null 2>&1";
//
// The user needs to define how plotter3d will plot objects, and it is up to the user to overload the << operator (more about this later)
//
// It is built to be used as
// plotter3d  << mesh << NOLABEL << COLOR("red") << function << LABEL(t_title) << COLOR("red") << PLOT;
// or equivalently
// plotter3d << mesh << NOLABEL << COLOR("red") ;
// plotter3d << function << LABEL(t_title) << COLOR("red") << PLOT;
//
// Notice that each plot can contain multiple surfaces or lines.
// The plot is finalised and executed only when the command PLOT is passed (in a sense similarly to the newline charracter in std::cout,
//  yet here the plot is created when the command PLOT is passed and there is no buffer like in std::cout)
//
// Some format commands are already defined. Format commands are to be specified AFTER the object they refer to
//   NOLABEL           : specifies that the current plotted object should not appear in the legend
//   LABEL("title")    : specifies that the current plotted object should appear in the legend with title "title"
//   COLOR(color)      : specifies that the current plotted object should be plotted with the given color
// Colors are passed either using their names or as an exadecimal string in the format  #AARRGGBB
// AA is the alpha channel (00 -> no transparency, FF -> completely transparent), the other ones represent the amount of red, green and blue
// Example COLOR("#BBFF0000")   transparent red
//   PLOTTITLE (title) : specifies the title of the window
//
// The user can overload the << operator to have custom objects to be plotted
//
// Example of Plotter3D << overloading
//Plotter3D& operator << (Plotter3D& plotter, const Region<2>& region){
//    if(plotter.numberplots++>0) plotter.plotCommand +=" , ";                                                              // Very important line
//    plotter.plotCommand += " '-' using 1:2:(0) with lines ";                                                              // passes GNUPlot command
//    Point<2> p(region.origin);        plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";       // Passes data
//    p+= region.gVec.col(0);           plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
//    p+= region.gVec.col(1);           plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
//    p-= region.gVec.col(0);           plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
//    p-= region.gVec.col(1);           plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
//    plotter.plotData += "e\n";                                                                                            // Concludes data
//    return plotter;                                                                                                       // Returns the plotter
//}


#ifndef Plotter_hpp
#define Plotter_hpp

#include <Geometry/GeometryCore/Geometry.hpp>

#include <Configuration>

#include <cassert>
#include <iostream>
#include <stdio.h>
#include <vector>


namespace Tortoise {

enum  PlotCommands { PLOT, NOLABEL, OPENSTREAM, CLOSESTREAM, LINEMODE, POLYGONMODE };
class COLOR { public: const std::string  color;  COLOR(const std::string& t_color): color(t_color){}; };
class LABEL { public: const std::string  label;  LABEL(const std::string& t_label): label(t_label){}; };
class PLOTTITLE { public: const std::string  title;  PLOTTITLE(const std::string& t_title): title(t_title){}; };
class GNUPLOTCOMMAND{ public: const std::string  command;  GNUPLOTCOMMAND(const std::string& t_command): command(t_command){}; };
class TEXT {
    public:
    const std::string  label;
    const std::vector<double> pos;
    
//    TEXT(const std::string& t_label, const std::vector<double>& t_pos): label(t_label), pos(t_pos){};
//    TEXT(const int t_label,          const std::vector<double>& t_pos): label(std::to_string(t_label)), pos(t_pos){};
//    TEXT(const double t_label,       const std::vector<double>& t_pos): label(std::to_string(t_label)), pos(t_pos){};
    
    TEXT(const std::string& t_label, const Point<1> t_pos): label(t_label), pos({t_pos(0)}){};
    TEXT(const int t_label,          const Point<1> t_pos): label(std::to_string(t_label)), pos({t_pos(0)}){};
    TEXT(const double t_label,       const Point<1> t_pos): label(std::to_string(t_label)), pos({t_pos(0)}){};
    TEXT(const std::string& t_label, const Point<2> t_pos): label(t_label), pos({t_pos(0),t_pos(1)}){};
    TEXT(const int t_label,          const Point<2> t_pos): label(std::to_string(t_label)), pos({t_pos(0),t_pos(1)}){};
    TEXT(const double t_label,       const Point<2> t_pos): label(std::to_string(t_label)), pos({t_pos(0),t_pos(1)}){};
    TEXT(const std::string& t_label, const Point<3> t_pos): label(t_label), pos({t_pos(0),t_pos(1),t_pos(2)}){};
    TEXT(const int t_label,          const Point<3> t_pos): label(std::to_string(t_label)), pos({t_pos(0),t_pos(1),t_pos(2)}){};
    TEXT(const double t_label,       const Point<3> t_pos): label(std::to_string(t_label)), pos({t_pos(0),t_pos(1),t_pos(2)}){};
};
class ASPOINTS { public: const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> vertices; template <typename Derived>ASPOINTS(const Eigen::MatrixBase<Derived>& t_vertices): vertices(t_vertices){}};
class ARROW{public: const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> arrow; const Real thickness; template <typename Derived>ARROW(const Eigen::MatrixBase<Derived>& t_arrow, const Real t_thickness): arrow(t_arrow), thickness(t_thickness){}};
class RANGE {public: const Point<2> plotLimits; template <typename Derived> RANGE(const Eigen::MatrixBase<Derived>& MinMax): plotLimits(MinMax){}; EIGEN_MAKE_ALIGNED_OPERATOR_NEW  };


class Plotter3D {
public:
    std::string         plotFormat;
    std::string         plotCommand;
    std::string         plotData;
    int                 numberplots = 0;
    int                 numberobjects = 0;

    const std::string   gnuplotcommand;
    const bool          interactive;
    bool                polygonModeActive = false;

private:
    int                 currentGnuplotPipe= 0;
    std::vector<FILE*>  gnuplotPipeList;
    bool                streamModeActive = false;

public:

    Plotter3D(const std::string& t_Gnuplotcommand = DEFAULTGNUPLOTCOMMAND, bool t_interactive = true);
    
    ~Plotter3D();
    
    Plotter3D& operator << (PlotCommands command);
    Plotter3D& operator << (COLOR color);
    Plotter3D& operator << (LABEL title);
    Plotter3D& operator << (TEXT text);
    Plotter3D& operator << (PLOTTITLE title);
    Plotter3D& operator << (GNUPLOTCOMMAND command);
    Plotter3D& operator << (ASPOINTS points);
    Plotter3D& operator << (ARROW arrow);
    
    template <typename Derived> Plotter3D& operator << (const Eigen::MatrixBase<Derived>& vertices){
        if( vertices.rows() == 2){
            if (streamModeActive) {
                auto siz = vertices.cols();
                for (int i=0; i<siz; i++){
                    plotData += std::to_string(vertices(0,i)) + " " + std::to_string(vertices(1,i)) + " 0. \n";
                }
                if (polygonModeActive){
                    plotData += std::to_string(vertices(0,0)) + " " + std::to_string(vertices(1,0)) + " 0. \n";
                }
                plotData += "\n\n";
//                if (firstStreamed){ plotData += "\n\n";} else {plotData += "e\n";}
//                firstStreamed = true;
            } else {
                if(numberplots++>0) plotCommand +=" , ";
                plotCommand += " '-' using 1:2:(0) with lines ";
                auto siz = vertices.cols();
                for (int i=0; i<siz; i++){
                    plotData += std::to_string(vertices(0,i)) + " " + std::to_string(vertices(1,i)) + "\n";
                }
                if (polygonModeActive){
                    plotData += std::to_string(vertices(0,0)) + " " + std::to_string(vertices(1,0)) + "\n";
                }
                plotData += "e\n";
            }
        } else if (vertices.rows() == 3){
            if (streamModeActive) {
                auto siz = vertices.cols();
                for (int i=0; i<siz; i++){
                    plotData += std::to_string(vertices(0,i)) + " " + std::to_string(vertices(1,i)) + " " + std::to_string(vertices(2,i)) + "\n";
                }
                if (polygonModeActive){
                    plotData += std::to_string(vertices(0,0)) + " " + std::to_string(vertices(1,0)) + " " + std::to_string(vertices(2,0))+ "\n";
                }
                plotData += "\n\n";
//                if (firstStreamed){ plotData += "\n\n";} else {plotData += "e\n";}
//                firstStreamed = true;
            } else {
                if(numberplots++>0) plotCommand +=" , ";
                plotCommand += " '-' with lines ";
                auto siz = vertices.cols();
                for (int i=0; i<siz; i++){
                    plotData += std::to_string(vertices(0,i)) + " " + std::to_string(vertices(1,i)) + " " + std::to_string(vertices(2,i)) + "\n";
                }
                if (polygonModeActive){
                    plotData += std::to_string(vertices(0,0)) + " " + std::to_string(vertices(1,0)) + " " + std::to_string(vertices(2,0))+ "\n";
                }
                plotData += "e\n";
            }
        }
        return *this;
    }
    
};

class Plotter2D {
public:
    std::string         plotFormat;
    std::string         plotCommand;
    std::string         plotData;
    int                 numberplots = 0;
    int                 numberobjects = 0;

    std::string         gnuplotcommand;
    const bool          interactive;
    bool                polygonModeActive = false;

private:
    int                 currentGnuplotPipe= 0;
    std::vector<FILE*>  gnuplotPipeList;
    bool                streamModeActive = false;

public:

    Plotter2D(const std::string& t_Gnuplotcommand = DEFAULTGNUPLOTCOMMAND, bool t_interactive=true);
    ~Plotter2D();
    
    Plotter2D& operator << (PlotCommands command);
    Plotter2D& operator << (COLOR color);
    Plotter2D& operator << (LABEL title);
    Plotter2D& operator << (TEXT text);
    Plotter2D& operator << (PLOTTITLE title);
    Plotter2D& operator << (GNUPLOTCOMMAND command);
    Plotter2D& operator << (ASPOINTS points);
    Plotter2D& operator << (ARROW arrow);
    Plotter2D& operator << (RANGE plotRange);
    
    template <typename Derived> Plotter2D& operator << (const Eigen::MatrixBase<Derived>& vertices){
        if( vertices.rows() == 1){
            if (streamModeActive) {
                auto siz = vertices.cols();
                for (int i=0; i<siz; i++){
                    plotData += std::to_string(vertices(0,i))  + " 0. \n";
                }
                if (polygonModeActive){
                    plotData += std::to_string(vertices(0,0))  + " 0. \n";
                }
                plotData += "\n\n";
//                if (firstStreamed){ plotData += "\n\n";} else {plotData += "e\n";}
//                firstStreamed = true;
            } else {
                if(numberplots++>0) plotCommand +=" , ";
                plotCommand += " '-' using 1:(0) with lines ";
                auto siz = vertices.cols();
                for (int i=0; i<siz; i++){
                    plotData += std::to_string(vertices(0,i))  + "\n";
                }
                if (polygonModeActive){
                    plotData += std::to_string(vertices(0,0))  + "\n";
                }
                plotData += "e\n";
            }
        } else if (vertices.rows() == 2){
            if (streamModeActive) {
                auto siz = vertices.cols();
                for (int i=0; i<siz; i++){
                    plotData += std::to_string(vertices(0,i)) + " " + std::to_string(vertices(1,i)) + "\n";
                }
                if (polygonModeActive){
                    plotData += std::to_string(vertices(0,0)) + " " + std::to_string(vertices(1,0)) + "\n";
                }
                plotData += "\n\n";
//                if (firstStreamed){ plotData += "\n\n";} else {plotData += "e\n";}
//                firstStreamed = true;
            } else {
                if(numberplots++>0) plotCommand +=" , ";
                plotCommand += " '-' with lines ";
                auto siz = vertices.cols();
                for (int i=0; i<siz; i++){
                    plotData += std::to_string(vertices(0,i)) + " " + std::to_string(vertices(1,i)) + "\n";
                }
                if (polygonModeActive){
                    plotData += std::to_string(vertices(0,0)) + " " + std::to_string(vertices(1,0)) + "\n";
                }
                plotData += "e\n";
            }
        }
        return *this;
    }
    
};

extern Plotter3D plotter3d;
extern Plotter2D plotter2d;

} // namespace Tortoise

#endif /* Plotter_hpp */
