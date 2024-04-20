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
//  Region.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//


#include <Geometry/Structured/DirectSpace/Region.hpp>

#include <cassert>
#include <fstream>

namespace Tortoise {


//=======================================================
// Constructors
//===================

template<int NDim> VectorPoint<NDim> calculate_gVec(const Vect_Points<NDim> &t_gVec){
    VectorPoint<NDim> toreturn(NDim,NDim);
    for (int i = 0; i < NDim ; ++i){
        toreturn.col(i)=t_gVec[i].transpose();
    }
    return toreturn;
}

template<int NDim> Region<NDim>::Region(const Point<NDim>& t_origin, const Vect_Points<NDim>& t_gVec): origin(t_origin), gVec(calculate_gVec<NDim>(t_gVec)) {
    // ******* Assertions
    assert(t_gVec.size()==NDim); // Check if the user passed vectors with correct dimensions
}
template<int NDim> Region<NDim>::Region(const Real t_origin, const Real t_gVec) requires (NDim==1): origin(Point<NDim>::Constant(t_origin)), gVec(ArrayPoint<NDim,NDim>::Constant(t_gVec)) {};

template<int NDim> Point<NDim> getOriginFromFile(const std::string &FileName){
    Point<NDim> origin;
    std::string linedata;
    std::ifstream inputfile;
    inputfile.open(FileName);
    std::getline(inputfile, linedata);
    std::stringstream linestream(linedata);
    for (int i=0; i<NDim; ++i){
        Real coord;
        linestream >> coord;
        origin(i) =coord;
    }
    inputfile.close();
    return origin;
}
template<int NDim> ArrayPoint<NDim,NDim> getgVecsFromFile(const std::string &FileName){
    ArrayPoint<NDim,NDim> gVec;
    std::string linedata;
    std::ifstream inputfile;
    inputfile.open(FileName);
    std::getline(inputfile, linedata);  // Skips first line
    for (int j=0; j<NDim; ++j){
        std::getline(inputfile, linedata);
        std::stringstream linestream(linedata);
        for (int i=0; i<NDim; ++i){
            Real coord;
            linestream >> coord;
            gVec(i,j) = coord;
        }
    }
    inputfile.close();
    return gVec;
}
template<int NDim> Region<NDim>::Region(const std::string &FileName): origin(getOriginFromFile<NDim>(FileName)), gVec(getgVecsFromFile<NDim>(FileName)){}

//=======================================================
// Geometry
//===================
template<int NDim> Point<NDim> Region<NDim>::toRelative(const Point<NDim>& absCoord) const {
    return gVec.inverse()*(absCoord-origin);
}
template<int NDim> Point<NDim> Region<NDim>::toAbsolute(const Point<NDim>& relCoord) const {
    return gVec*relCoord+origin;
}


//=======================================================
// Comparison
//===================
template<int NDim> bool Region<NDim>::equivalentTo(const Region<NDim>& mesh2) const{
    bool result;
    result = (origin == mesh2.origin);
    result = result && (gVec == mesh2.gVec);
    return result;
};
template<int NDim> bool Region<NDim>::equivalentTo(const std::string &FileName) const{
    Region<NDim> regionTemp(FileName);
    return equivalentTo(regionTemp);
};


//=======================================================
// I/O
//===================

//template<int NDim> std::ostream &operator<<(std::ostream &os, Region<NDim> const& region) {
//    os <<  "********************************************\n";
//    os <<  " Region: \n";
//    os <<  "origin  : " + to_string(region.origin.transpose()) + "\n";
//    for (int i=0; i<NDim; ++i) {
//        os <<  "gVec[" + std::to_string(i) + "] : " + to_string(region.gVec.col(i).transpose()) + "\n";
//    }
//    os <<  "********************************************\n";
//}

template<int NDim> void Region<NDim>::writeToTxtFile(const std::string &FileName) const {
    std::ofstream myfile (FileName);
    for (int i=0; i<NDim; ++i) { myfile << origin(i,0) << " ";}
    myfile << "\n";
    for (int j=0; j<NDim; ++j) {
        for (int i=0; i<NDim; ++i) { myfile << gVec(i,j) << " ";}
            myfile << "\n";
    }
    myfile.close();
}

template<> void Region<3>::plot(const std::string& t_title) const {
    plotter3d  << *this ;
    if (t_title!= "") plotter3d << PLOTTITLE(t_title);
    plotter3d << NOLABEL;
    plotter3d << COLOR("red") << PLOT;
}
template<> void Region<2>::plot(const std::string& t_title) const {
    plotter3d  << *this ;
    if (t_title!= "") plotter3d << PLOTTITLE(t_title);
    plotter3d << NOLABEL;
    plotter3d << COLOR("red") << PLOT;
}
template<> void Region<1>::plot(const std::string& t_title) const {
    plotter2d  << *this << LABEL("Region");
    if (t_title!= "") plotter2d << PLOTTITLE(t_title);
    plotter2d << PLOT;
}

Plotter3D& operator << (Plotter3D& plotter, const Region<3>& region){
    bool polyModeActive = plotter.polygonModeActive;
    if (!polyModeActive){ plotter << POLYGONMODE; }
    plotter << OPENSTREAM;
    ArrayPoint<3,4> vertices;
    vertices << region.origin ,
                region.origin + region.gVec.col(0),
                region.origin + region.gVec.col(0) + region.gVec.col(1),
                region.origin + region.gVec.col(1);
    plotter << vertices;
    vertices << region.origin ,
                region.origin + region.gVec.col(0),
                region.origin + region.gVec.col(0) + region.gVec.col(2),
                region.origin + region.gVec.col(2);
    plotter << vertices;
    vertices << region.origin ,
                region.origin + region.gVec.col(1),
                region.origin + region.gVec.col(1) + region.gVec.col(2),
                region.origin + region.gVec.col(2);
    plotter << vertices;
    vertices << region.origin + region.gVec.col(2),
                region.origin + region.gVec.col(0)+ region.gVec.col(2),
                region.origin + region.gVec.col(0) + region.gVec.col(1)+ region.gVec.col(2),
                region.origin + region.gVec.col(1)+ region.gVec.col(2);
    plotter << vertices;
    vertices << region.origin + region.gVec.col(1),
                region.origin + region.gVec.col(0)+ region.gVec.col(1),
                region.origin + region.gVec.col(0) + region.gVec.col(2)+ region.gVec.col(1),
                region.origin + region.gVec.col(2)+ region.gVec.col(1);
    plotter << vertices;
    vertices << region.origin + region.gVec.col(0),
                region.origin + region.gVec.col(1)+ region.gVec.col(0),
                region.origin + region.gVec.col(1) + region.gVec.col(2)+ region.gVec.col(0),
                region.origin + region.gVec.col(2)+ region.gVec.col(0);
    plotter << vertices;
    if (!polyModeActive){ plotter << LINEMODE; }
    return plotter << CLOSESTREAM;
}
Plotter3D& operator << (Plotter3D& plotter, const Region<2>& region){
    ArrayPoint<2,4> vertices;
    vertices << region.origin ,
                region.origin + region.gVec.col(0),
                region.origin + region.gVec.col(0) + region.gVec.col(1),
                region.origin + region.gVec.col(1);
    if (plotter.polygonModeActive){
        return plotter << vertices ;
    } else {
        return plotter << POLYGONMODE << vertices << LINEMODE;
    }
}
Plotter2D& operator << (Plotter2D& plotter, const Region<2>& region){
    ArrayPoint<2,4> vertices;
    vertices << region.origin ,
                region.origin + region.gVec.col(0),
                region.origin + region.gVec.col(0) + region.gVec.col(1),
                region.origin + region.gVec.col(1);
    if (plotter.polygonModeActive){
        return plotter << vertices ;
    } else {
        return plotter << POLYGONMODE << vertices << LINEMODE;
    }
}
Plotter2D& operator << (Plotter2D& plotter, const Region<1>& region){
    ArrayPoint<1,2> vertices;
    vertices << region.origin ,
                region.origin + region.gVec.col(0);
    return plotter << ASPOINTS(vertices) << COLOR("black") << NOLABEL << vertices  ;
}

// Explicit template instantiation (For more info see: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl)

template class Region<1>;
template class Region<2>;
template class Region<3>;

template VectorPoint<1> calculate_gVec<1>(const Vect_Points<1> &t_gVec);
template VectorPoint<2> calculate_gVec<2>(const Vect_Points<2> &t_gVec);
template VectorPoint<3> calculate_gVec<3>(const Vect_Points<3> &t_gVec);

} // namespace Tortoise
