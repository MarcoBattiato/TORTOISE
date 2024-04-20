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
//  Region.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//
//  A Region is a NDim-dimensional parallelogram. Different parts of a region can be meshed.
//  Certain integral operations between functions make sense only if the functions are defined over the same Region, even if defined on different meshes.
// The most common use is to construct a Brillouin zone.

#ifndef Region_hpp
#define Region_hpp

#include <Geometry/GeometryCore/Geometry.hpp>
#include <Geometry/Plotter/Plotter.hpp>

#include <vector>
#include <type_traits>

namespace Tortoise {

//*******************************
// Types definitions
//*******************************
template <int NDim> using Vect_Points               = std::vector<Eigen::Matrix<Real, NDim, 1>, Eigen::aligned_allocator<Eigen::Matrix<Real, NDim, 1> > > ; // Useful for input only

//*******************************
// Main class definition
//*******************************

template<int NDim> class Region {         // NDim = Number of Dimensions of the space
public:
    // Geometric Attributes
    Point<NDim> const                               origin;                                 // origin of the Region
    ArrayPoint<NDim,NDim> const                     gVec;                                   // Region sides vectors
    
public:
    //=======================================================
    // Constructors
    // !!! Notice !!! gVec are the lattice vectors and NOT the end points of the mesh.
    // The endpoints of the mesh are origin+gVec[1], origin+gVec[2], ...
    //===================
    Region(const Point<NDim>& t_origin, const Vect_Points<NDim>& t_gVec);
    Region(const Real t_origin, const Real t_gVec) requires (NDim==1);
    
    Region(const std::string &FileName);    // Constructs a region reading the data from file. The file must be built with the method writeToTxtFile declared below
//    Region(std::ostream &t_os);             // Constructs a region from an open stream
    //===================
    
    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& nodevalue){
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> toreturn = nodevalue;
        toreturn.rightCols(toreturn.cols()-1).colwise() -= toreturn.col(0);
        return toreturn;
    };
    
    //=======================================================
    // Geometry
    //===================
    Point<NDim> toRelative(const Point<NDim>& absCoord) const;
    Point<NDim> toAbsolute(const Point<NDim>& relCoord) const;
    
    //=======================================================
    // Comparison
    //===================
    bool equivalentTo(const Region<NDim>& mesh2) const;
    bool equivalentTo(const std::string &FileName) const;
    
    //=======================================================
    // I/O
    //===================
    void plot(const std::string& t_title = "") const;
    void writeToTxtFile(const std::string &FileName) const;
//    friend std::ostream &operator<<(std::ostream &t_os, Region<NDim> const& t_region);
    
public:
    //=======================================================
    // Implementation details
    //===================
    // Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
    // See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

// I/O operations
Plotter3D& operator << (Plotter3D& plotter, const Region<3>& region);
Plotter3D& operator << (Plotter3D& plotter, const Region<2>& region);
Plotter2D& operator << (Plotter2D& plotter, const Region<2>& region);
Plotter2D& operator << (Plotter2D& plotter, const Region<1>& region);

} // namespace Tortoise
#endif /* Region_hpp */
