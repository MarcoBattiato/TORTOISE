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
//  UnstructuredMesh.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 7/8/21.
//
//  INCOMPLETE

#ifndef UnstructuredMesh_hpp
#define UnstructuredMesh_hpp

#include <Geometry/GeometryCore/Geometry.hpp>
#include <functional>

namespace Tortoise {

class UnstructuredMesh1D {
//*******************************
// Dimensionality constants
//*******************************
    static constexpr int                            nVertElement = 2;                  // Number of vertices in Element
    static constexpr int                            nElementSection = 1;               // Number of element per reference Section
//    static const VectorElement<NDim>                vertElemInRefSec;                       // Vertices of Ref Elements (They all share the first vertex at the origin)

    // Geometric Attributes
public:
    const Real                                     origin;                                 // origin of the Mesh
    const Real                                     gVec;                                   // Mesh sides vectors
    const int                                      numberElements;                         // Number of sectioning
    const int                                      dataVectorDim;                          // Total number of basis functions
    
    const Eigen::Matrix<Real, 1, Eigen::Dynamic>   meshPoints;                             // Position of mesh points includint first and last
    
    // Methods
public:
    //=======================================================
    // Constructors
    // !!! Notice !!! gVec is the length of the and NOT the end points of the mesh.
    // The endpoint of the mesh is origin+gVec
    // The origin has to be passed in ABSOLUTE position
    //===================
    UnstructuredMesh1D(const Real origin, const Real gVec, const int numberElements);     // Creates a uniform mesh
    UnstructuredMesh1D(const Eigen::Matrix<Real, 1, Eigen::Dynamic> meshPoints);               // Creates a mesh with given meshpoints
    UnstructuredMesh1D(const Real origin, const Real gVec, const int numberElements,
                       const std::function<Real(const Real)>& density);                     // Creates a mesh with a mesh point density proportional to the given function
    //===================

    
};
    

} // namespace Tortoise

#endif /* UnstructuredMesh_hpp */
