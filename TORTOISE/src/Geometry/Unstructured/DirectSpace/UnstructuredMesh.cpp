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
//  UnstructuredMesh.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 7/8/21.
//
//  INCOMPLETE


#include <Geometry/Unstructured/DirectSpace/UnstructuredMesh.hpp>

#include <cassert>

namespace Tortoise {

//=======================================================
// Constructors
//===================
// Constructor helper functions
Eigen::Matrix<Real, 1, Eigen::Dynamic> construct_meshPoints_from_density(const Real origin, const Real gVec, const int numberElements,const std::function<Real(Real)>& density){
    Eigen::Array<Real, 1, Eigen::Dynamic> toreturn;
    toreturn.resize(Eigen::NoChange, numberElements+1);
    toreturn(0)=origin;
    Real rescaleDens = 1.;
    int iteration = 0;
    do {
        ++iteration;
        for (int i=0; i<numberElements; ++i){
            assert(density(toreturn(i))>0. );
            toreturn(i+1) = toreturn(i) + rescaleDens / density(toreturn(i));
        }
        rescaleDens = 0.5 * rescaleDens + 0.5 * rescaleDens * gVec/(toreturn(numberElements)-origin);
    } while ( (std::abs(toreturn(numberElements)-origin-gVec) > 0.0001*gVec/numberElements) and (iteration < 100));
    return ((toreturn - origin) * gVec / (toreturn(numberElements)-origin) + origin).matrix();
}
// Constructors
UnstructuredMesh1D::UnstructuredMesh1D(const Real t_origin, const Real t_gVec, const int t_numberElements): origin(t_origin), gVec(t_gVec), numberElements(t_numberElements), dataVectorDim(2*t_numberElements), meshPoints(Eigen::Matrix<Real, 1, Eigen::Dynamic>::LinSpaced(t_numberElements+1, t_origin, t_origin+t_gVec)){} ;
UnstructuredMesh1D::UnstructuredMesh1D(const Eigen::Matrix<Real, 1, Eigen::Dynamic> t_meshPoints): origin(t_meshPoints(0)), gVec(t_meshPoints(t_meshPoints.size()-1)-t_meshPoints(0)), numberElements(t_meshPoints.size()-1), dataVectorDim(2*(t_meshPoints.size()-1)), meshPoints(t_meshPoints){}
UnstructuredMesh1D::UnstructuredMesh1D(const Real t_origin, const Real t_gVec, const int t_numberElements,const std::function<Real(Real)>& t_f): origin(t_origin), gVec(t_gVec), numberElements(t_numberElements), dataVectorDim(2*t_numberElements), meshPoints(construct_meshPoints_from_density(t_origin,t_gVec,t_numberElements,t_f)) {};

} // namespace Tortoise
