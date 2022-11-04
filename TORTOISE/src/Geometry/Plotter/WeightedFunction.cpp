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
//  WeightedFunction.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 19/10/20.
//

#include <Geometry/Plotter/WeightedFunction.hpp>
#include <Generics/Utilities/StringUtilities.hpp>

namespace Tortoise {


template< typename T >
std::string int_to_hex( T i )
{
  std::stringstream stream;
  stream << std::setfill ('0') << std::setw(2)
         << std::hex << i;
  return stream.str();
}

Plotter3D& operator << (Plotter3D& plotter, const WeightedFunction<2>& function){
        if(plotter.numberplots++>0) plotter.plotCommand +=" , ";
        plotter.plotCommand += " '-' with lines lw 0.3";
//    plotter.plotFormat += "set hidden3d\n";
    
        for (auto elem = function.surface.mesh->elementIterator(); elem.unfinished; ++elem){
            for (auto vertex = function.surface.mesh->vertexInElemIterator(); vertex.unfinished; ++vertex){
                Point<2> pos((*function.surface.mesh)(elem,vertex));
                plotter.plotData += std::to_string(pos(0)) + " " + std::to_string(pos(1)) + " " + std::to_string(apply(function.surface(elem).vec,refNodes<2>.col(vertex.currentVertID))(0,0)) + "\n";
            }
            Point<2> pos ((*function.surface.mesh)(elem,function.surface.mesh->vertexInElemIterator()));
            plotter.plotData += std::to_string(pos(0)) + " " + std::to_string(pos(1)) + " " + std::to_string(apply(function.surface(elem).vec,refNodes<2>.col(0))(0,0)) + "\n";
            plotter.plotData += "\n\n";
            plotter.numberobjects++;
            plotter.plotFormat += "set obj " + std::to_string(plotter.numberobjects) + " polygon from "  ;
            for (auto vertex = function.surface.mesh->vertexInElemIterator(); vertex.unfinished; ++vertex){
                Point<2> pos((*function.surface.mesh)(elem,vertex));
                plotter.plotFormat += std::to_string(pos(0)) + "," + std::to_string(pos(1)) + "," + std::to_string(apply(function.surface(elem).vec,refNodes<2>.col(vertex.currentVertID))(0,0)) + " to ";
            }
             pos = ((*function.surface.mesh)(elem,function.surface.mesh->vertexInElemIterator()));
            plotter.plotFormat += std::to_string(pos(0)) + "," + std::to_string(pos(1)) + "," + std::to_string(apply(function.surface(elem).vec,refNodes<2>.col(0))(0,0)) ;
            plotter.plotFormat += " fc rgb '#" + int_to_hex(static_cast<int>((1.0 - apply(function.filling(elem).vec,refNodes<2>.col(0))(0,0))*255)) + int_to_hex(static_cast<int>(apply(function.filling(elem).vec,refNodes<2>.col(0))(0,0)*255)) + "0000' fillstyle solid noborder \n";
        }
        plotter.plotData += "e\n";
//    std::cout <<plotter.plotFormat ;
    return plotter;
}

} // namespace Tortoise
