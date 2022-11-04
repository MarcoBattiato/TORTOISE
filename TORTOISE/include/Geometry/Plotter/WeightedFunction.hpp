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
//  WeightedFunction.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 19/10/20.
//
// This class allows to join two Function, in which the first one is used to create a surface, and the second one is used as filling intensity
// This is only intended for plotting, and it is currently available only for 2D

#ifndef WeightedFunction_hpp
#define WeightedFunction_hpp

#include <Geometry/Structured/FunctionSpace/Function.hpp>

#include <stdio.h>

namespace Tortoise {

template<int NDim> class WeightedFunction {
public:
    const Function<NDim>   surface;
    const Function<NDim>   filling;
    
    WeightedFunction(const Function<NDim>& t_surface, const Function<NDim>& t_filling): surface(t_surface), filling(t_filling) {assert(t_surface.mesh == t_filling.mesh);};
    
};

Plotter3D& operator << (Plotter3D& plotter, const WeightedFunction<2>& function);

} // namespace Tortoise


#endif /* WeightedFunction_hpp */
