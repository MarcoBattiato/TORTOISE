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
//  AnalyticFunctionMultiSpace.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 24/9/21.
//

#include "AnalyticFunctionMultiSpace.hpp"

namespace Tortoise {

AnalyticFunctionMultiSpacePoint k0(0);
AnalyticFunctionMultiSpacePoint k1(1);
AnalyticFunctionMultiSpacePoint k2(2);
AnalyticFunctionMultiSpacePoint k3(3);
AnalyticFunctionMultiSpacePoint k4(4);

AnalyticFunctionMultiSpacePointCoordinate k0x(0,0);
AnalyticFunctionMultiSpacePointCoordinate k0y(0,1);
AnalyticFunctionMultiSpacePointCoordinate k0z(0,2);
AnalyticFunctionMultiSpacePointCoordinate k1x(1,0);
AnalyticFunctionMultiSpacePointCoordinate k1y(1,1);
AnalyticFunctionMultiSpacePointCoordinate k1z(1,2);
AnalyticFunctionMultiSpacePointCoordinate k2x(2,0);
AnalyticFunctionMultiSpacePointCoordinate k2y(2,1);
AnalyticFunctionMultiSpacePointCoordinate k2z(2,2);
AnalyticFunctionMultiSpacePointCoordinate k3x(3,0);
AnalyticFunctionMultiSpacePointCoordinate k3y(3,1);
AnalyticFunctionMultiSpacePointCoordinate k3z(3,2);
AnalyticFunctionMultiSpacePointCoordinate k4x(4,0);
AnalyticFunctionMultiSpacePointCoordinate k4y(4,1);
AnalyticFunctionMultiSpacePointCoordinate k4z(4,2);

} // namespace Tortoise
