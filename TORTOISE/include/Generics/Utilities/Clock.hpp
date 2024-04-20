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
//  Clock.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 1/6/19.
//

#ifndef Clock_hpp
#define Clock_hpp

#include <stdio.h>
#include <chrono>
    
namespace Tortoise {

namespace Utilities {

class Clock {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start, finish;
    std::chrono::duration<double> elapsed;
 
public:
    void tic() {start = std::chrono::high_resolution_clock::now();}
    double toc() {
        finish = std::chrono::high_resolution_clock::now();
        elapsed = finish - start;
        return elapsed.count();
    }
        
};
    
extern Clock stopWatch;

} // namespace Utilities 

} // namespace Tortoise

#endif /* Clock_hpp */
