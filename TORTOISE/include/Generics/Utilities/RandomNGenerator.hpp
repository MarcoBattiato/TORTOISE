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
//  RandomNGenerator.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 15/10/20.
//

#ifndef RandomNGenerator_hpp
#define RandomNGenerator_hpp

#include <stdio.h>
#include <random>
#include <vector>
#include <iterator>
#include <algorithm>

namespace Tortoise {

class RandomNGenerator {
    std::default_random_engine                      generator;
public:
    
    void seed (const int t_seed){generator.seed(t_seed);}
    template <class contType> contType randCont(const contType min = 0.0, const contType max = 1.0) {
        std::uniform_real_distribution<contType> distribution(min,max);
        return distribution(generator);
    }
    template <class discrType> discrType randDiscr(const discrType min = 0, const discrType max = 10) {
        std::uniform_int_distribution<discrType> distribution(min,max);
        return distribution(generator);
    }
    template <class contType> std::vector<contType> randSeqCont(const int numberofRN, const contType min = 0.0, const contType max = 1.0) {
        std::vector<contType> toreturn;
        std::uniform_real_distribution<contType> distribution(min,max);
        std::generate_n(std::back_inserter(toreturn), numberofRN, [this,&distribution] () mutable { return distribution(this->generator); });
        return toreturn;
    }
    template <class discrType> std::vector<discrType> randSeqDiscr(const int numberofRN, const discrType min = 0, const discrType max = 10) {
        std::vector<discrType> toreturn;
        std::uniform_int_distribution<discrType> distribution(min,max);
        std::generate_n(std::back_inserter(toreturn), numberofRN, [this,&distribution] () mutable{ return distribution(this->generator); });
        return toreturn;
    }
    
};


extern RandomNGenerator randGen;

} // namespace Tortoise


#endif /* RandomNGenerator_hpp */
