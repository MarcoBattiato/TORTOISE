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
//  VectorSpace.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 16/7/20.
//

#ifndef VectorSpace_hpp
#define VectorSpace_hpp

namespace Tortoise {

template <class groupType> class GroupSpace {
public:
//    virtual groupType& operator+=(const groupType& other) = 0;
//    virtual groupType& operator-=(const groupType& other) = 0;
    
    friend groupType operator+(groupType lhs, const groupType& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend groupType operator-(groupType lhs, const groupType& rhs) {
        lhs -= rhs;
        return lhs;
    }
};



template <class vecType, class scalarType> class VectorSpace {
public:
//    virtual vecType& operator+=(const vecType& other) = 0;
//    virtual vecType& operator-=(const vecType& other) = 0;
//    virtual vecType& operator*=(const scalarType scalar) = 0;
    
    
    friend vecType operator+(vecType lhs, const vecType& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vecType operator-(vecType lhs, const vecType& rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend vecType operator*(vecType lhs, const scalarType scalar){
        lhs *= scalar;
        return lhs;
    }
    friend vecType operator*(const scalarType scalar, vecType rhs){
        rhs *= scalar;
        return rhs;
    }
};

// Adds asymmetric operations with friendVecType which must be a smaller vector space
// The relationship is asymmetrisc since operations the between vecType and friendVecType will return vecType
template <class vecType, class friendVecType> class AsymmetricVectorSpace {
public:
//    virtual vecType& operator+=(const friendVecType& other) = 0;
//    virtual vecType& operator-=(const friendVecType& other) = 0;
//    virtual vecType operator-() const = 0;
    
    friend vecType operator+(vecType lhs, const friendVecType& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend vecType operator-(vecType lhs, const friendVecType& rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend vecType operator+(const friendVecType& lhs, vecType rhs) {
        rhs += lhs;
        return rhs;
    }
    friend vecType operator-(const friendVecType& lhs, vecType rhs) {
        rhs -= lhs;
        return -rhs;
    }
};

} // namespace Tortoise 

#endif /* VectorSpace_hpp */
