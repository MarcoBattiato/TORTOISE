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
//  MathFieldSpace.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 17/7/20.
//

#ifndef MathFieldSpace_hpp
#define MathFieldSpace_hpp

namespace Tortoise {

template <class fieldType> class MathFieldSpace {
public:
//    virtual fieldType& operator+=(const fieldType& other) = 0;
//    virtual fieldType& operator-=(const fieldType& other) = 0;
//    virtual fieldType& operator*=(const fieldType& other) = 0;
//    virtual fieldType& operator/=(const fieldType& other) = 0;
    
    friend fieldType operator+(fieldType lhs, const fieldType& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend fieldType operator-(fieldType lhs, const fieldType& rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend fieldType operator*(fieldType lhs, const fieldType& rhs){
        lhs *= rhs;
        return lhs;
    }
    friend fieldType operator/(fieldType lhs, const fieldType& rhs){
        lhs /= rhs;
        return lhs;
    }
};



template <class fieldType, class friendfieldType> class AsymmetricMathFieldSpace {
public:
//    virtual fieldType& operator+=(const friendfieldType& other) = 0;
//    virtual fieldType& operator-=(const friendfieldType& other) = 0;
//    virtual fieldType& operator*=(const friendfieldType& other) = 0;
//    virtual fieldType& operator/=(const friendfieldType& other) = 0;
//  The user has to define the following function as well
//    friend fieldType operator/(const friendfieldType& lhs, fieldType rhs) ;
//    virtual fieldType operator-() const = 0;
    
    
    friend fieldType operator+(fieldType lhs, const friendfieldType& rhs) {
        lhs += rhs;
        return lhs;
    }
    friend fieldType operator+(const friendfieldType& lhs, fieldType rhs) {
        rhs += lhs;
        return rhs;
    }
    friend fieldType operator-(fieldType lhs, const friendfieldType& rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend fieldType operator-(const friendfieldType& lhs, fieldType rhs) {
        rhs -= lhs;
        return -rhs;
    }
    friend fieldType operator*(fieldType lhs, const friendfieldType& rhs) {
        lhs *= rhs;
        return lhs;
    }
    friend fieldType operator*(const friendfieldType& lhs, fieldType rhs) {
        rhs *= lhs;
        return rhs;
    }
    friend fieldType operator/(fieldType lhs, const friendfieldType& rhs) {
        lhs /= rhs;
        return lhs;
    }
};

template <class fieldType, class friendfieldType> class AsymmetricMathFieldSpaceByValue {
    // Used when it is beneficial to pass the friendfieldType by value (for instance when friendfieldType is a scalar type, like double)
public:
//    virtual fieldType& operator+=(const friendfieldType other) = 0;
//    virtual fieldType& operator-=(const friendfieldType other) = 0;
//    virtual fieldType& operator*=(const friendfieldType other) = 0;
//    virtual fieldType& operator/=(const friendfieldType other) = 0;
//  The user has to define the following function as well
//    friend fieldType operator/(const friendfieldType& lhs, fieldType rhs) ;
//    virtual fieldType operator-() const = 0;
    
    
    friend fieldType operator+(fieldType lhs, const friendfieldType rhs) {
        lhs += rhs;
        return lhs;
    }
    friend fieldType operator+(const friendfieldType lhs, fieldType rhs) {
        rhs += lhs;
        return rhs;
    }
    friend fieldType operator-(fieldType lhs, const friendfieldType rhs) {
        lhs -= rhs;
        return lhs;
    }
    friend fieldType operator-(const friendfieldType lhs, fieldType rhs) {
        rhs -= lhs;
        return -rhs;
    }
    friend fieldType operator*(fieldType lhs, const friendfieldType rhs) {
        lhs *= rhs;
        return lhs;
    }
    friend fieldType operator*(const friendfieldType lhs, fieldType rhs) {
        rhs *= lhs;
        return rhs;
    }
    friend fieldType operator/(fieldType lhs, const friendfieldType rhs) {
        lhs /= rhs;
        return lhs;
    }
};

} // namespace Tortoise 

#endif /* MathFieldSpace_hpp */
