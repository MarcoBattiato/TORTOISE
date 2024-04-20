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
//  StringMap.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 4/11/20.
//
// >> EXAMPLE <<
//#include "StringMap.hpp"
//#include <deque>
//#include <iostream>
//
//class Bla : public StringMap<std::deque<double>> {
//
//};
//
//int main(int argc, const char * argv[]) {
//    Bla bla;
//    bla.emplace_back("uno",1);
//    bla.emplace_back("due",2);
//    bla.emplace_back("tre",3);
//    bla.emplace_back("trentatre",33);
//
//    std::cout << bla.matchList("tr*").size() << "\n=============\n";
//
//    for (auto i : bla.matchList("tr*")) { std::cout << *i << "\n";}
//    for (auto i : bla.matchListIndex("tr*")) { std::cout << bla[i] << "\n";}
//    for (auto i : bla.matchListId("tr*")) { std::cout << i << "\n";}
//
//    return 0;
//}


#ifndef StringMap_hpp
#define StringMap_hpp

#include <Generics/Utilities/StringUtilities.hpp>

#include <deque>
#include <cassert>

namespace Tortoise {

namespace Containers {

template <class DataType> class StringMap {
    //=======================================================
    // Properties
    //===================
    std::deque<std::string> ids;
    std::deque<DataType>    datas;
    
public:
    
    //=======================================================
    // Constructors
    //===================
    StringMap(){};                                                                              // Default contructor
    StringMap(const std::string& id_, const DataType& y): ids({id_}), datas({y}){};             // Single entry constructor
    StringMap(const StringMap& other): ids(other.ids), datas(other.datas){};                    // Copy constructor
    StringMap(StringMap&& other): ids(std::move(other.ids)), datas(std::move(other.datas)){};   // Move constructor
    // https://stackoverflow.com/questions/3413470/what-is-stdmove-and-when-should-it-be-used
    
    //=======================================================
    // Assignment
    //===================
    // Move assignment https://en.cppreference.com/w/cpp/language/move_assignment
    StringMap& operator=(const StringMap& other){ids = other.ids; datas = other.datas; return *this;}                   // Copy assignment
    StringMap& operator=(StringMap&& other){ids = std::move(other.ids); datas = std::move(other.datas); return *this;}  // Move assignment
    
    //=======================================================
    // Status
    //===================
    bool                        empty() const {return ids.empty();};
    unsigned long               size() const {return ids.size();};

    //=======================================================
    // Array subscripting
    //===================
    // Returns the argument of the function (i.e. the x variable)
    DataType&                   operator[](const unsigned long i) {return datas[i];};                  // Array Subscripting
    const DataType&             operator[](const unsigned long i) const {return datas[i];};            // Array Subscripting
    DataType&                   operator[](const std::string &id_)                                     // Array Subscripting
        {auto i = findID(id_); assert(i>=0 && i <size()); return datas[i];};
    const DataType&             operator[](const std::string &id_) const                               // Array Subscripting
        {auto i = findID(id_); assert(i>=0 && i <size()); return datas[i];};
    DataType&                   data(const unsigned long i) {return datas[i];};                        // Array Subscripting
    const DataType&             data(const unsigned long i) const {return datas[i];};                  // Array Subscripting
        
    // Returns the value of the function (i.e. the y variable)
    std::string&                id(const unsigned long i) {return ids[i];};                 // Array Subscripting ID
    const std::string&          id(const unsigned long i) const {return ids[i];};           // Array Subscripting ID
//    unsigned long               id(const std::string &id_) const {return findID(id_);}
    int                         id(const std::string &id_) const {return findID(id_);}
    //=======================================================
    // Filling of new elements
    //===================
    // It prevents the adding of data with an already existing id
    StringMap&     emplace_back(const std::string& id_, const DataType& y) {
        assert(findID(id_)==-1);
        ids.emplace_back(id_);
        datas.emplace_back(y);
        return *this;
    };
    template< class... Args > StringMap& emplace_back(const std::string& id_, Args&&... args) {
        assert(findID(id_)==-1);
        ids.emplace_back(id_);
        datas.emplace_back(args...);
        return *this;
    };
    // Notice that emplace_back is used when constructing objects. See https://stackoverflow.com/questions/4303513/push-back-vs-emplace-back
    
    //=======================================================
    // Erase
    //===================
    void            erase(const unsigned long i) {
        ids.erase(ids.begin() + i);
        datas.erase(datas.begin() + i);
    };
    void            pop_back() {
        ids.pop_back(); datas.pop_back();
    };
    
    //=======================================================
    // Iterator
    //===================
    std::vector<DataType*> matchList(const std::string& idWildCard) const {
        std::vector<DataType*> toreturn;
        for (int i=0; i<ids.size(); i++){ if (Utilities::matchString(idWildCard,ids[i])) {toreturn.emplace_back(&data(i));} }
        return toreturn;
    }
    std::vector<int> matchListIndex(const std::string& idWildCard) const {
        std::vector<int> toreturn;
        for (int i=0; i<ids.size(); i++){ if (Utilities::matchString(idWildCard,ids[i])) {toreturn.emplace_back(i);} }
        return toreturn;
    }
    std::vector<std::string> matchListId(const std::string& idWildCard) const {
        std::vector<std::string> toreturn;
        for (unsigned long i=0; i<ids.size(); i++){ if (Utilities::matchString(idWildCard,ids[i])) {toreturn.emplace_back(ids[i]);} }
        return toreturn;
    }
    
    int findID(const std::string &id_) const {
        for (int i=0; i<ids.size(); i++) { if (ids[i]==id_) return i;}
        return -1;
    };
};

} // namespace Containers

} // namespace Tortoise

#endif /* StringMap_hpp */
