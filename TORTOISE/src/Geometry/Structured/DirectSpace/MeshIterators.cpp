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
//  MeshIterators.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 23/6/21.
//

#include <Geometry/Structured/DirectSpace/MeshIterators.hpp>

namespace Tortoise {

//*******************************
// Section Iterator
//*******************************

//=======================================================
// Constructors
//===================
// Constructor from mesh
//template<int NDim> SectionIterator<NDim>::SectionIterator(const Mesh<NDim>& mesh): maxcoord(mesh.nSecSplits), lcoord(CartIndex<NDim>::Zero()), hcoord(mesh.nSecSplits), currentcoord(CartIndex<NDim>::Zero()), unfinished((mesh.nSecSplits.array() != 0).all()){};
//template<int NDim> SectionIterator<NDim>::SectionIterator(const Mesh<NDim>* mesh): SectionIterator(*mesh){};

// Explicit constructor
template<int NDim> SectionIterator<NDim>::SectionIterator(const CartIndex<NDim> t_maxcoord, const CartIndex<NDim> t_lcoord, const CartIndex<NDim> t_hcoord): maxcoord(t_maxcoord), lcoord(t_lcoord), hcoord(t_hcoord), currentcoord(t_lcoord), unfinished((t_lcoord.array() != t_hcoord.array()).all())
{
    assert( (lcoord.array() >= 0).all() && (lcoord.array() <= maxcoord.array()).all());
    assert( (hcoord.array() >= 0).all() && (lcoord.array() <= maxcoord.array()).all());
};

// Copy and move constructors
template<int NDim> SectionIterator<NDim>::SectionIterator(const SectionIterator<NDim>& other): lcoord(other.lcoord), hcoord(other.hcoord), maxcoord(other.maxcoord), currentcoord(other.currentcoord), unfinished(other.unfinished){};
template<int NDim> SectionIterator<NDim>::SectionIterator(SectionIterator<NDim>&& other): lcoord(std::move(other.lcoord)), hcoord(std::move(other.hcoord)), maxcoord(std::move(other.maxcoord)), currentcoord(std::move(other.currentcoord)), unfinished(std::move(other.unfinished)){};

//=======================================================
// Increment and assignment
//===================
// Prefix increment
template<> SectionIterator<1>& SectionIterator<1>::operator++ () {
    assert(unfinished);
    ++currentcoord(0);
    if (currentcoord(0)==hcoord(0)) unfinished=false;
    currentcoord(0) %= maxcoord(0);
    return *this;
};
template<> SectionIterator<2>& SectionIterator<2>::operator++ () {
    assert(unfinished);
    ++currentcoord(0);
    if (currentcoord(0)==hcoord(0)) {
        ++currentcoord(1);
        if (currentcoord(1)==hcoord(1)) unfinished=false;
        currentcoord(0) = lcoord(0);
    }
    currentcoord(0) %= maxcoord(0);
    currentcoord(1) %= maxcoord(1);
    return *this;
};
template<> SectionIterator<3>& SectionIterator<3>::operator++ () {
    assert(unfinished);
    ++currentcoord(0);
    if (currentcoord(0)==hcoord(0)) {
        ++currentcoord(1);
        if (currentcoord(1)==hcoord(1)) {
            ++currentcoord(2);
            if (currentcoord(2)==hcoord(2)) unfinished=false;
            currentcoord(1) = lcoord(1);
        }
        currentcoord(0) = lcoord(0);
    }
    currentcoord(0) %= maxcoord(0);
    currentcoord(1) %= maxcoord(1);
    currentcoord(2) %= maxcoord(2);
    return *this;
};
// Assignment
template<int NDim> SectionIterator<NDim>& SectionIterator<NDim>::operator=(const SectionIterator<NDim>& other) {
    lcoord = other.lcoord;
    hcoord = other.hcoord;
    maxcoord = other.maxcoord;
    currentcoord = other.currentcoord;
    unfinished = other.unfinished;
    return *this;
}
template<int NDim> SectionIterator<NDim>& SectionIterator<NDim>::operator=(SectionIterator<NDim>&& other) {
    lcoord = std::move(other.lcoord);
    hcoord = std::move(other.hcoord);
    maxcoord = std::move(other.maxcoord);
    currentcoord = std::move(other.currentcoord);
    unfinished = std::move(other.unfinished);
    return *this;
}

//=======================================================
// I/O
//===================
template<int NDim> std::ostream &operator<<(std::ostream &os, SectionIterator<NDim> const& iter){
    std::string output;
    os << "structure : [" << iter.lcoord.transpose() << "][" << iter.hcoord.transpose() << "][" << iter.maxcoord.transpose() << "]  ";
    os << "status : [" << iter.currentcoord.transpose() << "][" << iter.unfinished <<  "]";
    return os;
}





//*******************************
// Element Iterators
//*******************************

//*******************************
// Element Iterators : meshes
//*******************************

//=======================================================
// Constructors
//===================
// Constructor from mesh
//template<int NDim> MeshElementIterator<NDim>::MeshElementIterator(const Mesh<NDim>& mesh): maxElemID(mesh.numberElements), currentElementID(0), unfinished(mesh.numberElements != 0){};
//template<int NDim> MeshElementIterator<NDim>::MeshElementIterator(const Mesh<NDim>* mesh): MeshElementIterator(*mesh){};

// Explicit constructor
template<int NDim> MeshElementIterator<NDim>::MeshElementIterator(const int t_maxElemID, const int t_currentElementID): maxElemID(t_maxElemID), currentElementID(t_currentElementID), unfinished(t_currentElementID < t_maxElemID){};

template<int NDim> MeshElementIterator<NDim>::MeshElementIterator(){};

// Copy and move constructors
template<int NDim> MeshElementIterator<NDim>::MeshElementIterator(const MeshElementIterator<NDim>& other): maxElemID(other.maxElemID), currentElementID(other.currentElementID), unfinished(other.unfinished){};
template<int NDim> MeshElementIterator<NDim>::MeshElementIterator(MeshElementIterator<NDim>&& other): maxElemID(other.maxElemID), currentElementID(other.currentElementID), unfinished(other.unfinished){};

//=======================================================
// Increment and assignment
//===================
// Prefix increment
template<int NDim> MeshElementIterator<NDim>& MeshElementIterator<NDim>::operator++() {
    assert(unfinished);
    ++currentElementID;
    if (currentElementID==maxElemID) unfinished=false;
    return *this;
};

// Assignment
template<int NDim> MeshElementIterator<NDim>&  MeshElementIterator<NDim>::operator=(const MeshElementIterator<NDim>& other) {
    maxElemID = other.maxElemID; currentElementID = other.currentElementID; unfinished = other.unfinished;
    return *this;
}
template<int NDim> MeshElementIterator<NDim>&  MeshElementIterator<NDim>::operator=(MeshElementIterator<NDim>&& other){
    maxElemID = other.maxElemID; currentElementID = other.currentElementID; unfinished = other.unfinished;
    return *this;
}
template<int NDim> std::ostream &operator<<(std::ostream &os, MeshElementIterator<NDim> const& iter){
    std::string output;
    os << "structure : [" << iter.maxElemID << "]  status : [" << iter.currentElementID << "][" << iter.unfinished <<  "]";
    return os;
}




//*******************************
// Element Iterators : mesh subsets
//*******************************

//=======================================================
// Constructors
//===================
// Constructor from mesh subsets

int calculateCurrentElementID(const CartIndex<1>& maxcoord, const CartIndex<1>& currentSeccoord, const int currentElementinSectionID){
    return currentSeccoord(0);
}
int calculateCurrentElementID(const CartIndex<2>& maxcoord, const CartIndex<2>& currentSeccoord, const int currentElementinSectionID){
    return 2*(currentSeccoord(0) + currentSeccoord(1)*maxcoord(0)) + currentElementinSectionID;
}
int calculateCurrentElementID(const CartIndex<3>& maxcoord, const CartIndex<3>& currentSeccoord, const int currentElementinSectionID){
    return 6*(currentSeccoord(0) + (currentSeccoord(1) + currentSeccoord(2)*maxcoord(1))*maxcoord(0)) + currentElementinSectionID;
}

// Explicit constructor
template<int NDim>
MeshSubsetElementIterator<NDim>::MeshSubsetElementIterator(const CartIndex<NDim> t_maxcoord, const CartIndex<NDim> t_lcoord, const CartIndex<NDim> t_hcoord)
: maxcoord(t_maxcoord), lcoord(t_lcoord), hcoord(t_hcoord), currentSeccoord(t_lcoord), currentElementinSectionID(0),
currentElementID(calculateCurrentElementID(t_maxcoord, t_lcoord, 0)), unfinished((t_lcoord.array() != t_hcoord.array()).all())
{
//    assert( (lcoord.array() >= 0).all() && (lcoord.array() < maxcoord.array()).all());
    assert( (lcoord.array() >= 0).all() && (lcoord.array() <= maxcoord.array()).all());
    assert( (hcoord.array() >= 0).all() && (hcoord.array() <= maxcoord.array()).all());
};

// Copy and move constructors
template<int NDim>
MeshSubsetElementIterator<NDim>::MeshSubsetElementIterator(const MeshSubsetElementIterator<NDim>& other)
: lcoord(other.lcoord), hcoord(other.hcoord), maxcoord(other.maxcoord), currentSeccoord(other.currentSeccoord), currentElementinSectionID(other.currentElementinSectionID), currentElementID(other.currentElementID), unfinished(other.unfinished){};

template<int NDim>
MeshSubsetElementIterator<NDim>::MeshSubsetElementIterator(MeshSubsetElementIterator<NDim>&& other): lcoord(std::move(other.lcoord)), hcoord(std::move(other.hcoord)), maxcoord(std::move(other.maxcoord)), currentSeccoord(std::move(other.currentSeccoord)), currentElementinSectionID(other.currentElementinSectionID), currentElementID(other.currentElementID), unfinished(std::move(other.unfinished)){};

//=======================================================
// Increment and assignment
//===================
// Prefix increment
template<> MeshSubsetElementIterator<1>& MeshSubsetElementIterator<1>::operator++ () {
    assert(unfinished);
    ++currentSeccoord(0);
    if (currentSeccoord(0)==hcoord(0)) {unfinished=false; }
    currentSeccoord(0) %= maxcoord(0);
    currentElementID = calculateCurrentElementID(maxcoord, currentSeccoord, currentElementinSectionID);
    return *this;
};
template<> MeshSubsetElementIterator<2>& MeshSubsetElementIterator<2>::operator++ () {
    assert(unfinished);
    ++currentElementinSectionID;
    if(currentElementinSectionID==2){
        ++currentSeccoord(0);
        if (currentSeccoord(0)==hcoord(0)) {
            ++currentSeccoord(1);
            if (currentSeccoord(1)==hcoord(1)) unfinished=false;
            currentSeccoord(0) = lcoord(0);
        }
        currentSeccoord(0) %= maxcoord(0);
        currentSeccoord(1) %= maxcoord(1);
        currentElementinSectionID = 0;
    }
    currentElementID = calculateCurrentElementID(maxcoord, currentSeccoord, currentElementinSectionID);
    return *this;
};
template<> MeshSubsetElementIterator<3>& MeshSubsetElementIterator<3>::operator++ () {
    assert(unfinished);
    ++currentElementinSectionID;
    if(currentElementinSectionID==6){
        ++currentSeccoord(0);
        if (currentSeccoord(0)==hcoord(0)) {
            ++currentSeccoord(1);
            if (currentSeccoord(1)==hcoord(1)) {
                ++currentSeccoord(2);
                if (currentSeccoord(2)==hcoord(2)) unfinished=false;
                currentSeccoord(1) = lcoord(1);
            }
            currentSeccoord(0) = lcoord(0);
        }
        currentSeccoord(0) %= maxcoord(0);
        currentSeccoord(1) %= maxcoord(1);
        currentSeccoord(2) %= maxcoord(2);
        currentElementinSectionID = 0;
    }
    currentElementID = calculateCurrentElementID(maxcoord, currentSeccoord, currentElementinSectionID);
    return *this;
};

// Assignment
template<int NDim> MeshSubsetElementIterator<NDim>& MeshSubsetElementIterator<NDim>::operator=(const MeshSubsetElementIterator<NDim>& other) {
    lcoord = other.lcoord;
    hcoord = other.hcoord;
    maxcoord = other.maxcoord;
    currentSeccoord = other.currentSeccoord;
    currentElementinSectionID = other.currentElementinSectionID;
    currentElementID = other.currentElementID;
    unfinished = other.unfinished;
    return *this;
}
template<int NDim> MeshSubsetElementIterator<NDim>& MeshSubsetElementIterator<NDim>::operator=(MeshSubsetElementIterator<NDim>&& other) {
    lcoord = std::move(other.lcoord);
    hcoord = std::move(other.hcoord);
    maxcoord = std::move(other.maxcoord);
    currentSeccoord = std::move(other.currentSeccoord);
    currentElementinSectionID = other.currentElementinSectionID;
    currentElementID = other.currentElementID;
    unfinished = std::move(other.unfinished);
    return *this;
}

//template<> int MeshSubsetElementIterator<1>::currentElementID() const {
//    return currentSeccoord(0);
//}
//template<> int MeshSubsetElementIterator<2>::currentElementID() const{
//    return 2*(currentSeccoord(0) + currentSeccoord(1)*maxcoord(0)) + currentElementinSectionID;
//}
//template<> int MeshSubsetElementIterator<3>::currentElementID() const{
//    return 6*(currentSeccoord(0) + (currentSeccoord(1) + currentSeccoord(2)*maxcoord(1))*maxcoord(0)) + currentElementinSectionID;
//}



//=======================================================
// I/O
//===================
template<int NDim> std::ostream &operator<<(std::ostream &os, MeshSubsetElementIterator<NDim> const& iter){
    std::string output;
    os << "structure : [" << iter.lcoord.transpose() << "][" << iter.hcoord.transpose() << "][" << iter.maxcoord.transpose() << "]  ";
    os << "status : [" << iter.currentSeccoord.transpose() << " (" << iter.currentElementinSectionID << ")][" << iter.unfinished <<  "]";
    return os;
}




//*******************************
// Vertex Iterator
//*******************************

template<int Ndim> VertexInElemIterator<Ndim>::VertexInElemIterator (const VertexInElemIterator<Ndim>& other) : currentVertID(other.currentVertID), unfinished(other.unfinished) {};
template<int Ndim> VertexInElemIterator<Ndim>::VertexInElemIterator (VertexInElemIterator<Ndim>&& other) : currentVertID(std::move(other.currentVertID)), unfinished(std::move(other.unfinished)) {};
template<int Ndim> VertexInElemIterator<Ndim>::VertexInElemIterator () : currentVertID(0), unfinished(true) {};
template<int Ndim> VertexInElemIterator<Ndim>& VertexInElemIterator<Ndim>::operator++() {
    assert(unfinished);
    ++currentVertID;
    if (currentVertID == maxcoord) unfinished=false;
    return *this;
}
template<int Ndim> VertexInElemIterator<Ndim>& VertexInElemIterator<Ndim>::operator=(const VertexInElemIterator<Ndim>& other) {
    currentVertID = other.currentVertID; unfinished = other.unfinished;
    return *this;
}
template<int Ndim> VertexInElemIterator<Ndim>& VertexInElemIterator<Ndim>::operator=(VertexInElemIterator<Ndim>&& other) {
    currentVertID = std::move(other.currentVertID); unfinished = std::move(other.unfinished);
    return *this;
}
//=======================================================
// I/O
//===================
template<int NDim> std::ostream &operator<<(std::ostream &os, VertexInElemIterator<NDim> const& iter){
    std::string output;
    os << "structure : [" << iter.maxcoord << "status : [" << iter.currentVertID << "][" << iter.unfinished <<  "]";
    return os;
}




// ===========================================
// ===========================================
// ===========================================
// Explicit template instantiation (For more info see: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl)
template class SectionIterator<1>;
template class SectionIterator<2>;
template class SectionIterator<3>;

template  std::ostream &operator<<(std::ostream &os, const SectionIterator<1> & iter);
template  std::ostream &operator<<(std::ostream &os, const SectionIterator<2> & iter);
template  std::ostream &operator<<(std::ostream &os, const SectionIterator<3> & iter);

template class MeshElementIterator<1>;
template class MeshElementIterator<2>;
template class MeshElementIterator<3>;

template  std::ostream &operator<<(std::ostream &os, const MeshElementIterator<1> & iter);
template  std::ostream &operator<<(std::ostream &os, const MeshElementIterator<2> & iter);
template  std::ostream &operator<<(std::ostream &os, const MeshElementIterator<3> & iter);

template class MeshSubsetElementIterator<1>;
template class MeshSubsetElementIterator<2>;
template class MeshSubsetElementIterator<3>;

template  std::ostream &operator<<(std::ostream &os, const MeshSubsetElementIterator<1> & iter);
template  std::ostream &operator<<(std::ostream &os, const MeshSubsetElementIterator<2> & iter);
template  std::ostream &operator<<(std::ostream &os, const MeshSubsetElementIterator<3> & iter);

template class VertexInElemIterator<1>;
template class VertexInElemIterator<2>;
template class VertexInElemIterator<3>;

template  std::ostream &operator<<(std::ostream &os, const VertexInElemIterator<1> & iter);
template  std::ostream &operator<<(std::ostream &os, const VertexInElemIterator<2> & iter);
template  std::ostream &operator<<(std::ostream &os, const VertexInElemIterator<3> & iter);

} // namespace Tortoise
