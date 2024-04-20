//
//  BlockTensor.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 7/11/22.
//
// The structure of a BlockTensor, for instance in 2D,
//  x \ y | 0  1  2  3  4  5  6  7
// ----------------------------------
//    0   | 0  1  1  0  1  1  1  0
//    1   | 0  1  1  0  0  1  1  0
//    2   | 0  0  0  0  0  0  0  0
//    3   | 1  1  1  1  1  0  0  0
//    4   | 1  1  1  0  1  1  1  1
//    5   | 0  0  0  0  0  0  0  0
//
// is saved as a number of horizontal stripes of stored length
// 0 :
//      1  : (2, 0)
//      4  : (3, 2)
// 1 :
//      1  : (2, 5)
//      5  : (2, 7)
// 3 :
//      0  : (5, 9)
// 4 :
//      0  : (3, 14)
//      4  : (4, 17)
//
// The first index is the index of the section row
// The second index is the  index of the section column where the stripe begins
// The first stored integer is the length of the stripe
// The second stored index is the position of the cumulative number of sections until then
// (it will be used to calculate the position of data in vectors storing information on the
// PolygonalMesh, for instance the vector stored in a Function)

#ifndef BlockTensor_hpp
#define BlockTensor_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include <array>
#include <vector>
#include <map>
#include <iostream>
#include <cassert>

namespace Tortoise {

namespace Containers {

template <int NDim> class BlockTensorStructure;
template <> class BlockTensorStructure<1>: public std::map<int,std::array<int,2>> {};
template <> class BlockTensorStructure<2>: public std::map<int,std::map<int,std::array<int,2>>> {};
template <> class BlockTensorStructure<3>: public std::map<int,std::map<int,std::map<int,std::array<int,2>>>>{};


template<int NDim> class BlockTensor {
    // NDim = Number of Dimensions of the tensor
public:
    //=======================================================
    // Constructors
    //===================
    template<typename Derived>
    BlockTensor(const array<int,NDim>& dimensions, const Eigen::MatrixBase<Derived>& whichSections, bool reverseIndexing = false);
 
    template<typename Derived>
    BlockTensor(int dimensions, const Eigen::MatrixBase<Derived>& whichSections, bool reverseIndexing = false) requires (NDim ==1);

    template<typename DerivedA, typename DerivedB>
    BlockTensor(const Eigen::MatrixBase<DerivedA>& dimensions, const Eigen::MatrixBase<DerivedB>& whichSections, bool reverseIndexing = false);
    // whichSections must be a vector with dimension equal to the product of the dimensions.
    // reverseIndexing influences how the entries of the whichSections vector have been consttructed, i.e. which is the
    // outermost index x->y->z or z->y->x
    // Notice that reverseIndexing has no impact on the 1D case
    
    //=======================================================
    // Structure
    //===================
    int size() const;
    int sizeAtDim(int dim) const;
    auto fullTensor() const;

    //=======================================================
    // Explicit Indexing (inefficient)
    //===================
    // returns true if the coordinates are present
    template<typename Derived>
    bool isPresent(const Eigen::MatrixBase<Derived>& coord) const;
    bool isPresent(const array<int,NDim>& coord) const;

    Eigen::Matrix<int, NDim, 1> coord(int id) const;

    // In the following, if debug mode is on, an assert fails if the coord in not present
    // If not in debug mode, then the result is undefined if the coord in not present
    template<typename Derived>
    int id(const Eigen::MatrixBase<Derived>& coord) const;
    int id(const array<int,NDim>& coord) const;
    
    //=======================================================
    // Iterators
    //===================
    class Iterator;
    
    Iterator begin(); // Points at the first element
    Iterator end();   // Points at eh element after the last
    
    //=======================================================
    // I/O
    //===================
    template<int n>
    friend std::ostream& operator<<(std::ostream& os,const BlockTensor<n>& btensor);

    // https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html#StructHavingEigenMembers_othersolutions
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
protected:
    BlockTensorStructure<NDim>   structure;
};

template<int NDim> std::ostream &operator<<(std::ostream& os, BlockTensor<NDim> const& btensor);



//=======================================================
//=======================================================
// Implementation
//=======================================================
//=======================================================


//=======================================================
// Indices and Iterators
//=======================================================
//==========
// NDim = 1
//==========
template<> class BlockTensor<1>::Iterator {
public:
    Iterator(const BlockTensorStructure<1>& structure, const BlockTensorStructure<1>::const_iterator& x, int dx): structure(structure), x(x), dx(dx){};
    Iterator(const BlockTensorStructure<1>& structure, const BlockTensorStructure<1>::const_iterator& x): Iterator(structure, x, 0){};
    
    Iterator& operator++() {
        ++dx;
        if (dx == x->second[0]) {
            dx = 0;
            ++x ;
        }
        return *this;
    }
    Iterator operator++(int){
        Iterator retval = *this;
        ++(*this);
        return retval;
    }
    bool operator==(Iterator other) const {
        return x == other.x && dx == other.dx;
    }
    bool operator!=(Iterator other) const {
        return !(*this == other);
    }
    Eigen::Matrix<int,1,1> operator*() const {
        return Eigen::Matrix<int,1,1> {x->first + dx};
    }
    int id() const {
        return x->second[1]+dx;
    }
    friend std::ostream& operator<<(std::ostream& os,const typename BlockTensor<1>::Iterator& iter);
    
private:
    const BlockTensorStructure<1>&                                      structure;
    BlockTensorStructure<1>::const_iterator                             x;
    int                                                                 dx;
};
std::ostream& operator<<(std::ostream& os,const typename BlockTensor<1>::Iterator& iter){
    os << iter.x->first << " " << iter.dx ;
    return  os;
}
//==========
// NDim = 2
//==========
template<> class BlockTensor<2>::Iterator {
public:
    Iterator(const BlockTensorStructure<2>& structure, const BlockTensorStructure<2>::const_iterator& x, const BlockTensorStructure<2>::value_type::second_type::const_iterator& y, int dy): structure(structure), x(x), y(y), dy(dy){};
    Iterator(const BlockTensorStructure<2>& structure, const BlockTensorStructure<2>::const_iterator& x, const BlockTensorStructure<2>::value_type::second_type::const_iterator& y): Iterator(structure, x,y,0){};
    
    Iterator& operator++() {
        if ( ++dy == y->second[0]) {
            dy = 0;
            if (++y == x->second.end()) {
                if( ++x != structure.end()) {y = x->second.begin();}
            }
        }
        return *this;
    }
    Iterator operator++(int){
        Iterator retval = *this;
        ++(*this);
        return retval;
    }
    bool operator==(Iterator other) const {
        return x == other.x && y == other.y && dy == other.dy;
    }
    bool operator!=(Iterator other) const {
        return !(*this == other);
    }
    Eigen::Matrix<int,2,1> operator*() const {
        return Eigen::Matrix<int,2,1> {x->first, y->first + dy};
    }
    int id() const {
        return y->second[1]+dy;
    }
    friend std::ostream& operator<<(std::ostream& os,const typename BlockTensor<2>::Iterator& iter);
    
private:
    const BlockTensorStructure<2>&                                      structure;
    BlockTensorStructure<2>::const_iterator                             x;
    BlockTensorStructure<2>::value_type::second_type::const_iterator    y;
    int                                                                 dy;
};
std::ostream& operator<<(std::ostream& os,const typename BlockTensor<2>::Iterator& iter){
    os << iter.x->first << " " << iter.y->first << " " << iter.dy ;
    return  os;
}
//==========
// NDim = 3
//==========
template<> class BlockTensor<3>::Iterator {
public:
    Iterator(const BlockTensorStructure<3>& structure, const BlockTensorStructure<3>::const_iterator& x, const BlockTensorStructure<3>::value_type::second_type::const_iterator& y, const BlockTensorStructure<3>::value_type::second_type::value_type::second_type::const_iterator& z, int dz): structure(structure), x(x), y(y), z(z), dz(dz){};
    Iterator(const BlockTensorStructure<3>& structure, const BlockTensorStructure<3>::const_iterator& x, const BlockTensorStructure<3>::value_type::second_type::const_iterator& y, const BlockTensorStructure<3>::value_type::second_type::value_type::second_type::const_iterator& z): Iterator(structure, x,y,z,0){};
    
    Iterator& operator++() {
        if (++dz == z->second[0]) {
            dz = 0;
            if (++z == y->second.end()) {
                if( ++y != x->second.end()) {
                    z = y->second.begin();
                } else {
                    if( ++x != structure.end()) {
                        y = x->second.begin();
                        z = y->second.begin();
                    }
                }
            }
        }
        return *this;
    }
    Iterator operator++(int){
        Iterator retval = *this;
        ++(*this);
        return retval;
    }
    bool operator==(Iterator other) const {
        return x == other.x && y == other.y && z == other.z && dz == other.dz;
    }
    bool operator!=(Iterator other) const {
        return !(*this == other);
    }
    Eigen::Matrix<int,3,1> operator*() const {
        return Eigen::Matrix<int,3,1> {x->first, y->first, z->first + dz};
    }
    int id() const {
        return z->second[1]+dz;
    }
    friend std::ostream& operator<<(std::ostream& os,const typename BlockTensor<3>::Iterator& iter);
    
private:
    const BlockTensorStructure<3>&                                                              structure;
    BlockTensorStructure<3>::const_iterator                                                     x;
    BlockTensorStructure<3>::value_type::second_type::const_iterator                            y;
    BlockTensorStructure<3>::value_type::second_type::value_type::second_type::const_iterator   z;
    int                                                                                         dz;
};

std::ostream& operator<<(std::ostream& os,const typename BlockTensor<3>::Iterator& iter){
    os << iter.x->first << " " << iter.y->first << " " << iter.z->first << " " << iter.dz ;
    return  os;
}


//=======================================================
// Constructors
//===================
//==========
// NDim = 1
//==========
template<int NDim> template<typename Derived>
BlockTensor<NDim>::BlockTensor(int dimensions, const Eigen::MatrixBase<Derived>& whichSections, bool reverseIndexing) requires (NDim ==1) :
BlockTensor(Eigen::Matrix<int,1,1>{dimensions}, whichSections, reverseIndexing){}
template<> template<typename Derived> BlockTensor<1>::BlockTensor(const array<int,1>& dimensions, const Eigen::MatrixBase<Derived>& whichSections, bool reverseIndexing) :
BlockTensor(Eigen::Matrix<int,1,1>{dimensions[0]}, whichSections, reverseIndexing){}
template<> template<typename DerivedA, typename DerivedB> BlockTensor<1>::BlockTensor(const Eigen::MatrixBase<DerivedA>& dimensions, const Eigen::MatrixBase<DerivedB>& whichSections, bool reverseIndexing){
    assert(dimensions.size() == 1);
    assert(whichSections.size() == dimensions[0]);
    int accumulated = 0;
    int j = 0;
    while (j < dimensions[0]) {
        while ( (j < dimensions[0]) && (whichSections(j) == false)) {
            ++j;
        }
        int length = 1;
        if ( j != dimensions[0]) { // Then a 1 has been found
            while ( ((j+length) < dimensions[0]) && (whichSections(j+length) != false)) { ++length;}
            structure[j] = {length, accumulated };
            accumulated += length;
        }
        j += length;
    }
    if (structure.size() == 0) {
        structure[0] = {0,0};
    }
}
//==========
// NDim = 2
//==========
template<> template<typename Derived> BlockTensor<2>::BlockTensor(const array<int,2>& dimensions, const Eigen::MatrixBase<Derived>& whichSections, bool reverseIndexing) :
BlockTensor(Eigen::Matrix<int,2,1>{dimensions[0],dimensions[1]}, whichSections, reverseIndexing){}
template<> template<typename DerivedA, typename DerivedB> BlockTensor<2>::BlockTensor(const Eigen::MatrixBase<DerivedA>& dimensions, const Eigen::MatrixBase<DerivedB>& whichSections, bool reverseIndexing){
    assert(dimensions.size() == 2);
    assert(whichSections.size() == dimensions[0]*dimensions[1]);
    int accumulated = 0;
    for(int i = 0; i < dimensions[0]; ++i){
        int j = 0;
        while (j < dimensions[1]) {
            if (reverseIndexing){
                while ( (j < dimensions[1]) && (whichSections(i+j*dimensions[0]) == false)) {
                    ++j;
                }
            } else {
                while ( (j < dimensions[1]) && (whichSections(i*dimensions[1]+j) == false)) {
                    ++j;
                }
            }
            int length = 1;
            if ( j != dimensions[1]) { // Then a 1 has been found
                if (reverseIndexing){
                    while ( ((j+length) < dimensions[1]) && (whichSections(i+(j+length)*dimensions[0]) != false)) { ++length;}
                } else {
                    while ( ((j+length) < dimensions[1]) && (whichSections(i*dimensions[1]+(j+length)) != false)) { ++length;}
                }
                structure[i][j] = {length, accumulated };
                accumulated += length;
            }
            j += length;
        }
    }
    if (structure.size() == 0) {
        structure[0][0] = {0,0};
    }
}
//==========
// NDim = 3
//==========
template<> template<typename Derived> BlockTensor<3>::BlockTensor(const array<int,3>& dimensions, const Eigen::MatrixBase<Derived>& whichSections, bool reverseIndexing) :
BlockTensor(Eigen::Matrix<int,3,1>{dimensions[0],dimensions[1],dimensions[2]}, whichSections, reverseIndexing){}
template<> template<typename DerivedA, typename DerivedB> BlockTensor<3>::BlockTensor(const Eigen::MatrixBase<DerivedA>& dimensions, const Eigen::MatrixBase<DerivedB>& whichSections, bool reverseIndexing){
    assert(dimensions.size() == 3);
    assert(whichSections.size() == dimensions[0]*dimensions[1]*dimensions[2]);
    int accumulated = 0;
    for(int i = 0; i < dimensions[0]; ++i){
        for(int j = 0; j < dimensions[1]; ++j){
            int l = 0;
            while (l < dimensions[2]) {
                if (reverseIndexing){
                    while ( (l < dimensions[2]) && (whichSections(i+j*dimensions[0]+l*dimensions[0]*dimensions[1]) == false)) {
                        ++l;
                    }
                } else {
                    while ( (l < dimensions[2]) && (whichSections(i*dimensions[2]*dimensions[1]+j*dimensions[2]+l) == false)) {
                        ++l;
                    }
                }
                int length = 1;
                if ( l != dimensions[2]) { // Then a 1 has been found
                    if (reverseIndexing){
                        while ( ((l+length) < dimensions[2]) && (whichSections(i+j*dimensions[0]+(l+length)*dimensions[0]*dimensions[1]) != false)) {
                            ++length;
                        }
                    } else {
                        while ( ((l+length) < dimensions[2]) && (whichSections(i*dimensions[2]*dimensions[1]+j*dimensions[2]+(l+length)) != false)) {
                            ++length;
                        }
                    }
                    structure[i][j][l] = {length, accumulated };
                    accumulated += length;
                }
                l += length;
            }
        }
    }
    if (structure.size() == 0) {
        structure[0][0][0] = {0,0};
    }
}

//=======================================================
// Structure
//===================
//==========
// NDim = 1
//==========
template<> int BlockTensor<1>::size() const {
    if (structure.begin() != structure.end()){
        auto   lastPoint = structure.rbegin();
        return lastPoint->second[0] + lastPoint->second[1];
    }
    return 0;
}
template<> auto BlockTensor<1>::fullTensor() const{
    auto last = structure.rbegin();
    int maxx = last->first+last->second[0];
    if (size() == 0) { maxx = 1;}
    Eigen::Matrix<int, Eigen::Dynamic, 1> matrix = Eigen::Matrix<int, Eigen::Dynamic, 1>::Zero(maxx,1) ;
    for (const auto& [x, info] : structure) {
        for (int i = x; i < x + info[0]; ++i){
            matrix(i,0) = 1;
        }
    }
    return matrix;
}
//==========
// NDim = 2
//==========
template<> int BlockTensor<2>::size() const {
    if (structure.begin() != structure.end()){
        auto   lastPoint = structure.rbegin()->second.rbegin();
        return lastPoint->second[0] + lastPoint->second[1];
    }
    return 0;
}
template<> auto BlockTensor<2>::fullTensor() const{
    int maxy = 0;
    for (const auto& [x, ys] : structure) {
        auto last = ys.rbegin();
        int maxInStripe = last->first + last->second[0];
        if (maxInStripe > maxy) { maxy = maxInStripe;}
    }
    int maxx = size() != 0 ? structure.rbegin()->first + 1 : 0;
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> matrix = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>::Zero(maxx, maxy) ;
    for (const auto& [x, ys] : structure) {
        for (const auto& [y, info] : ys) {
            for (int i = y; i < y + info[0]; ++i){
                matrix(x,i) = 1;
            }
        }
    }
    return matrix;
}
//==========
// NDim = 3
//==========
template<> int BlockTensor<3>::size() const {
    if (structure.begin() != structure.end()){
        auto   lastPoint = structure.rbegin()->second.rbegin()->second.rbegin();
        return lastPoint->second[0] + lastPoint->second[1];
    }
    return 0;
}
template<> auto BlockTensor<3>::fullTensor() const{
    int maxz = 0;
    int maxy = 0;
    for (const auto& [x, ys] : structure) {
        for (const auto& [y, zs] : ys) {
            auto last = zs.rbegin();
            int maxInStripe = last->first + last->second[0];
            if (maxInStripe > maxz) { maxz = maxInStripe;}
        }
        int maxInStripe = ys.rbegin()->first + 1;
        if (maxInStripe > maxy) { maxy = maxInStripe;}
    }
    int maxx = size() != 0 ? structure.rbegin()->first + 1 : 0;
    std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>> tensor;
    for (int i=0; i< maxx; ++i){
        tensor.emplace_back(Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>::Zero(maxy, maxz));
    }
    
    for (const auto& [x, ys] : structure) {
        for (const auto& [y, zs] : ys) {
            for (const auto& [z, info] : zs) {
                for (int i = z; i < z + info[0]; ++i){
                    tensor[x](y,i) = 1;
                }
            }
        }
    }
    return tensor;
}
//=======================================================
// Explicit Indexing (inefficient)
//===================
//==========
// NDim = 1
//==========
template<> template<typename Derived>
bool BlockTensor<1>::isPresent(const Eigen::MatrixBase<Derived>& coord) const{
    assert(coord.size() == 1);
    // The logic below is meant to find the key which is either equal or lower than the passed coodinate
    // Since map does not directly provide that function we need to use upper_bound and then eventually decrement
    auto x = structure.upper_bound(coord[0]);
    if (x == structure.begin()) return false;
    --x;
    return x->first+x->second[0] > coord[0];
}
template<>
bool BlockTensor<1>::isPresent(const array<int,1>& coord) const {
    return isPresent(Eigen::Matrix<int, 1,1>{coord[0]});
}
template<>
Eigen::Matrix<int, 1,1> BlockTensor<1>::coord(int id) const {
    assert(id < size());
    for (const auto& [x, info] : structure) {
        if (id < info[0]) {
            return Eigen::Matrix<int,1,1> {x+id};
        }
        else id -= info[0];
    }
    return Eigen::Matrix<int,1,1> {-1};
}
template<> template<typename Derived>
int BlockTensor<1>::id(const Eigen::MatrixBase<Derived>& coord) const{
    assert(isPresent(coord));
    assert(coord.size() == 1);
    // The logic below is meant to find the key which is either equal or lower than the passed coodinate
    // Since map does not directly provide that function we need to use upper_bound and then eventually decrement
    auto x = structure.upper_bound(coord[0]);
    if (x == structure.begin()) return -1;
    --x;
    return x->second[1]+coord[0]-x->first;
}
template<> int BlockTensor<1>::id(const array<int,1>& coord) const {
    return id(Eigen::Matrix<int, 1,1>{coord[0]});
}
//==========
// NDim = 2
//==========
template<> template<typename Derived>
bool BlockTensor<2>::isPresent(const Eigen::MatrixBase<Derived>& coord) const{
    assert(coord.size() == 2);
    auto x = structure.find(coord[0]);
    if (x == structure.end()) return false;
    // The logic below is meant to find the key which is either equal or lower than the passed coodinate
    // Since map does not directly provide that function we need to use upper_bound and then eventually decrement
    auto y = x->second.upper_bound(coord[1]);
    if (y == x->second.begin()) return false;
    --y;
    return y->first+y->second[0] > coord[1];
}
template<>
bool BlockTensor<2>::isPresent(const array<int,2>& coord) const {
    return isPresent(Eigen::Matrix<int, 2,1>{coord[0],coord[1]});
}
template<>
Eigen::Matrix<int, 2,1> BlockTensor<2>::coord(int id) const {
    assert(id < size());
    for (const auto& [x, ys] : structure) {
        for (const auto& [y, info] : ys) {
            if (id < info[0]) {return Eigen::Matrix<int,2,1> {x, y+id }; }
            else id -= info[0];
        }
    }
    return Eigen::Matrix<int,2,1> {-1, -1};
}
template<> template<typename Derived>
int BlockTensor<2>::id(const Eigen::MatrixBase<Derived>& coord) const{
    assert(isPresent(coord));
    assert(coord.size() == 2);
    auto x = structure.find(coord[0]);
    if (x == structure.end()) return -1;
    // The logic below is meant to find the key which is either equal or lower than the passed coodinate
    // Since map does not directly provide that function we need to use upper_bound and then eventually decrement
    auto y = x->second.upper_bound(coord[1]);
    if (y == x->second.begin()) return -1;
    --y;
    return y->second[1]+coord[1]-y->first;
}
template<> int BlockTensor<2>::id(const array<int,2>& coord) const {
    return id(Eigen::Matrix<int, 2,1>{coord[0],coord[1]});
}
//==========
// NDim = 3
//==========
template<> template<typename Derived>
bool BlockTensor<3>::isPresent(const Eigen::MatrixBase<Derived>& coord) const{
    assert(coord.size() == 3);
    auto x = structure.find(coord[0]);
    if (x == structure.end()) return false;
    auto y = x->second.find(coord[1]);
    if (y == x->second.end()) return false;
    // The logic below is meant to find the key which is either equal or lower than the passed coodinate
    // Since map does not directly provide that function we need to use upper_bound and then eventually decrement
    auto z = y->second.upper_bound(coord[2]);
    if (z == y->second.begin()) return false;
    --z;
    return z->first+z->second[0] > coord[2];
}
template<>
bool BlockTensor<3>::isPresent(const array<int,3>& coord) const {
    return isPresent(Eigen::Matrix<int, 3,1>{coord[0],coord[1],coord[2]});
}
template<>
Eigen::Matrix<int, 3,1> BlockTensor<3>::coord(int id) const {
    assert(id < size());
    for (const auto& [x, ys] : structure) {
        for (const auto& [y, zs] : ys) {
            for (const auto& [z, info] : zs) {
                if (id < info[0]) {return Eigen::Matrix<int,3,1> {x, y, z+id }; }
                else id -= info[0];
            }
        }
    }
    return Eigen::Matrix<int,3,1> {-1, -1, -1};
}
template<> template<typename Derived>
int BlockTensor<3>::id(const Eigen::MatrixBase<Derived>& coord) const{
    assert(isPresent(coord));
    assert(coord.size() == 3);
    auto x = structure.find(coord[0]);
    if (x == structure.end()) return -1;
    auto y = x->second.find(coord[1]);
    if (y == x->second.end()) return -1;
    // The logic below is meant to find the key which is either equal or lower than the passed coodinate
    // Since map does not directly provide that function we need to use upper_bound and then eventually decrement
    auto z = y->second.upper_bound(coord[2]);
    if (z == y->second.begin()) return -1;
    --z;
    return z->second[1]+coord[2]-z->first;
}
template<> int BlockTensor<3>::id(const array<int,3>& coord) const {
    return id(Eigen::Matrix<int, 3,1>{coord[0],coord[1],coord[2]});
}
//=======================================================
// Iterators
//===================
//==========
// NDim = 1
//==========
template<> BlockTensor<1>::Iterator BlockTensor<1>::begin() {
    return Iterator(structure, structure.begin());
}
template<> BlockTensor<1>::Iterator BlockTensor<1>::end() {
    if (size() == 0){
        return Iterator(structure, structure.begin());
    }
    return Iterator(structure, structure.end());
}
//==========
// NDim = 2
//==========
template<> BlockTensor<2>::Iterator BlockTensor<2>::begin() {
    return Iterator(structure, structure.begin(), structure.begin()->second.begin());
}
template<> BlockTensor<2>::Iterator BlockTensor<2>::end() {
    if (size() == 0){
        return Iterator(structure, structure.begin(), structure.begin()->second.begin());
    }
    return Iterator(structure, structure.end(), structure.rbegin()->second.end());
}
//==========
// NDim = 3
//==========
template<> BlockTensor<3>::Iterator BlockTensor<3>::begin() {
    return Iterator(structure, structure.begin(), structure.begin()->second.begin(), structure.begin()->second.begin()->second.begin());
}
template<> BlockTensor<3>::Iterator BlockTensor<3>::end() {
    if (size() == 0){
        return Iterator(structure, structure.begin(), structure.begin()->second.begin(), structure.begin()->second.begin()->second.begin());
    }
    return Iterator(structure, structure.end(), structure.rbegin()->second.end(), structure.rbegin()->second.rbegin()->second.end());
}

//=======================================================
// I/O
//===================
//==========
// NDim = 1
//==========
template<> std::ostream& operator<< <1>(std::ostream& os, BlockTensor<1> const& btensor){
    for (const auto& [x, info] : btensor.structure) {
        os << "   " << x << " : (" << info[0] << ", " << info[1] << ")\n"  ;
    }
    return  os;
}
//==========
// NDim = 2
//==========
template<> std::ostream& operator<< <2>(std::ostream& os, BlockTensor<2> const& btensor){
    for (const auto& [x, ys] : btensor.structure) {
        os << x << " :\n";
        for (const auto& [y, info] : ys) {
            os << "   " << y << " : (" << info[0] << ", " << info[1] << ")\n"  ;
        }
    }
    return  os;
}
//==========
// NDim = 3
//==========
template<> std::ostream& operator<< <3>(std::ostream& os, BlockTensor<3> const& btensor){
    for (const auto& [x, ys] : btensor.structure) {
        os << x << " :\n";
        for (const auto& [y, zs] : ys) {
            os << "   " << y << " :\n";
            for (const auto& [z, info] : zs) {
                os << "      " << z << " : (" << info[0] << ", " << info[1] << ")\n"  ;
            }
        }
    }
    return  os;
}


} // namespace Containers

} // namespace Tortoise


#endif /* BlockTensor_hpp */
    
