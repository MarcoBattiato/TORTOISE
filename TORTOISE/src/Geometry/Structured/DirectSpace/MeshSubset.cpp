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
//  MeshSubset.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 23/6/21.
//

#include <Geometry/Structured/DirectSpace/MeshSubset.hpp>
#include <Generics/Utilities/RandomNGenerator.hpp>

#include <cmath>        // floor
#include <algorithm>    // std::max

namespace Tortoise {

//=======================================================
// Constructors
//===================
// Subset covering the whole mesh
template<int NDim>
MeshSubset<NDim>::MeshSubset(const Mesh<NDim>& t_mesh)
: mesh(&t_mesh), lSecCoord(CartIndex<NDim>::Zero()), hSecCoord(t_mesh.nSecSplits){}

template<int NDim>
MeshSubset<NDim>::MeshSubset(const Mesh<NDim>* t_mesh)
: MeshSubset(*t_mesh){}


// Subset with specified coverage
template<int NDim>
MeshSubset<NDim>::MeshSubset(const Mesh<NDim>& t_mesh, const CartIndex<NDim>& t_lSecCoord, const CartIndex<NDim>& t_hSecCoord)
: mesh(&t_mesh), lSecCoord(t_lSecCoord), hSecCoord(t_hSecCoord){
    assert( (t_lSecCoord.array() >= 0).all() && (t_lSecCoord.array() <= mesh->nSecSplits.array()).all() );
    assert( (t_hSecCoord.array() >= 0).all() && (t_hSecCoord.array() <= mesh->nSecSplits.array()).all() );
    assert( (t_hSecCoord.array() <= t_hSecCoord.array()).all() );
};

template<int NDim>
MeshSubset<NDim>::MeshSubset(const Mesh<NDim>* t_mesh, const CartIndex<NDim>& t_lSecCoord, const CartIndex<NDim>& t_hSecCoord)
: MeshSubset(*t_mesh, t_lSecCoord, t_hSecCoord){};


// Subset covering the section pointed at by the SectionIterator
template<int NDim>
MeshSubset<NDim>::MeshSubset(const Mesh<NDim>& t_mesh, const SectionIterator<NDim>& section)
: mesh(&t_mesh), lSecCoord(section.currentcoord), hSecCoord(section.currentcoord + CartIndex<NDim>::Ones()){
    assert(section.unfinished);
}

template<int NDim>
MeshSubset<NDim>::MeshSubset(const Mesh<NDim>& t_mesh, const CartIndex<NDim>& section)
: mesh(&t_mesh), lSecCoord(section), hSecCoord(section + CartIndex<NDim>::Ones()){}

template<int NDim>
MeshSubset<NDim>::MeshSubset(const Mesh<NDim>* t_mesh, const SectionIterator<NDim>& section)
: MeshSubset(*t_mesh, section){}

template<int NDim>
MeshSubset<NDim>::MeshSubset(const Mesh<NDim>* t_mesh, const CartIndex<NDim>& section)
: MeshSubset(*t_mesh, section){}

template<int NDim>
MeshSubset<NDim>::MeshSubset(const Mesh<NDim>& t_mesh, const TrapezoidalSubset<NDim>& subset)
: mesh(&t_mesh) {
    lSecCoord = (t_mesh.invdGVec*(subset.lcoord-t_mesh.origin)).array().floor().template cast<int>().max(0).min(t_mesh.nSecSplits.array());
    hSecCoord = ((t_mesh.invdGVec*(subset.hcoord-t_mesh.origin)).array().floor().template cast<int>() + 1).max(0).min(t_mesh.nSecSplits.array());
}

template<int NDim>
MeshSubset<NDim>::MeshSubset(const Mesh<NDim>* t_mesh, const TrapezoidalSubset<NDim>& subset)
: MeshSubset(*t_mesh, subset){}

//=======================================================
// Converter to UnitSquareSubset
//===================
template<int NDim>
TrapezoidalSubset<NDim> MeshSubset<NDim>::trapezoidalSubset() const {
    assert((hSecCoord.array()>=lSecCoord.array()).all()
           && "This MeshSubset has low and high SecCoord inverted, this should have not happened! Please let the developer know.");
    return TrapezoidalSubset<NDim>(mesh->origin + mesh->dGVec * (lSecCoord.template cast<Real>().array() + 0.00000001).matrix(), mesh->origin + mesh->dGVec * (hSecCoord.template cast<Real>().array() - 0.00000001).matrix());
}

//=======================================================
// Iterators
//===================
template<int NDim>
SectionIterator<NDim> MeshSubset<NDim>::sectionIterator() const {
    return SectionIterator<NDim>(mesh->nSecSplits, lSecCoord, hSecCoord);
}
template<int NDim>
MeshSubsetElementIterator<NDim> MeshSubset<NDim>::elementIterator() const {
    return MeshSubsetElementIterator<NDim>(mesh->nSecSplits, lSecCoord, hSecCoord);
}

template<int NDim>
int MeshSubset<NDim>::numberCoveredSections() const {
    return (hSecCoord-lSecCoord).prod();
}

//********************************
//* Arithmetic
//********************************

template<int NDim>
MeshSubset<NDim>& MeshSubset<NDim>::operator= (const MeshSubset<NDim>& other) {
    assert(mesh == other.mesh);
    lSecCoord = other.lSecCoord;
    hSecCoord = other.hSecCoord;
    return *this;
}

template<int NDim>
MeshSubset<NDim>& MeshSubset<NDim>::operator= (const TrapezoidalSubset<NDim>& subset) {
    lSecCoord = (mesh->invdGVec*(subset.lcoord-mesh->origin)).array().floor().template cast<int>().max(0).min(mesh->nSecSplits.array());
    hSecCoord = ((mesh->invdGVec*(subset.hcoord-mesh->origin)).array().floor() + 1).template cast<int>().max(0).min(mesh->nSecSplits.array());
    return *this;
}
                 
//=======================================================
// I/O
//===================

Plotter3D& addsquaretoplot(Plotter3D& plotter, Point<2>& corner, Point<2>& side1, Point<2>& side2) {
    if(plotter.numberplots++>0) plotter.plotCommand +=" , ";
    plotter.plotCommand += " '-' using 1:2:(0) with lines ";
    Point<2> p = corner;    plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
    p += side1;             plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
    p += side2;             plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
    p -= side1;             plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
    p -= side2;             plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
    plotter.plotData += "e\n";
    return plotter << COLOR("green") << NOLABEL;
}
Plotter3D& operator << (Plotter3D& plotter, const MeshSubset<2>& subset){
//    plotter << *(subset.mesh->region) << NOLABEL << COLOR("black") ;
    plotter << *subset.mesh << NOLABEL << COLOR("#BBFF0000") ;
        Point<2> corner = subset.mesh->origin + subset.mesh->dGVec*subset.lSecCoord.template cast<Real>();
        Point<2> side1  = subset.mesh->dGVec.col(0)*(subset.hSecCoord.template cast<Real>()(0)-subset.lSecCoord.template cast<Real>()(0));
        Point<2> side2  = subset.mesh->dGVec.col(1)*(subset.hSecCoord.template cast<Real>()(1)-subset.lSecCoord.template cast<Real>()(1));
        addsquaretoplot(plotter,corner,side1,side2);
    return plotter;
}

template<int NDim>
std::ostream &operator<<(std::ostream &os, MeshSubset<NDim> const& meshSubset){
    os << "[" << meshSubset.lSecCoord.transpose() << "][" << meshSubset.hSecCoord.transpose()<<"]" ;
    return os;
}
template<int NDim>
std::ostream &operator<<(std::ostream &os, TrapezoidalSubset<NDim> const& meshSubset){
    os << "[" << meshSubset.lcoord.transpose() << "][" << meshSubset.hcoord.transpose()<<"]" ;
    return os;
}


//=======================================================
//=======================================================
// UnitSquareSubset
//===================
//===================

template <int NDim>
TrapezoidalSubset<NDim>::TrapezoidalSubset(const Point<NDim>& t_lcoord, const Point<NDim>& t_hcoord)
: lcoord(t_lcoord), hcoord(t_hcoord){}

template <int NDim>
TrapezoidalSubset<NDim> TrapezoidalSubset<NDim>::operator-() const {
    return TrapezoidalSubset(-hcoord,-lcoord);}

template <int NDim>
TrapezoidalSubset<NDim>& TrapezoidalSubset<NDim>::operator=(TrapezoidalSubset<NDim> other) {
    lcoord.swap(other.lcoord);
    hcoord.swap(other.hcoord);
    return *this;}

template <int NDim>
TrapezoidalSubset<NDim>& TrapezoidalSubset<NDim>::operator+=(const TrapezoidalSubset<NDim>& other) {
    lcoord += other.lcoord;
    hcoord += other.hcoord;
    return *this;}

template <int NDim>
TrapezoidalSubset<NDim>& TrapezoidalSubset<NDim>::operator-=(const TrapezoidalSubset<NDim>& other) {
    lcoord -= other.hcoord;
    hcoord -= other.lcoord;
    return *this;}

template <int NDim>
TrapezoidalSubset<NDim>& TrapezoidalSubset<NDim>::operator*=(const Real scalar) {
    if (scalar>=0.0)  {lcoord *= scalar; hcoord *= scalar;}
    else {lcoord.swap(hcoord); lcoord *= scalar; hcoord *= scalar;}
    return *this;}

template <int NDim>
TrapezoidalSubset<NDim>& TrapezoidalSubset<NDim>::operator+=(const Point<NDim>& other){
    lcoord += other;
    hcoord += other;
    return *this;}

template <int NDim>
TrapezoidalSubset<NDim>& TrapezoidalSubset<NDim>::operator-=(const Point<NDim>& other) {
    lcoord -= other;
    hcoord -= other;
    return *this;}

  
//=======================================================
// Debug utilities
//===================
//template <int NDim> class UnitSquareSubsetToPlot {
//public:
//    const Region<NDim>*                 region;
//    const UnitSquareSubset<NDim>*       unitSquare;
//    UnitSquareSubsetToPlot(const Region<NDim>& t_region, const UnitSquareSubset<NDim>& t_unitSquare): region(&t_region), unitSquare(&t_unitSquare){};
//};
//
//Plotter3D& operator << (Plotter3D& plotter, const UnitSquareSubsetToPlot<2>& t_subsetToPlot){
//    plotter << *(t_subsetToPlot.region) << NOLABEL << COLOR("black") ;
//    if(plotter.numberplots++>0) plotter.plotCommand +=" , ";
//    plotter.plotCommand += " '-' using 1:2:(0) with lines ";
//    Point<2> p(t_subsetToPlot.region->origin + t_subsetToPlot.region->gVec * (t_subsetToPlot.unitSquare->lcoord));
//    plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
//    p += t_subsetToPlot.region->gVec.col(0) * (t_subsetToPlot.unitSquare->hcoord(0)-t_subsetToPlot.unitSquare->lcoord(0));
//    plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
//    p += t_subsetToPlot.region->gVec.col(1) * (t_subsetToPlot.unitSquare->hcoord(1)-t_subsetToPlot.unitSquare->lcoord(1));
//    plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
//    p -= t_subsetToPlot.region->gVec.col(0) * (t_subsetToPlot.unitSquare->hcoord(0)-t_subsetToPlot.unitSquare->lcoord(0));
//    plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
//    p -= t_subsetToPlot.region->gVec.col(1) * (t_subsetToPlot.unitSquare->hcoord(1)-t_subsetToPlot.unitSquare->lcoord(1));
//    plotter.plotData += std::to_string(p(0)) + " " + std::to_string(p(1)) + "\n";
//    plotter.plotData += "e\n";
//    return plotter;
//}




// ===========================================
// ===========================================
// ===========================================
// Explicit template instantiation (For more info see: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl)

template class MeshSubset<1>;
template class MeshSubset<2>;
template class MeshSubset<3>;

template std::ostream &operator<<(std::ostream &os, MeshSubset<1> const& meshSubset);
template std::ostream &operator<<(std::ostream &os, MeshSubset<2> const& meshSubset);
template std::ostream &operator<<(std::ostream &os, MeshSubset<3> const& meshSubset);

template class TrapezoidalSubset<1>;
template class TrapezoidalSubset<2>;
template class TrapezoidalSubset<3>;

template std::ostream &operator<<(std::ostream &os, TrapezoidalSubset<1> const& meshSubset);
template std::ostream &operator<<(std::ostream &os, TrapezoidalSubset<2> const& meshSubset);
template std::ostream &operator<<(std::ostream &os, TrapezoidalSubset<3> const& meshSubset);

} // namespace Tortoise
