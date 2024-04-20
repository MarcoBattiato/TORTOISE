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
//  Function.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//
//
//  TO ADD:
//      - Write to File
//      - Convolution
//      - Change mesh
//
//  TO ADD LATER:
//      - Derivative

#ifndef Function_hpp
#define Function_hpp

#include <Geometry/Structured/FunctionSpace/FunctionElement.hpp>
#include <Generics/Features/VectorSpace.hpp>
#include <Geometry/Structured/DirectSpace/MeshSubset.hpp>
#include <Geometry/Structured/DirectSpace/MeshIterators.hpp>
#include <Geometry/Plotter/Plotter.hpp>

#include <functional>
#include <math.h>
#include <concepts>

namespace Tortoise {

template <typename functType> concept boolValuedFunctOfReal = requires (functType f, Real val) { {f(val)} -> std::same_as<bool>;};
enum RequirementType {any = -1, all = 1};

//*******************************
// Main class definition
//*******************************

template<int NDim> class Function :
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Inheritance
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
// Inherits all the operations of a field from MathFieldSpace and AsymmetricMathFieldSpace
// Math field operations are defined (+,-,*,/).
// The same operations are defined with scalars.
// Only + and - are defined on FunctionElement. This is because it would be unclear what the user might expect to happen on the remaining
// elements. The algebraic operations are done assuming that the FunctionElement is non zero only on the specified element. However if the same
// interpretation were to be used for * and /, problems arise (especially in the case of /). If the user wants to do those operations he/she should first
// construct a Function from the FunctionElement and explicitly set the value everywhere else to avoid confusion.
// This inheritance constructs all the possible combinations of operators
public Features::MathFieldSpace<Function<NDim>> ,
public Features::AsymmetricMathFieldSpaceByValue<Function<NDim>,Real> ,
//    public Features::AsymmetricMathFieldSpace<FunctionElement<NDim>, std::function<Real(Point<NDim>)>> ,
public Features::AsymmetricVectorSpace< Function<NDim>, FunctionElement<NDim>> ,
public Features::AsymmetricMathFieldSpace< Function<NDim>, std::function<Real(Point<NDim>)>>{

        
        
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main properties
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    Mesh<NDim> const * const    mesh;
    DataVector                  vec;
    
    static GeometryCore::LinTransform<NDim> const  fromLinFormToModalTransf;
    static GeometryCore::LinTransform<NDim> const  fromModalToLinFormTransf;
        
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Main methods
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~

public:
    //=======================================================
    // Constructors
    //===================
    explicit Function(const Mesh<NDim>& t_mesh);                                                    // Constructs a constant 0 function
    Function(const Mesh<NDim>& t_mesh, const Real t_value);                                         // Constructs a constant function
    Function(const Mesh<NDim>& t_mesh, const DataVector& t_nodalValues);                            // Constructs using the values at nodes
    Function(const Mesh<NDim>& t_mesh, const std::function<Real(Point<NDim>)>& t_f);                // Constructs using the function passed
    explicit Function(const Mesh<NDim>* t_mesh);                                                    // Constructs a constant 0 function
    Function(const Mesh<NDim>* t_mesh, const Real t_value);                                         // Constructs a constant function
    Function(const Mesh<NDim>* t_mesh, const DataVector& t_nodalValues);                            // Constructs using the values at nodes
    Function(const Mesh<NDim>* t_mesh, const std::function<Real(Point<NDim>)>& t_f);                // Constructs using the function passed
    Function(const std::function<Real(Real)>& t_f, const Function<NDim>& t_other);                  // Constructs a function of a function f(other)
    Function(const std::function<Real(std::vector<Real>)>& t_f, const std::vector<Function<NDim>>& t_others); // Constructs a function of a function f(other[0], other[1], ...)
    Function(const Function<NDim>& t_other);                                                        // Copy constructor
    Function(Function<NDim>&& t_other);                                                             // Move constructor
    
    //********************************
    //* Arithmetic
    //********************************
    // ATTENTION!!!!!!! Multiplications and divisions by anything except a scalar do not conserve precision!
    // If you don't know what this means, do not use those operations. You will get unexpected results in many cases
    Function& operator=(const Function<NDim>& other);
    Function& operator=(Function<NDim>&& other);
    Function& operator+=(const Function<NDim>& other);
    Function& operator-=(const Function<NDim>& other);
    Function& operator*=(const Function<NDim>& other);                          // !!!! Loses precision. DO NOT USE if you do not know what this means!!!!
    Function& operator/=(const Function<NDim>& other);                          // !!!! Loses precision. DO NOT USE if you do not know what this means!!!!

    Function<NDim> operator-() const;
        
    Function& operator=(const Real scalar);
    Function& operator+=(const Real scalar);
    Function& operator-=(const Real scalar);
    Function& operator*=(const Real scalar);
    Function& operator/=(const Real scalar);
//    friend Function<NDim> operator/(const Real lhs, Function<NDim> rhs);
        
    Function& operator+=(const FunctionElement<NDim>& other);
    Function& operator-=(const FunctionElement<NDim>& other);
        
    Function& operator=(const std::function<Real(Point<NDim>)>& other);
    Function& operator+=(const std::function<Real(Point<NDim>)>& other);
    Function& operator-=(const std::function<Real(Point<NDim>)>& other);
    Function& operator*=(const std::function<Real(Point<NDim>)>& other);
    Function& operator/=(const std::function<Real(Point<NDim>)>& other);
//    friend Function<NDim> operator/(const std::function<Real(Point<NDim>)>& lhs, Function<NDim> rhs);
      
    //********************************
    //* Mesh structure
    //********************************
    //=======================================================
    // Mesh Iterators
    // Iterating through the mesh should be done as
    // >>  for (auto elem = functexmpl.elemIterator(); elem.unfinished; ++elem) {
    // >>       something[elem]
    // >>  }
    // This allows for abstraction and also works for MeshSubsets that run only over some of the elements
    // Notice that the iterator is always initialised to the first section, element, or vertex
    //===================
    // Mesh iterators
    SectionIterator<NDim>       sectionIterator() const;
    MeshElementIterator<NDim>   elementIterator() const;
    VertexInElemIterator<NDim>  vertexInElemIterator() const;
        
    // Random accessors
    CartIndex<NDim>             randomSection() const;
    int                         randomElement() const;

    //********************************
    //* Local Evaluation
    //********************************
    FunctionElement<NDim> operator()(const MeshElementIterator<NDim>& t_elementID) const;                        // Function in a given element
    FunctionElement<NDim> operator()(const MeshSubsetElementIterator<NDim>& t_elementID) const;                  // Function in a given element
    FunctionElement<NDim> operator()(const int t_elementID) const;                                               // Function in a given element
        
    Real operator()(const MeshElementIterator<NDim>& t_elementID, const VertexInElemIterator<NDim>& t_nVertex) const;     // Value at a given vertex in a given element
    Real operator()(const MeshElementIterator<NDim>& t_elementID, const Point<NDim>& t_point_ref) const;     // Value at a given point in relative coord in a given element
    Real operator()(const MeshElementIterator<NDim>& t_elementID, const Real t_point_ref) const requires (NDim == 1);     // Value at a given point in relative coord in a given element
    Real operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const VertexInElemIterator<NDim>& t_nVertex) const;     // Value at a given vertex in a given element
    Real operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const Point<NDim>& t_point_ref) const;     // Value at a given point in relative coord in a given element
    Real operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const Real t_point_ref) const requires (NDim == 1);     // Value at a given point in relative coord in a given element
    Real operator()(const int t_elementID, const Point<NDim>& t_point_ref) const;                       // Value at a given point in relative coord in a given element
    Real operator()(const int t_elementID, const Real t_point_ref) const requires (NDim == 1);          // Value at a given point in relative coord in a given element
    Real operator()(const int t_elementID, const int t_nVertex) const;                                  // Value at a given vertex in a given element
    // Lower level methods
    auto elementVec(const MeshElementIterator<NDim> t_elementID) const;                 // part of the vecvtor corresponding to a given element: returns LinearForm<NDim>
    auto elementVec(const MeshSubsetElementIterator<NDim> t_elementID) const;           // part of the vecvtor corresponding to a given element: returns LinearForm<NDim>
    auto elementVec(const int t_elementID) const;                                       // part of the vecvtor corresponding to a given element: returns LinearForm<NDim>
    auto elementVec(const MeshElementIterator<NDim> t_elementID);                       // part of the vecvtor corresponding to a given element: returns LinearForm<NDim>
    auto elementVec(const MeshSubsetElementIterator<NDim> t_elementID);                 // part of the vecvtor corresponding to a given element: returns LinearForm<NDim>
    auto elementVec(const int t_elementID);                                             // part of the vecvtor corresponding to a given element: returns LinearForm<NDim>

    //********************************
    //* Other operations
    //********************************
                
    Real  elementMax(const MeshElementIterator<NDim>& t_elementID) const;                      // Max within an element
    Real  elementMin(const MeshElementIterator<NDim>& t_elementID) const;                      // Max within an element
    Real  elementMax(const MeshSubsetElementIterator<NDim>& t_elementID) const;                // Max within an element
    Real  elementMin(const MeshSubsetElementIterator<NDim>& t_elementID) const;                // Max within an element
    Real  elementMax(const int t_elementID) const;                                             // Max within an element
    Real  elementMin(const int t_elementID) const;                                             // Max within an element
    Real  max() const;                                                                         // Max of the whole function
    Real  min() const;                                                                         // Min of the whole function
    std::pair<std::vector<int>,std::vector<Point<NDim>>> zeros() const requires (NDim == 1);   // Find zeros of a function
        
    Real  integrate() const;                                                                   // Calculates the integral of the function
    Real  integrate(const Function<NDim>& rhs) const;                                          // Calculates the integral of the function multiplied by another one
    Real  integrate(const Function<NDim>& rhs1, const Function<NDim>& rhs2) const;             // Calculates the integral of the function multiplied other two
    Real  integrateDiracDelta(const Function<NDim>& functInDiracDelta) const;                  // Calculates the integral of the function multiplied another Dirac-delta function
    Mesh<NDim>     convolutionMesh() const;                                                    // Returns the mesh over which the associated convolution function should be defined
    Function<NDim> convolve(const Function<NDim>& convolveFunction) const;
            
    Function<NDim> derivative(const int direction) const;                                      // Calculates the derivative
                
    Function&      applyInverseMass();                                                         // Multiplies by inverse mass (it modifies the function!!!)
                    
    // The method below returns an Eigen vector of bool specifying which sections satisfy a certain constraint on any or all of the vertices of the section
//    template <typename functType>
//    auto evaluateBoolCondOnSections(const functType& boolFunctionToEvaluate, const RequirementType requirementType) requires boolValuedFunctOfReal<functType>;
    auto evaluateBoolCondOnSections(const std::function<bool(Real)>& boolFunctionToEvaluate, const RequirementType requirementType) const;
        
    //********************************
    //* I/O
    //********************************
    
    void plot(const std::string& t_title = "") const;
    void plot(const VectorPoint<NDim>& line, const std::vector<std::string>& names, const std::string& t_title = "") const;
    void writeToTxtFile(const std::string &FileName) const;
    // Added by Xavier, Dec 2023 : Begin
    void writeToTxtFilePlot(const std::string &FileName) const; // Creates file to be used by plotting programs like GNUPlot
    void writeToTxtFilePlotPy(std::string &FileNameP)  const; // Creates file to be used by plotting programs like GNUPlot
    void writeToTxtFilePlotAppend(std::string &FileName)  const; // Creates file to be used by plotting programs like GNUPlot
    void writeToTxtFilePlotColour(std::string &FileName)  const;
    // Added by Xavier, Dec 2023 : End
    
    //=======================================================
    // Technicalities
    //=======================================================
public:
    auto elementview();         // Returns a map. Can be used to modify directly vec, but it is accessed in element view
    auto elementview() const ;  // Returns a map. Cannot be used to modify directly vec, since it is constant
    void toNodalFromLinForm() ;      // Modifies the vec (the user should remember to eventually modify it back into linear form representation)
    void toLinFormFromNodal() ;      // Modifies the vec
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> toNodalFromLinForm() const ; // When the object is const it cannot modify the vec, but instead it returns a matrix. It is up to the user to store it
    void toModalFromLinForm();
    void toLinFormFromModal();

    // Remove the possibility of passing temporaries to the constructor
    explicit Function(const Mesh<NDim>&& t_mesh) = delete;
    Function(const Mesh<NDim>&& t_mesh, const Real t_value) = delete;
    Function(const Mesh<NDim>&& t_mesh, const DataVector& t_nodalValues) = delete;
    Function(const Mesh<NDim>&& t_mesh, const std::function<Real(Point<NDim>)>& t_f) = delete;

};

template<int NDim> class Projected{
public:
    Function<NDim> const&  fun;
    VectorPoint<NDim> const& line;
    std::vector<std::string> names;
    Projected(const Function<NDim>& fun, const VectorPoint<NDim>& line, const std::vector<std::string>& names): fun(fun), line(line), names(names){}
    Projected(const Function<NDim>& fun, const VectorPoint<NDim>& line): fun(fun), line(line){
        for (int i=0; i<line.cols(); ++i){ names.emplace_back("");}
    }
};

template<int NDim> class ProjectedShifted{
public:
    Function<NDim> const&  fun;
    VectorPoint<NDim> const& line;
    std::vector<std::string> names;
    Real const shift;
    ProjectedShifted(const Function<NDim>& fun, Real shift, const VectorPoint<NDim>& line, const std::vector<std::string>& names): fun(fun), shift(shift), line(line), names(names){}
    ProjectedShifted(const Function<NDim>& fun, Real shift, const VectorPoint<NDim>& line): fun(fun), shift(shift), line(line){
        for (int i=0; i<line.cols(); ++i){ names.emplace_back("");}
    }
};

Plotter3D& operator << (Plotter3D& plotter, const Function<2>& function);
Plotter2D& operator << (Plotter2D& plotter, const Function<1>& function);
Plotter2D& operator << (Plotter2D& plotter, const Projected<2>& projFunction);
Plotter2D& operator << (Plotter2D& plotter, const ProjectedShifted<2>& projFunction);

template<int NDim> Function<NDim> operator/(const Real lhs, Function<NDim> rhs);
template<int NDim> Function<NDim> operator/(const std::function<Real(Point<NDim>)>& lhs, Function<NDim> rhs);

template<int NDim> void swap(Function<NDim>& first, Function<NDim>& second);                                // Swap functionality

template<int NDim> auto sin(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::sin(z);}, other); }
template<int NDim> auto cos(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::cos(z);}, other); }
template<int NDim> auto tan(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::tan(z);}, other); }
template<int NDim> auto asin(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::asin(z);}, other); }
template<int NDim> auto acos(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::acos(z);}, other); }
template<int NDim> auto atan(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::atan(z);}, other); }
template<int NDim> auto cosh(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::cosh(z);}, other); }
template<int NDim> auto sinh(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::sinh(z);}, other); }
template<int NDim> auto tanh(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::tanh(z);}, other); }
template<int NDim> auto pow(const Function<NDim>& other, Real exponent) {return Function<NDim>([exponent](Real z){return std::pow(z, exponent);}, other); }
template<int NDim> auto sqrt(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::sqrt(z);}, other); }
template<int NDim> auto log(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::log(z);}, other); }
template<int NDim> auto floor(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::floor(z);}, other); }
template<int NDim> auto ceil(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::ceil(z);}, other); }
template<int NDim> auto abs(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::fabs(z);}, other); }
template<int NDim> auto fabs(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::fabs(z);}, other); }
template<int NDim> auto exp(const Function<NDim>& other) { return Function<NDim>([](Real z){return std::exp(z);}, other); }

template<int NDim> auto min(const Function<NDim>& f0, const Function<NDim>& f1) { return Function<NDim>([](const std::vector<Real>& z){return std::min(z[0],z[1]);}, {f0, f1});}
template<int NDim> auto max(const Function<NDim>& f0, const Function<NDim>& f1) { return Function<NDim>([](const std::vector<Real>& z){return std::max(z[0],z[1]);}, {f0, f1});}
template<int NDim> auto max(const Function<NDim>& f0, Real val) { return Function<NDim>([val](const Real z){return std::max(z,val);}, f0);}
template<int NDim> auto max(Real val, const Function<NDim>& f0) { return Function<NDim>([val](const Real z){return std::max(z,val);}, f0);}
template<int NDim> auto min(const Function<NDim>& f0, Real val) { return Function<NDim>([val](const Real z){return std::min(z,val);}, f0);}
template<int NDim> auto min(Real val, const Function<NDim>& f0) { return Function<NDim>([val](const Real z){return std::min(z,val);}, f0);}

//=======================================================
// Implementation Details
//===================

template<int NDim> auto Function<NDim>::elementview() {
    return Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vec.data(),  mesh->numberElements, mesh->nVertElement);
}
template<int NDim> auto Function<NDim>::elementview() const {
    return Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vec.data(), mesh->numberElements ,mesh->nVertElement );
}
template<int NDim> void Function<NDim>::toNodalFromLinForm() {
    auto elementviewvar = elementview();
    elementviewvar.rightCols(elementviewvar.cols()-1).colwise() += elementviewvar.col(0);
}
template<int NDim> void Function<NDim>::toLinFormFromNodal() {
    auto elementviewvar = elementview();
    elementviewvar.rightCols(elementviewvar.cols()-1).colwise() -= elementviewvar.col(0);
}
template<int NDim> Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> Function<NDim>::toNodalFromLinForm() const {
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> elementviewvar = elementview();
    elementviewvar.rightCols(elementviewvar.cols()-1).colwise() += elementviewvar.col(0);
    return elementviewvar;
}
template<int NDim> void Function<NDim>::toModalFromLinForm(){
    elementview() = elementview() * fromLinFormToModalTransf;
};
template<int NDim> void Function<NDim>::toLinFormFromModal(){
    elementview() = elementview() * fromModalToLinFormTransf;
}

template<int NDim> auto Function<NDim>::elementVec(const MeshElementIterator<NDim> t_elementID) const {
    return vec.template segment<NDim+1>(t_elementID.currentElementID*(mesh->nVertElement));
}
template<int NDim> auto Function<NDim>::elementVec(const MeshSubsetElementIterator<NDim> t_elementID) const {
    return vec.template segment<NDim+1>(t_elementID.currentElementID*(mesh->nVertElement));
}
template<int NDim> auto Function<NDim>::elementVec(const int t_elementID) const{
    return vec.template segment<NDim+1>(t_elementID*(mesh->nVertElement));
}
template<int NDim> auto Function<NDim>::elementVec(const MeshElementIterator<NDim> t_elementID) {
    return vec.template segment<NDim+1>(t_elementID.currentElementID*(mesh->nVertElement));
}
template<int NDim> auto Function<NDim>::elementVec(const MeshSubsetElementIterator<NDim> t_elementID) {
    return vec.template segment<NDim+1>(t_elementID.currentElementID*(mesh->nVertElement));
}
template<int NDim> auto Function<NDim>::elementVec(const int t_elementID) {
    return vec.template segment<NDim+1>(t_elementID*(mesh->nVertElement));
}


template<> const inline GeometryCore::LinTransform<1> Function<1>::fromModalToLinFormTransf = [](){
    GeometryCore::LinTransform<1> tmp ;
    tmp << 1.,0. , std::sqrt(3.), -2.*std::sqrt(3.);
    return tmp;
}();
// ({{1.,0.},{std::sqrt(3.), -2.*std::sqrt(3.)}});
template<> const inline GeometryCore::LinTransform<1> Function<1>::fromLinFormToModalTransf (Function<1>::fromModalToLinFormTransf.inverse());

template<> const inline GeometryCore::LinTransform<2> Function<2>::fromModalToLinFormTransf = [](){
    GeometryCore::LinTransform<2> tmp ;
    tmp << std::sqrt(2.),0.,0. ,2., -6., 0.,2., 0., -6. ;
    return tmp;
}();
 //({{std::sqrt(2.),0.,0.},{2., -6., 0.},{2., 0., -6.}});
template<> const inline GeometryCore::LinTransform<2> Function<2>::fromLinFormToModalTransf (Function<2>::fromModalToLinFormTransf.inverse());

template<> const inline GeometryCore::LinTransform<3> Function<3>::fromModalToLinFormTransf 
= [](){
    GeometryCore::LinTransform<3> tmp ;
    tmp << std::sqrt(2.), 0., 0. ,0., std::sqrt(10.), -4.*std::sqrt(10.), 0., 0.,
   std::sqrt(10.), 0., -4.*std::sqrt(10.), 0., std::sqrt(10.), 0., 0., -4.*std::sqrt(10.);
    return tmp;
}();
 //({{std::sqrt(2.), 0., 0. ,0.},{std::sqrt(10.), -4.*std::sqrt(10.), 0., 0.},
 //   {std::sqrt(10.), 0., -4.*std::sqrt(10.), 0.},{std::sqrt(10.), 0., 0., -4.*std::sqrt(10.)}});
template<> const inline GeometryCore::LinTransform<3> Function<3>::fromLinFormToModalTransf (Function<3>::fromModalToLinFormTransf.inverse());


template<int NDim>
auto Function<NDim>::evaluateBoolCondOnSections(const std::function<bool(Real)>& boolFunctionToEvaluate, const RequirementType requirementType) const {
    Eigen::Matrix<bool,1,Eigen::Dynamic> toreturn;
    if (requirementType == any) { toreturn = Eigen::Matrix<bool,1,Eigen::Dynamic>::Constant(1, mesh->numberSections, false);}
    else { toreturn = Eigen::Matrix<bool,1,Eigen::Dynamic>::Constant(1, mesh->numberSections, true);}
    int i = 0;
    for(auto sec = sectionIterator(); sec.unfinished; ++sec){
        MeshSubset<NDim> subset(mesh, sec);
        for (auto elem = subset.elementIterator(); elem.unfinished; ++elem){
            for ( VertexInElemIterator<NDim> vert; vert.unfinished; ++vert) {
                if (requirementType == any)
                    toreturn(0,i) = toreturn(0,i) || boolFunctionToEvaluate(operator()(elem, vert));
                else
                    toreturn(0,i) = toreturn(0,i) && boolFunctionToEvaluate(operator()(elem, vert));
            }
        }
        ++i;
    }
    return toreturn;
}

} // namespace Tortoise

#endif /* Function_hpp */

