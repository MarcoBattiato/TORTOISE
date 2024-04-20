//
//  PolygonalFunction.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 4/11/22.
//

#ifndef PolygonalFunction_hpp
#define PolygonalFunction_hpp

#include <Geometry/Structured/DirectSpace/Polygonal/PolygonalMeshIterators.hpp>
#include <Generics/Features/VectorSpace.hpp>

namespace Tortoise {

template<int NDim> class PolygonalFunction :
public Features::MathFieldSpace<PolygonalFunction<NDim>> ,
public Features::AsymmetricMathFieldSpaceByValue<PolygonalFunction<NDim>,Real> ,
public Features::AsymmetricVectorSpace< PolygonalFunction<NDim>, FunctionElement<NDim>> ,
public Features::AsymmetricMathFieldSpace< PolygonalFunction<NDim>, std::function<Real(Point<NDim>)>>{

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Main properties
    //~~~~~~~~~~~~~~~~~~
    //~~~~~~~~~~~~~~~~~~

public:
    const PolygonalMesh<NDim>*   mesh;
    DataVector                   vec;

    static const GeometryCore::LinTransform<NDim>  fromLinFormToModalTransf;
    static const GeometryCore::LinTransform<NDim>  fromModalToLinFormTransf;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Main methods
    //~~~~~~~~~~~~~~~~~~
    //~~~~~~~~~~~~~~~~~~
public:
    //=======================================================
    // Constructors
    //===================
    explicit PolygonalFunction(const PolygonalMesh<NDim>& t_mesh);                                                    // Constructs a constant 0 function
    PolygonalFunction(const PolygonalMesh<NDim>& t_mesh, const Real t_value);                                // Constructs a constant function
    PolygonalFunction(const PolygonalMesh<NDim>& t_mesh, const DataVector& t_nodalValues);                            // Constructs using the values at nodes
    PolygonalFunction(const PolygonalMesh<NDim>& t_mesh, const std::function<Real(Point<NDim>)>& t_f);       // Constructs using the function passed
    PolygonalFunction(const PolygonalMesh<NDim>& t_mesh, const Function<NDim>& t_f);       // Constructs using the function passed
    explicit PolygonalFunction(const PolygonalMesh<NDim>* t_mesh);                                                    // Constructs a constant 0 function
    PolygonalFunction(const PolygonalMesh<NDim>* t_mesh, const Real t_value);                                // Constructs a constant function
    PolygonalFunction(const PolygonalMesh<NDim>* t_mesh, const DataVector& t_nodalValues);                            // Constructs using the values at nodes
    PolygonalFunction(const PolygonalMesh<NDim>* t_mesh, const std::function<Real(Point<NDim>)>& t_f);       // Constructs using the function passed
    PolygonalFunction(const PolygonalMesh<NDim>* t_mesh, const Function<NDim>& t_f);       // Constructs using the function passed

    PolygonalFunction(const std::function<Real(Real)>& t_f, const PolygonalFunction<NDim>& t_other);// Constructs a function of a function f(other)
    PolygonalFunction(const std::function<Real(std::vector<Real>)>& t_f, const std::vector<PolygonalFunction<NDim>>& t_others);// Constructs a function of a function f(other[0], other[1], ...)

    //********************************
    //* Arithmetic
    //********************************
    // ATTENTION!!!!!!! Multiplications and divisions by anything except a scalar do not conserve precision!
    // If you don't know what this means, do not use those operations. You will get unexpected results in many cases
    PolygonalFunction& operator=(const PolygonalFunction<NDim>& other);
    PolygonalFunction& operator=(PolygonalFunction<NDim>&& other);
    PolygonalFunction& operator+=(const PolygonalFunction<NDim>& other);
    PolygonalFunction& operator-=(const PolygonalFunction<NDim>& other);
    PolygonalFunction& operator*=(const PolygonalFunction<NDim>& other);                          // !!!! Loses precision. DO NOT USE if you do not know what this means!!!!
    PolygonalFunction& operator/=(const PolygonalFunction<NDim>& other);                          // !!!! Loses precision. DO NOT USE if you do not know what this means!!!!

    PolygonalFunction<NDim> operator-() const;
        
    PolygonalFunction& operator=(const Real scalar);
    PolygonalFunction& operator+=(const Real scalar);
    PolygonalFunction& operator-=(const Real scalar);
    PolygonalFunction& operator*=(const Real scalar);
    PolygonalFunction& operator/=(const Real scalar);
//    friend PolygonalFunction<NDim> operator/(const Real lhs, PolygonalFunction<NDim> rhs);
        
    PolygonalFunction& operator+=(const FunctionElement<NDim>& other);
    PolygonalFunction& operator-=(const FunctionElement<NDim>& other);
        
    PolygonalFunction& operator=(const std::function<Real(Point<NDim>)>& other);
    PolygonalFunction& operator+=(const std::function<Real(Point<NDim>)>& other);
    PolygonalFunction& operator-=(const std::function<Real(Point<NDim>)>& other);
    PolygonalFunction& operator*=(const std::function<Real(Point<NDim>)>& other);
    PolygonalFunction& operator/=(const std::function<Real(Point<NDim>)>& other);
//    friend PolygonalFunction<NDim> operator/(const std::function<Real(Point<NDim>)>& lhs, PolygonalFunction<NDim> rhs);

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
    PolygonalSectionIterator<NDim>       sectionIterator() const;
    PolygonalElementIterator<NDim>       elementIterator() const;
    VertexInElemIterator<NDim>           vertexInElemIterator() const;
        
    // Random accessors
    CartIndex<NDim>             randomSection() const;
    int                         randomElement() const;

    //=======================================================
    // Technicalities
    //=======================================================
public:
    auto elementview();         // Returns a map. Can be used to modify directly vec, but it is accessed in element view
    auto elementview() const ;  // Returns a map. Cannot be used to modify directly vec, since it is constant
    void toNodalFromLinForm() ;      // Modifies the vec (the user should remember to eventually modify it back into linear form representation)
    void toLinFormFromNodal() ;    // Modifies the vec
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> toNodalFromLinForm() const ; // When the object is const it cannot modify the vec, but instead it returns a matrix. It is up to the user to store it
    void toModalFromLinForm();
    void toLinFormFromModal();

    
    // Remove the possibility of passing temporaries to the constructor
    explicit PolygonalFunction(const PolygonalMesh<NDim>&& t_mesh) = delete;
    PolygonalFunction(const PolygonalMesh<NDim>&& t_mesh, const Real t_value) = delete;
    PolygonalFunction(const PolygonalMesh<NDim>&& t_mesh, const DataVector& t_nodalValues) = delete;
    PolygonalFunction(const PolygonalMesh<NDim>&& t_mesh, const std::function<Real(Point<NDim>)>& t_f) = delete;
    PolygonalFunction(const PolygonalMesh<NDim>&& t_mesh, const Function<NDim>& t_f)= delete;
    
};

Plotter3D& operator << (Plotter3D& plotter, const PolygonalFunction<2>& function);
Plotter2D& operator << (Plotter2D& plotter, const PolygonalFunction<1>& function);

template<int NDim> PolygonalFunction<NDim> operator/(const Real lhs, PolygonalFunction<NDim> rhs);
template<int NDim> PolygonalFunction<NDim> operator/(const std::function<Real(Point<NDim>)>& lhs, PolygonalFunction<NDim> rhs);

template<int NDim> void swap(PolygonalFunction<NDim>& first, PolygonalFunction<NDim>& second);                                // Swap functionality

template<int NDim> auto sin(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::sin(z);}, other); }
template<int NDim> auto cos(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::cos(z);}, other); }
template<int NDim> auto tan(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::tan(z);}, other); }
template<int NDim> auto asin(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::asin(z);}, other); }
template<int NDim> auto acos(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::acos(z);}, other); }
template<int NDim> auto atan(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::atan(z);}, other); }
template<int NDim> auto cosh(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::cosh(z);}, other); }
template<int NDim> auto sinh(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::sinh(z);}, other); }
template<int NDim> auto tanh(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::tanh(z);}, other); }
template<int NDim> auto pow(const PolygonalFunction<NDim>& other, Real exponent) {return PolygonalFunction<NDim>([exponent](Real z){return std::pow(z, exponent);}, other); }
template<int NDim> auto sqrt(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::sqrt(z);}, other); }
template<int NDim> auto log(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::log(z);}, other); }
template<int NDim> auto floor(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::floor(z);}, other); }
template<int NDim> auto ceil(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::ceil(z);}, other); }
template<int NDim> auto abs(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::fabs(z);}, other); }
template<int NDim> auto fabs(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::fabs(z);}, other); }
template<int NDim> auto exp(const PolygonalFunction<NDim>& other) { return PolygonalFunction<NDim>([](Real z){return std::exp(z);}, other); }

template<int NDim> auto min(const PolygonalFunction<NDim>& f0, const PolygonalFunction<NDim>& f1) { return PolygonalFunction<NDim>([](const std::vector<Real>& z){return std::min(z[0],z[1]);}, {f0, f1});}
template<int NDim> auto max(const PolygonalFunction<NDim>& f0, const PolygonalFunction<NDim>& f1) { return PolygonalFunction<NDim>([](const std::vector<Real>& z){return std::max(z[0],z[1]);}, {f0, f1});}
template<int NDim> auto max(const PolygonalFunction<NDim>& f0, double val) { return PolygonalFunction<NDim>([val](const Real z){return std::max(z,val);}, f0);}
template<int NDim> auto max(double val, const PolygonalFunction<NDim>& f0) { return PolygonalFunction<NDim>([val](const Real z){return std::max(z,val);}, f0);}
template<int NDim> auto min(const PolygonalFunction<NDim>& f0, double val) { return PolygonalFunction<NDim>([val](const Real z){return std::min(z,val);}, f0);}
template<int NDim> auto min(double val, const PolygonalFunction<NDim>& f0) { return PolygonalFunction<NDim>([val](const Real z){return std::min(z,val);}, f0);}

//=======================================================
// Implementation Details
//===================

template<int NDim> auto PolygonalFunction<NDim>::elementview() {
    return Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vec.data(),  mesh->numberElements, mesh->nVertElement);
}
template<int NDim> auto PolygonalFunction<NDim>::elementview() const {
    return Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> (vec.data(), mesh->numberElements ,mesh->nVertElement );
}
template<int NDim> void PolygonalFunction<NDim>::toNodalFromLinForm() {
    auto elementviewvar = elementview();
    elementviewvar.rightCols(elementviewvar.cols()-1).colwise() += elementviewvar.col(0);
}
template<int NDim> void PolygonalFunction<NDim>::toLinFormFromNodal() {
    auto elementviewvar = elementview();
    elementviewvar.rightCols(elementviewvar.cols()-1).colwise() -= elementviewvar.col(0);
}
template<int NDim> Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> PolygonalFunction<NDim>::toNodalFromLinForm() const {
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> elementviewvar = elementview();
    elementviewvar.rightCols(elementviewvar.cols()-1).colwise() += elementviewvar.col(0);
    return elementviewvar;
}
template<int NDim> void PolygonalFunction<NDim>::toModalFromLinForm(){
    elementview() = elementview() * fromLinFormToModalTransf;
};
template<int NDim> void PolygonalFunction<NDim>::toLinFormFromModal(){
    elementview() = elementview() * fromModalToLinFormTransf;
}


} // namespace Tortoise
#endif /* PolygonalFunction_hpp */
