//
//  LinearForm.hpp
//  TORTOISE_Development
//
//  Created by Marco Battiato on 2/11/22.
//

#ifndef LinearForm_hpp
#define LinearForm_hpp

#include <Geometry/GeometryCore/Point.hpp>
#include <Geometry/GeometryCore/ReferenceElement.hpp>

namespace Tortoise {

namespace GeometryCore {

//*******************************
// Non-homogeneous linear forms  ( Functionals of Points   f(P) )
//*******************************
// ===  Types
// ================
template <int NDim> using LinearForm                        = Eigen::Matrix<Real, 1, NDim+1>;
// It is defines such that
// LinearForm<2> exmplLinForm = {-0.3,2.3,1.5};
// applying exmplLinForm to exmplPoint returns   exmplLinForm(0) + exmplLinForm(1)+exmplPoint(0) + exmplLinForm(2)+exmplPoint(1)
template <int NDim> using VectorLinearForm                  = Eigen::Matrix<Real, Eigen::Dynamic, NDim+1>;
// In VectorLinearForm the dimension related to the number of linear forms is left to be Dynamic since there is a hard limit to the size of statically stored Matrices
template <int NDim, int NForms> using ArrayLinearForm       = Eigen::Matrix<Real, NForms, NDim+1>;
// A vector of linear forms with statically assigned dimension (to be used for small number of linear forms)

// ===  Creators
// ================
template <int NDim> LinearForm<NDim> createLinearForm(const std::array<Real,NDim+1> nodevalue){
    LinearForm<NDim> toreturn;
    toreturn(0) = nodevalue[0];
    for (int i=1; i<NDim+1; ++i){ toreturn(i) = nodevalue[i]-nodevalue[0];}
    return toreturn;
};
// Creates a linear form that takes the given values at the reference tetrahedron's node coordinate. The nodes are ordered {0,0,0}, {1,0,0}, {0,1,0}. {0,0,1}
// nodevalue is in the form Eigen::Matrix<Real, nLinForms, NDim+1>
template <typename Derived> auto createLinearForm(const Eigen::MatrixBase<Derived>& nodevalue){
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> toreturn = nodevalue;
    toreturn.rightCols(toreturn.cols()-1).colwise() -= toreturn.col(0);
    return toreturn;
};
// Creates a linear form that takes the given values at the reference tetrahedron's node coordinate. The nodes are ordered {0,0,0}, {1,0,0}, {0,1,0}. {0,0,1}
template <int NDim> LinearForm<NDim> createLinearForm(const std::function<Real(Point<NDim>)>& t_f){
    LinearForm<NDim> values;
    for (int i = 0; i<NDim+1; ++i){ values(i) = t_f(refNodes<NDim>.col(i));}
    return createLinearForm(values);
}
// Creates a linear form that approximates the given function. The linear form has the same values of the function at the reference tetrahedron's node coordinate.

// ===  Operations
// ================
// This function conveniently works for all the following combinations
//  ->  apply(LinForm, Point)           ->  apply(LinForm, PointVector)             ->  apply(LinForm, ArrayPoint)
//  ->  apply(VectorLinForm, Point)     ->  apply(VectorLinForm, VectorPoint)       ->  apply(VectorLinForm, ArrayPoint)
//  ->  apply(ArrayLinearForm, Point)   ->  apply(ArrayLinearForm, VectorPoint)     ->  apply(ArrayLinearForm, ArrayPoint)
template <typename DerivedA,typename DerivedB> auto apply(const Eigen::MatrixBase<DerivedA>& linearForm, const Eigen::MatrixBase<DerivedB>& point){
    return (linearForm.rightCols(linearForm.cols()-1)*point).colwise()+linearForm.col(0);
}

//*******************************
// Linear operators on the non-homogeneous linear forms
//*******************************
template <int NDim> using LinTransform   = Eigen::Matrix<Real, NDim+1, NDim+1>;

} // namespace GeometryCore 

} // namespace Tortoise

#endif /* LinearForm_hpp */
