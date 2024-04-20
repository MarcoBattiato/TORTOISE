//
//  Polynomial.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 30/9/23.
//

#ifndef Polynomial_hpp
#define Polynomial_hpp

#include "Generics/MachineLearning/SparseFunction.hpp"
#include "Generics/MachineLearning/BasisFunctionsSets.hpp"

namespace Tortoise {
namespace SparseFunctionals {

template <int order, int nVar, typename ScalarType = double>
struct Polynomial: public SparseFunction<order, nVar, TaylorBasis, ScalarType> {

    using SparseFunction<order, nVar, TaylorBasis, ScalarType>::coeffPtr;
    
public:
    //~~~~ Substituition
    template<typename Derived1, typename Derived2>
    Polynomial<order, 1> substitute (const Eigen::ArrayBase<Derived1>& vectorA, const Eigen::ArrayBase<Derived2>& vectorB) {
        assert(vectorA.rows() == nVar && vectorA.cols() == 1); assert(vectorB.rows() == nVar && vectorB.cols() == 1);
        
        //~~~ Necessary variables to be constructed
        Eigen::Array<ScalarType, order + 1, 1> wValues = Eigen::Array<ScalarType, order + 1, 1>::LinSpaced(order + 1, 0.0, 1.0);
        Polynomial<order, 1, ScalarType> smallPoly;
        
        //~~~ Calculations
        auto bigXValue = (vectorA.matrix() * wValues.matrix().transpose()).array().colwise() + vectorB;
        smallPoly.fitData(wValues, (*this)(bigXValue.transpose()));
        
        return smallPoly;
    }
    
    auto roots() const requires (nVar==1) {
        Eigen::EigenSolver<Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>> solver(buildCompanionMatrix());
        return solver.eigenvalues();
    }
    auto realRoots() const requires (nVar==1) {
        
        auto roots = (*this).roots();
        const int size = roots.size();
        int nReal = 0;
        for (int i = 0; i < size; ++i){
            if(roots(i).imag() == 0) {
                if (nReal != i) roots(nReal) = roots(i);
                ++nReal;
            }
        }
        roots.conservativeResize(nReal);
    }
    
private:
    //~~~~ Companion Matrix
    auto buildCompanionMatrix() const requires (nVar==1) {
        
        int skipped = 0;
        while ( std::fabs(*(coeffPtr+skipped)) < 1.e-30 ) {     //TODO: check
            ++skipped;
        }
        size_t degree = order-skipped;
        
        Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> companionMatrix(degree, degree);
        companionMatrix.setZero();
        // Fill the subdiagonal with ones
        for (size_t i = 1; i != degree; ++i) {
            companionMatrix(i, i - 1) = 1.0;
        }
        // Set the first row of the companion matrix using SparseFunction::Polynomial coefficient
        for (size_t i{}; i != degree; ++i) {
            companionMatrix(i , degree - 1) = - (coeffPtr[order- i]) / (coeffPtr[skipped]);
        }
        
        return companionMatrix;
    }
    
};


} // namespace SparseFunction
} // namespace Tortoise

#endif /* Polynomial_hpp */
