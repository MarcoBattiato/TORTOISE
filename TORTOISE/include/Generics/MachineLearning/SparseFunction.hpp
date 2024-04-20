//
//  SparseFunction.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 30/9/23.
//

#ifndef SparseFunction_hpp
#define SparseFunction_hpp

#include "Generics/MachineLearning/SparseFunctionStruct.hpp"

#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>

namespace Tortoise {
namespace SparseFunctionals {

template <int order, int nVar, template<int> typename basisType, typename ScalarType>
struct SparseFunction: public Internal::sparseFunctionStruct<order, order, nVar, basisType, ScalarType> {
    
    ScalarType * const coeffPtr;   // Pointer to the first coefficient
    
    //~~~~~ Constant & type definition
    using ParentType = Internal::sparseFunctionStruct<order, order, nVar, basisType, ScalarType>;
    constexpr static int numBasisFunction(){
        int n = 0;
        return ParentType::numBasisFunction(n);
    }
    SparseFunction(): coeffPtr(ParentType::firstCoeff()) {}
    
    //~~~ Construction/Assignment
    template<typename Derived1, typename Derived2>
    ScalarType fitData(const Eigen::ArrayBase<Derived1>& x, const Eigen::ArrayBase<Derived2>& y){
        assert(x.cols() == nVar); assert(x.rows() == y.rows()); assert(y.cols() == 1);
        
        Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic> psyx = applyAllBasisFunctions(x).matrix();
        Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic> psyxTy = psyx.transpose()*y.matrix();
        Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic> psyTPsy = psyx.transpose()*psyx;
        
        Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic> solution = psyTPsy.llt().solve(psyxTy);
        
        for (size_t i{}; i != numBasisFunction(); ++i) {
            coeffPtr[i] = solution(i);
        }
        
        return ((operator()(x)-y)/(y.maxCoeff() - y.minCoeff())).matrix().squaredNorm()/y.size();
    }
    template <typename Derived>
    auto operator()(const Eigen::ArrayBase<Derived>& x) const {
        assert(x.cols() == nVar);
        return ParentType::operator()(x);
    }
    template <typename Derived>
    void assignCoefficients(const Eigen::DenseBase<Derived>& RowVector) {
        assert(RowVector.cols() == 1);
        for (size_t i{}; i != order; ++i) {
            *(coeffPtr+i) = RowVector(i);
        }
    }
    template <typename Derived>
    static auto applyAllBasisFunctions(const Eigen::ArrayBase<Derived>& x) {
        assert(x.cols() == nVar);
        using returnType = Eigen::Array<typename Eigen::ArrayBase<Derived>::Scalar, Eigen::ArrayBase<Derived>::RowsAtCompileTime, Eigen::Dynamic>;
        returnType result;
        result.resize(x.rows(), numBasisFunction());
        int currentBasisFunction = 0;
        ParentType::applyAllBasisFunctions(x, result, currentBasisFunction);
        return result;
    }
    
    //~~~~ I/O
    void printCoeff() const {
        std::vector<int> orders;
        orders.resize(nVar);
        ParentType::printCoeff(orders);
    }
    void saveToFile(const std::string& fileName) const {
        std::ofstream file(fileName, std::ios::out | std::ios::binary | std::ios::trunc);
        file.write((char*) coeffPtr, numBasisFunction()*sizeof(ScalarType) );
        file.close();
    }
    void loadFromFile(const std::string& fileName) const {
        std::ifstream file(fileName, std::ios::in | std::ios::binary);
        file.read((char*) coeffPtr, numBasisFunction()*sizeof(ScalarType) );
        file.close();
    }
    
};

} // namespace SparseFunction
} // namespace Tortoise

#endif /* SparseFunction_hpp */
