//
//  SparseFunctionStruct.h
//
//  Created by Marco Battiato on 26/8/23.
//

#ifndef SparseFunctionStruct_hpp
#define SparseFunctionStruct_hpp

#include <vector>
#include <Eigen/Dense>
#include <cassert>
#include <iostream>
//#include <ranges>  // used for std::views::drop


namespace Tortoise {
namespace SparseFunctionals {
namespace Utilities {
template <typename ScalarType> inline void printOrderVectors(const std::vector<int>& orders, ScalarType coeff){
    std::cout << "[" << orders[0];
    // Unfortunately std::views doesn't seem to work on MacOS  https://stackoverflow.com/questions/73628848/is-the-stdviews-namespace-not-available-in-xcodes-c?rq=1 the code below is kept for future reference
    //        for (auto ord : orders | std::views::drop(1)){}
    for (size_t i=1; i != orders.size(); ++i){
        std::cout << "," << orders[i];
    }
    std::cout << "](" << coeff << ")\n";
}

} // namespace Utilities


namespace Internal {

// ====== sparseFunctionStruct =======
template <int totalOrder, int maxOrder, int nVar, template<int> typename basisType, typename ScalarType>
class sparseFunctionStruct;

// ======
template <template<int> typename basisType, typename ScalarType>
struct sparseFunctionStruct<0, 0, 1, basisType, ScalarType>{
    //===================
    ScalarType coeff{};     // c_[0,...] where ... is calling order list
    //===================
    //~~~~ Constants
    constexpr static int       numBasisFunction(int n) {
        return n+1;
    }
    //~~~~ Construction/Assignment
    ScalarType * const firstCoeff() {
        return &coeff;
    }
    //~~~~ Evaluation
    template <typename Derived>
    auto                operator()(const Eigen::ArrayBase<Derived>& x) const {
        return basisType<0>::basisF(x.col(0)) * coeff;
    };
    template <typename Derived1, typename Derived2>
    static void         applyAllBasisFunctions(const Eigen::ArrayBase<Derived1>& x, Eigen::ArrayBase<Derived2>& result, int& counter) {
        result.col(counter++) = basisType<0>::basisF(x.col(0)) ;
    }
    template <typename Derived1, typename Derived2, typename Derived3>
    static void         applyAllBasisFunctions(const Eigen::ArrayBase<Derived1>& x, Eigen::ArrayBase<Derived2>& result, int& counter, const Eigen::ArrayBase<Derived3>& restOfExpression) {
        result.col(counter++) = basisType<0>::basisF(x.col(0)) * restOfExpression;
    }
    //~~~~ I/O
    void                print() const {
        std::cout << "f0[x1]";
    }
    std::vector<int>&   printCoeff(std::vector<int>& orders) const {
        orders[0] = 0;
        Utilities::printOrderVectors(orders, coeff);
        return orders;
    }
    
};

// ======
template <int totalOrder, int maxOrder, template<int> typename basisType, typename ScalarType>
struct sparseFunctionStruct<totalOrder, maxOrder, 1, basisType, ScalarType>{
    using RestType = sparseFunctionStruct<totalOrder-1, maxOrder-1, 1, basisType, ScalarType>;
    //===================
    ScalarType coeff{};      // c_[maxOrder, ...]
    RestType rest;
    //===================
    //~~~~ Constants
    constexpr static int       numBasisFunction(int n) {
        return RestType::numBasisFunction(n)+1;
    }
    //~~~~ Construction/Assignment
    ScalarType * const firstCoeff() const {
        return &coeff;
    }
    //~~~~ Evaluation
    template <typename Derived>
    auto                operator()(const Eigen::ArrayBase<Derived>& x) const {
        return basisType<maxOrder>::basisF(x.col(0)) * coeff + rest(x);
    };
    template <typename Derived1, typename Derived2>
    static void         applyAllBasisFunctions(const Eigen::ArrayBase<Derived1>& x, Eigen::ArrayBase<Derived2>& result, int& counter) {
        result.col(counter++) = basisType<maxOrder>::basisF(x.col(0));
        RestType::applyAllBasisFunctions(x, result, counter);
    }
    template <typename Derived1, typename Derived2, typename Derived3>
    static void         applyAllBasisFunctions(const Eigen::ArrayBase<Derived1>& x, Eigen::ArrayBase<Derived2>& result, int& counter, const Eigen::ArrayBase<Derived3>& restOfExpression) {
        result.col(counter++) = basisType<maxOrder>::basisF(x.col(0)) * restOfExpression;
        RestType::applyAllBasisFunctions(x, result, counter, restOfExpression);
    }
    //~~~~ I/O
    void                print() const {
        std::cout << "(f" << maxOrder << "[x1]";
        std::cout << "+";
        rest.print();
        std::cout << ")";
    }
    std::vector<int>&   printCoeff(std::vector<int>& orders) const {
        orders[0] = totalOrder;
        Utilities::printOrderVectors(orders, coeff);
        rest.printCoeff(orders);
        return orders;
    }
};

// ======
template <int totalOrder, int nVar, template<int> typename basisType, typename ScalarType>
struct sparseFunctionStruct<totalOrder, 0, nVar, basisType, ScalarType>{
    using TopType = sparseFunctionStruct<totalOrder, totalOrder , nVar-1, basisType, ScalarType>;
    //===================
    TopType top;
    //===================
    //~~~~ Constants
    constexpr static int       numBasisFunction(int n) {
        return TopType::numBasisFunction(n);
    }
    //~~~~ Construction/Assignment
    ScalarType * const firstCoeff() {
        return top.firstCoeff();
    }
    
    //~~~~ Evaluation
    template <typename Derived>
    auto                operator()(const Eigen::ArrayBase<Derived>& x) const {
        return basisType<0>::basisF(x.col(nVar-1)) * top(x);
    };
    template <typename Derived1, typename Derived2>
    static void         applyAllBasisFunctions(const Eigen::ArrayBase<Derived1>& x, Eigen::ArrayBase<Derived2>& result, int& counter) {
        TopType::applyAllBasisFunctions(x, result, counter, basisType<0>::basisF(x.col(nVar-1)));
    }
    template <typename Derived1, typename Derived2, typename Derived3>
    static void         applyAllBasisFunctions(const Eigen::ArrayBase<Derived1>& x, Eigen::ArrayBase<Derived2>& result, int& counter, const Eigen::ArrayBase<Derived3>& restOfExpression) {
        TopType::applyAllBasisFunctions(x, result, counter, basisType<0>::basisF(x.col(nVar-1)) * restOfExpression);
    }
    //~~~~ I/O
    void                print() const {
        std::cout << "(f0[x" << nVar << "]";
        top.print();
        std::cout << ")";
    }
    std::vector<int>&   printCoeff(std::vector<int>& orders) const {
        orders[nVar-1] = 0;
        top.printCoeff(orders);
        return orders;
    }
    
};

// ======
template <int totalOrder, int maxOrder, int nVar, template<int> typename basisType, typename ScalarType>
struct sparseFunctionStruct {
    using TopType = sparseFunctionStruct<totalOrder - maxOrder, totalOrder - maxOrder, nVar-1, basisType, ScalarType>;
    using RestType = sparseFunctionStruct<totalOrder, maxOrder-1, nVar, basisType, ScalarType>;
    //===================
    TopType top;
    RestType rest;
    //===================
    //~~~~ Constants
    constexpr static int       numBasisFunction(int n) {
        return TopType::numBasisFunction(RestType::numBasisFunction(n));
    }
    //~~~~ Construction/Assignment
    ScalarType * const firstCoeff() {
        return top.firstCoeff();
    }
    
    //~~~~ Evaluation
    template <typename Derived>
    auto                operator()(const Eigen::ArrayBase<Derived>& x) const {
        return basisType<maxOrder>::basisF(x.col(nVar-1)) * top(x) + rest(x);
    };
    template <typename Derived1, typename Derived2>
    static void         applyAllBasisFunctions(const Eigen::ArrayBase<Derived1>& x, Eigen::ArrayBase<Derived2>& result, int& counter) {
        TopType::applyAllBasisFunctions(x, result, counter, basisType<maxOrder>::basisF(x.col(nVar-1)));
        RestType::applyAllBasisFunctions(x, result, counter);
    }
    template <typename Derived1, typename Derived2, typename Derived3>
    static void         applyAllBasisFunctions(const Eigen::ArrayBase<Derived1>& x, Eigen::ArrayBase<Derived2>& result, int& counter, const Eigen::ArrayBase<Derived3>& restOfExpression) {
        TopType::applyAllBasisFunctions(x, result, counter, basisType<maxOrder>::basisF(x.col(nVar-1)) * restOfExpression);
        RestType::applyAllBasisFunctions(x, result, counter, restOfExpression);
    }
    //~~~~ I/O
    void                print() const {
        std::cout << "(f" << maxOrder << "[x" << nVar << "]*";
        top.print();
        std::cout << "+";
        rest.print();
        std::cout << ")";
    }
    std::vector<int>&   printCoeff(std::vector<int>& orders) const {
        orders[nVar-1] = maxOrder;
        top.printCoeff(orders);
        rest.printCoeff(orders);
        return orders;
    }
};

} // namespace Internal
} // namespace SparseFunctionals
} // namespace Tortoise

#endif /* SparseFunctionStruct_hpp */
