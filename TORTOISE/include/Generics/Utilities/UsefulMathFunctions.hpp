//
//  UsefulMathFunctions.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 20/9/22.
//

#ifndef UsefulMathFunctions_hpp
#define UsefulMathFunctions_hpp

#include <Eigen/Dense>
#include <cassert>

namespace Tortoise {

namespace Utilities {

constexpr inline unsigned int factorial(unsigned int t_n){
    if (t_n <= 0)
        return 1;
    return t_n * factorial(t_n - 1);
};

// The following function is taken from https://stackoverflow.com/questions/1505675/power-of-an-integer-in-c
// The author is Matthieu M. https://stackoverflow.com/users/147192/matthieu-m
template <unsigned int p> int constexpr intPower(const int x) {
  if constexpr (p == 0) return 1;
  if constexpr (p == 1) return x;

  int tmp = intPower<p / 2>(x);
  if constexpr ((p % 2) == 0) { return tmp * tmp; }
  else { return x * tmp * tmp; }
}

//int intPower(const int x, const int p) {
//  if (p == 0) return 1;
//  if (p == 1) return x;
//
//  int tmp = intPower(x,p/2);
//  if ((p % 2) == 0) { return tmp * tmp; }
//  else { return x * tmp * tmp; }
//}



// Calculates the binomial coefficient
constexpr int binomial(const int n, const int k){
    assert(n>=0 &&k>=0);
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    int result = n - k + 1;
    for (int i = 1; i < k; ++i) { result = result * (n - k + 1 + i) / (i + 1);}
    return result;
}


// Constructs an upper triangular matrix using a vector containing the non-zero elements
// diagonalShift =  0 (upper diagonal starting at the diagonal)   = 1 (upper diagonal starting one row above the diagonal)  etc
template <int NDim, typename Derived> auto upperTriangular(const Eigen::MatrixBase<Derived>& coeff, const int nUpperRowsSkipped){
    const int numCoeff = (NDim-nUpperRowsSkipped)*(NDim-nUpperRowsSkipped+1)/2;
    assert(nUpperRowsSkipped < NDim && coeff.size() == numCoeff);
    Eigen::Matrix<typename Derived::Scalar, NDim, NDim> toReturn(Eigen::Matrix<typename Derived::Scalar, NDim, NDim>::Zero());
    for (int row = nUpperRowsSkipped; row < NDim; ++row) {
        toReturn.block(row, row, 1, NDim-row) = coeff.reshaped(1, numCoeff).block(0, (row-nUpperRowsSkipped)*(2*NDim-nUpperRowsSkipped+1-row)/2, 1, NDim-row);
    }
    return toReturn;
}

} // namespace Utilities

} // namespace Tortoise

#endif /* UsefulMathFunctions_hpp */
