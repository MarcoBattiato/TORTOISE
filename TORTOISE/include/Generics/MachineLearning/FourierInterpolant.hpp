//
//  FourierInterpolant.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 30/9/23.
//

#ifndef FourierInterpolant_hpp
#define FourierInterpolant_hpp

#include "Generics/MachineLearning/SparseFunction.hpp"
#include "Generics/MachineLearning/BasisFunctionsSets.hpp"

namespace Tortoise {
namespace SparseFunctionals {

template <int order, int nVar, typename ScalarType = double>
using FourierInterpolant = SparseFunction<order, nVar, FourierBasis, ScalarType>;

} // namespace SparseFunction
} // namespace Tortoise


#endif /* FourierInterpolant_hpp */
