//
//  Levenberg-Marquardt.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 3/11/22.
//
// INCOMPLETE

#ifndef Levenberg_Marquardt_hpp
#define Levenberg_Marquardt_hpp

#include <stdio.h>
#include <array>
#include <cassert>
#include <Eigen/Dense>

namespace Tortoise {

namespace Algorithms {

template<typename FittingFunctType,
         typename GradientFittingFunctType,
         typename TargetDataType,
         typename SourceDataType,
         typename ParameterType>
auto LevenbergMarquardt(const FittingFunctType& f,
                        const GradientFittingFunctType& df,
                        const std::array<TargetDataType>& targetArray,
                        const std::array<SourceDataType>& sourceArray,
                        const Eigen::EigenBase<Derived>& initGuess){  // should be Eigen::Matrix<"", >
    assert(targetArray.size() >= 1 && "Dataset is empty!");
    assert(targetArray.size() == sourceArray.size() && "Data is inconsistent!");
    
    // Calculation error
    auto deviation = targetArray[0] - FittingFunctType(sourceArray[0], guess);
    auto error = (deviation * deviation).evaluate();
    for (int i=1; i<sourceArray.size(); ++i){
        auto deviation = targetArray[i] - FittingFunctType(sourceArray[i], guess);
        error += deviation * deviation;
    }
    
    auto grad    = GradientFittingFunctType(sourceArray[0], guess).evaluate();
    auto hessian = (grad.transpose() * grad).evaluate();
    
    
    
    
}

} // namespace Algorithms 

} // namespace Tortoise

#endif /* Levenberg_Marquardt_hpp */
