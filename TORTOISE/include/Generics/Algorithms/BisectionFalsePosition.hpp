//
//  BisectionFalsePosition.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 13/12/22.
//

#ifndef BisectionFalsePosition_hpp
#define BisectionFalsePosition_hpp

#include <tuple>
#include <cassert>
#include <cmath>

namespace Tortoise {

namespace Algorithms {

template <typename timeOpFunctType, typename xType, typename yType>
xType bisectionFalsePosition(timeOpFunctType f,                 // Function
                             xType left,                        // Starting point
                             xType right,                       // Step sizes on each vaariable for first iteration
                             int maxIterations,                 // Max number of iterations
                             xType positionTolerance,           // accepted error in the position of the root
                             yType valueTolerance,              // accepted distance from 0 of the value of the function at the numerical root
                             int& nIterations){                 // RETURN VALUE: number of iterations performed
    
    positionTolerance = std::abs(positionTolerance); // makes sure that no negative numbers are passed
    valueTolerance = std::abs(valueTolerance);       // makes sure that no negative numbers are passed
    
    yType fLeft = f(left);
    yType fRight = f(right);
    
    assert( fRight * fRight >= 0. );
    
    auto secant = (left * fRight - right * fLeft)/(fRight - fLeft);
    auto fSecant = f(secant);
    auto middle = (left + right)/2.;
    auto fMiddle = f(middle);
    
    auto error = 10. * positionTolerance + 1.;
    nIterations = 0;
    
    while ( error > positionTolerance && std::abs(fSecant) > valueTolerance && nIterations < maxIterations){
        ++nIterations;
        
        // We want to make sure that the secant is between middle and secondEdge
        if ( (secant - middle)*(right - left) < 0. ){
            auto temp = left; left = right; right = temp;
            auto fTemp = fLeft; fLeft = fRight; fRight = fTemp;
        }
        
        if (fRight * fSecant <= 0.){            // then the zero is between secant and secondEdge
            left = secant;
            fLeft = fSecant;
        } else if (fMiddle * fSecant <= 0.) {   // then the zero is between middle and secant
            left = middle;
            fLeft = fMiddle;
            right = secant;
            fRight = fSecant;
        } else {                                 // then the zero is between firstEdge and middle
            right = middle;
            fRight = fMiddle;
        }
        
        auto oldSolution = secant;
        secant = (left * fRight - right * fLeft)/(fRight - fLeft);
        fSecant = f(secant);
        middle = (left + right)/2.;
        fMiddle = f(middle);
        error = std::abs(secant - oldSolution);
    }
    
    return secant;
    
}

} // namespace Algorithms

} // namespace Tortoise

#endif /* BisectionFalsePosition_hpp */
