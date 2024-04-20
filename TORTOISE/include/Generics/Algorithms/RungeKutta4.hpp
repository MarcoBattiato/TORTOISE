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
//  RungeKutta4.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 23/6/21.
//
// It calculates the time propagation with Runge Kutta 4
// It needs to work on a container that describes a function of time into a vector space
// The easiest way to ensure that the container has all the necessary methods is to inherit UnivariateRelationContainer
// and ensuring that containerArg in UnivariateRelationContainer is a container of real type values like, for instance, double
// and that containerVal is a container of vectors on the containerArg.
// The easiest way to ensure that is that the type containerVal::value_type inherits from VectorSpace
//
// IMPORTANT: remeber to inehrit everything as >> public << otherwise its members might be considered
// as private (if you are defining a class and not a struct)
//
// If the user does not want to inherit from the suggested base classes, the user must ensure that the following type/methods are defined
//  ContainerType must define public:
//      typename    ->  ContainerType::ArgType
//      typename    ->  ContainerType::ValType
//      constructor ->  ContainerType(const ContainerType::ArgType (&), const ContainerType::ValType (&));
//      accessor    ->  ContainerType::ValType     back()
//      accessor    ->  ContainerType::ArgType     timesback();
//      add element ->  (void) emplace_back(const ContainerType::ArgType (&), const ContainerType::ValType (&));
//
//   ContainerType::ArgType must be a continuous scalar type like double and must have a constructor from double
//
//   ContainerType::ValType must define
//      copy constr ->  ContainerType::ValType(const ContainerType::ValType (&));
//      +, - between ContainerType::ValType
//      * by the scalar ContainerType::ArgType
//
//
// =====================================
// Input description
// -> timeContainer:    RK4 will start the propagation from the last entry. If the timeContainer is empty if will stop. The calculated dense time steps will be appendend
//                          to this container. This makes this parameter both input and output
// -> finalTime:        RK4 will propagate until the time reaches finalTime. Notice that this is not the propagation time (on top of the time of the last entry in timeContainer). If
//                          the time of the last entry in the timeContainer is larger than finalTime, RK4 will do nothing.
// -> timeOp:           operator that gives the time derivative.
//                          It must be > timeOp(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> ContainerType::ValType
// -> toExecuteFunct:   this function is executed on the solution at every step. It is used to enforce conditions on the solution that go beyond the differenntial equation
//                          It must be > toExecuteFunct(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> ContainerType::ValType
// -> outputFunct:      used to produce user defined output. It will attached to the standard output of DP853. It is suggested to keep the user-defined output short
//                          and within a single line for readibility. Also it is better not to add newlines, again to improve readability.
//                          It must be in the form > outputFunct(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> void


#ifndef RungeKutta4_hpp
#define RungeKutta4_hpp

#include <stdio.h>

namespace Tortoise {

namespace Algorithms {

template <typename ContainerType, typename timeOpFunctType> void rungeKutta4
 (ContainerType&                                timeContainer,      // input/output
  const typename ContainerType::ArgType         finalTime,          // end time of propagation: notice that this is the final time and not the propagation time. If the last time in the container is larger that this, no time propagation will be done.
  const typename ContainerType::ArgType         timeStepSize,       // time step size
  timeOpFunctType                               timeOp              // operator that gives the derivative
//  It must be > timeOp(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> ContainerType::ValType
) {
    typedef typename ContainerType::ArgType  RealType;       // Short name for the real type (for instance double)
    typedef typename ContainerType::ValType  vecType;        // Short name for the vector type (for instance

    vecType k0(timeContainer.back()),
            k1(timeContainer.back()),
            k2(timeContainer.back()),
            k3(timeContainer.back());                        // Will hold the RK4 intermediate steps
    RealType currentTime = timeContainer.argback();          // The time propagation starts with the last time already present in the container
    
    while ( timeContainer.argback() + timeStepSize < finalTime ){
        k0 = timeStepSize * timeOp(currentTime, timeContainer.back());
        k1 = timeStepSize * timeOp(currentTime + 0.5 * timeStepSize, timeContainer.back() + 0.5 * k0);
        k2 = timeStepSize * timeOp(currentTime + 0.5 * timeStepSize, timeContainer.back() + 0.5 * k1);
        k3 = timeStepSize * timeOp(currentTime + timeStepSize, timeContainer.back() + k2);
        currentTime += timeStepSize;
        timeContainer.emplace_back(currentTime, timeContainer.back() + (1./6.) * k0 + (1./3.) * k1 + (1./3.) * k2 + (1./6.) * k3);
    }
    
    RealType finalStepSize = finalTime - timeContainer.argback();
    k0 = finalStepSize * timeOp(currentTime, timeContainer.back());
    k1 = finalStepSize * timeOp(currentTime + 0.5 * finalStepSize, timeContainer.back() + 0.5 * k0);
    k2 = finalStepSize * timeOp(currentTime + 0.5 * finalStepSize, timeContainer.back() + 0.5 * k1);
    k3 = finalStepSize * timeOp(currentTime + finalStepSize, timeContainer.back() + k2);
    currentTime += finalStepSize;
    timeContainer.emplace_back(currentTime, timeContainer.back() + (1./6.) * k0 + (1./3.) * k1 + (1./3.) * k2 + (1./6.) * k3);
}

template <typename ContainerType, typename timeOpFunctType, typename outputFunctType> void rungeKutta4
 (ContainerType&                                timeContainer,      // input/output
  const typename ContainerType::ArgType         finalTime,          // end time of propagation: notice that this is the final time and not the propagation time. If the last time in the container is larger that this, no time propagation will be done.
  const typename ContainerType::ArgType         timeStepSize,       // time step size
  timeOpFunctType                               timeOp,             // operator that gives the derivative
//  It must be > timeOp(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> ContainerType::ValType
  outputFunctType                               outputFunct         // Defines some output operations
//  It must be in the form > outputFunct(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> void
) {
    typedef typename ContainerType::ArgType  RealType;       // Short name for the real type (for instance double)
    typedef typename ContainerType::ValType  vecType;        // Short name for the vector type (for instance

    vecType k0(timeContainer.back()),
            k1(timeContainer.back()),
            k2(timeContainer.back()),
            k3(timeContainer.back());                        // Will hold the RK4 intermediate steps
    RealType currentTime = timeContainer.argback();          // The time propagation starts with the last time already present in the container
    
    while ( timeContainer.argback() + timeStepSize < finalTime ){
        k0 = timeStepSize * timeOp(currentTime, timeContainer.back());
        k1 = timeStepSize * timeOp(currentTime + 0.5 * timeStepSize, timeContainer.back() + 0.5 * k0);
        k2 = timeStepSize * timeOp(currentTime + 0.5 * timeStepSize, timeContainer.back() + 0.5 * k1);
        k3 = timeStepSize * timeOp(currentTime + timeStepSize, timeContainer.back() + k2);
        currentTime += timeStepSize;
        timeContainer.emplace_back(currentTime, timeContainer.back() + (1./6.) * k0 + (1./3.) * k1 + (1./3.) * k2 + (1./6.) * k3);
        outputFunct(currentTime, timeContainer.back());
    }
    
    RealType finalStepSize = finalTime - timeContainer.argback();
    k0 = finalStepSize * timeOp(currentTime, timeContainer.back());
    k1 = finalStepSize * timeOp(currentTime + 0.5 * finalStepSize, timeContainer.back() + 0.5 * k0);
    k2 = finalStepSize * timeOp(currentTime + 0.5 * finalStepSize, timeContainer.back() + 0.5 * k1);
    k3 = finalStepSize * timeOp(currentTime + finalStepSize, timeContainer.back() + k2);
    currentTime += finalStepSize;
    timeContainer.emplace_back(currentTime, timeContainer.back() + (1./6.) * k0 + (1./3.) * k1 + (1./3.) * k2 + (1./6.) * k3);
    outputFunct(currentTime, timeContainer.back());
}


template <typename ContainerType, typename timeOpFunctType, typename toExecuteFunctType, typename outputFunctType> void rungeKutta4
 (ContainerType&                                timeContainer,      // input/output
  const typename ContainerType::ArgType         finalTime,          // end time of propagation: notice that this is the final time and not the propagation time. If the last time in the container is larger that this, no time propagation will be done.
  const typename ContainerType::ArgType         timeStepSize,       // time step size
  timeOpFunctType                               timeOp,             // operator that gives the derivative
  toExecuteFunctType                            toExecuteFunct,     // function that will be executed at every time step (used to enforce certain properties to solution, beyond RK)
  outputFunctType                               outputFunct         // Defines some output operations
) {
    typedef typename ContainerType::ArgType  RealType;       // Short name for the real type (for instance double)
    typedef typename ContainerType::ValType  vecType;        // Short name for the vector type (for instance

    vecType k0(timeContainer.back()),
            k1(timeContainer.back()),
            k2(timeContainer.back()),
            k3(timeContainer.back());                        // Will hold the RK4 intermediate steps
    RealType currentTime = timeContainer.argback();          // The time propagation starts with the last time already present in the container
    
    while ( timeContainer.argback() + timeStepSize < finalTime ){
        k0 = timeStepSize * timeOp(currentTime, toExecuteFunct(currentTime, timeContainer.back()));
        k1 = timeStepSize * timeOp(currentTime + 0.5 * timeStepSize, toExecuteFunct(currentTime + 0.5 * timeStepSize, timeContainer.back() + 0.5 * k0));
        k2 = timeStepSize * timeOp(currentTime + 0.5 * timeStepSize, toExecuteFunct(currentTime + 0.5 * timeStepSize, timeContainer.back() + 0.5 * k1));
        k3 = timeStepSize * timeOp(currentTime + timeStepSize, toExecuteFunct(currentTime + timeStepSize, timeContainer.back() + k2));
        currentTime += timeStepSize;
        timeContainer.emplace_back(currentTime, toExecuteFunct(currentTime, timeContainer.back() + (1./6.) * k0 + (1./3.) * k1 + (1./3.) * k2 + (1./6.) * k3));
        outputFunct(currentTime, timeContainer.back());
    }
    RealType finalStepSize = finalTime - timeContainer.argback();
    k0 = finalStepSize * timeOp(currentTime, toExecuteFunct(currentTime, timeContainer.back()));
    k1 = finalStepSize * timeOp(currentTime + 0.5 * finalStepSize, toExecuteFunct(currentTime + 0.5 * timeStepSize, timeContainer.back() + 0.5 * k0));
    k2 = finalStepSize * timeOp(currentTime + 0.5 * finalStepSize, toExecuteFunct(currentTime + 0.5 * timeStepSize, timeContainer.back() + 0.5 * k1));
    k3 = finalStepSize * timeOp(currentTime + finalStepSize, toExecuteFunct(currentTime + timeStepSize, timeContainer.back() + k2));
    currentTime += finalStepSize;
    timeContainer.emplace_back(currentTime, toExecuteFunct(currentTime, timeContainer.back() + (1./6.) * k0 + (1./3.) * k1 + (1./3.) * k2 + (1./6.) * k3));
    outputFunct(currentTime, timeContainer.back());
}


} // namespace Algorithms 

} // namespace Tortoise
#endif /* RungeKutta4_hpp */
