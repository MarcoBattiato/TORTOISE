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
//  DormandPrince54.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 7/9/21.
//
// It calculates the time propagation with adaptive time step Dormand Prince 54
// It needs to work on a container that describes a function of time into a vector space
// The easiest way to ensure that the container has all the necessary methods is to inherit UnivariateRelationContainer
// and ensuring that containerArg in UnivariateRelationContainer is a container of Real type values like for instance double
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
// -> timeContainer: DP54 will start the propagation from the last entry. If the timeContainer is empty if will stop. The calculated dense time steps will be appendend
//                      to this container. This makes this parameter both input and output
// -> finalTime:     DP54 will propagate until the time reaches finalTime. Notice that this is not the time on top of the time of the last entry in timeContainer. If
//                      the time of the last entry in the timeContainer is larger than finalTime, DP853 will do nothing
// -> denseTimeStep: time step of the dense output (see description of the numerical method)
// -> timeOp:        operator that gives the time derivative.
//                      It must be > timeOp(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> ContainerType::ValType
// -> sqrErrorForm:  give the square of the normalised global error given a local error and the y status.
//                      It must be > sqrErrorForm(ContainerType::ValType localErrY, ContainerType::ValType currentY) -> ContainerType::ArgType
// -> outputFunct:   used to produce user defined output. It will attached to the standard output of DP853. It is suggested to keep the user-defined output short
//                      and within a single line for readibility. Also it is better not to add newlines, again to improve readability.
//                      It must be in the form > outputFunct(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> void
//
// The function returns a timeContainer with all the calculated adaptive steps. Usually these are not important, but are returned for testing reasons
//
// =====================================
// Example
// =====================================
//
//#include "DiscretizationCore.hpp"
//#include "MiscUtilities.hpp"
//#include "PhysicsCore.hpp"
//#include <iostream>
//
//using namespace Tortoise;
//using std::cout;
//
//Real timeOp(const Real time, const Real stat){return cos(time);}
//Real errorForm(const Real err, const Real stat) { return 10000000000000.*err*err;}
//void outputFunct(const Real time, const Real stat){ }
//
//int main(int argc, const char * argv[]) {
//
//    UnivariateRelationContainer<std::vector<Real>,std::vector<Real>>   solution(0.0,0.0);
//    dormandPrince54(solution, 10.0, 1., timeOp, errorForm, outputFunct);
//    cout << "\n";
//    for (int i=0; i<solution.size(); ++i){
//        cout << solution.arg(i) << " " << solution[i] << " [" << sin(solution.arg(i))-solution[i]<< "] \n";
//    }
//
//
//    return 0;
//}

#ifndef DormandPrince54_hpp
#define DormandPrince54_hpp

#include <stdio.h>
#include <iomanip>

namespace Tortoise {

template <typename ContainerType, typename timeOpFunctType, typename errorFunctType, typename outputFunctType> void dormandPrince54
 (ContainerType&                                timeContainer,      // input/output
  const typename ContainerType::ArgType         finalTime,          // end time of propagation
  const typename ContainerType::ArgType         denseTimeStep,      // step for dense grid output
  timeOpFunctType                               timeOp,             // operator that gives the derivative
  // It must be >   timeOp(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> ContainerType::ValType
  errorFunctType                                normErrorForm,       //  normalised error
  // It must be >    normErrorForm(ContainerType::ValType localErrY, ContainerType::ValType currentY) -> ContainerType::ArgType
  outputFunctType                               outputFunct         // Defines some output operations
  // It must be in the form >     outputFunct(ContainerType::ArgType currentTime, ContainerType::ValType currentY)
) {
  
    // Check if the container contains at least an element to be used as an initial contition
    assert(!timeContainer.empty());                          // If the container is empty there is no usable starting point, it quits

    typedef typename ContainerType::ArgType  RealType;       // Short name for the Real type (for instance double)
    typedef typename ContainerType::ValType  vecType;        // Short name for the vector type (for instance
    
    // ==================
    // Butcher tableau
    const RealType c2 = 1./5., c3 = 3./10., c4 = 4./5., c5 = 8./9. /*, c6 = 1., c7 = 1. */ ;
    const RealType  a21 = 1./5.,
                    a31 = 3./40.,       a32 = 9./40.,
                    a41 = 44./45.,      a42 = -56./15.,         a43 = 32./9.,
                    a51 = 19372./6561., a52 = -25360./2187.,    a53 = 64448./6561., a54 = -212./729.,
                    a61 = 9017./3168.,  a62 = -355./33.,        a63 = 46732./5247., a64 = 49./176.,     a65 = -5103./18656.;
                    // a71 = b1, a72 = b2, a73 = b3, a74 = b4, a75 = b5, a76 = b6;     // Notice the FSAL property
    const RealType  b1 = 35./384.,     /* b2 = 0., */           b3 = 500./1113.,    b4 = 125./192.,     b5 = -2187./6784.,     b6 = 11./84.; // b7 = 0.;
    const RealType  bh1 = 5179./57600., /* bh2 = 0. */          bh3 = 7571./16695., bh4 = 393./640.,    bh5 = -92097./339200., bh6 = 187./2100., bh7=1./40.;

    // ==================
    // Allocation of the Y on adaptive time
    vecType currentY(timeContainer.back());
    // Allocation of the RK stages, attempted solution, and error
    // They are initialised by copy, but their initial value is irrelevant.
    // The initialisation by copy is done in order to allow for non-default constructible types,
    // since even non-default constructible types, will most likely have a copy constructor.
    // In principle one could ask the user for a user-defined constructor,
    // however the gain in efficiency is minimal and not worth the increase in technicality
    vecType k1(currentY), k2(k1), k3(k1), k4(k1), k5(k1), k6(k1), k7(k1);
    vecType y5th(k1);   // Will contain the solution at the 5th order
    vecType err(k1);    // Will contain the solution at the 5th order minus the solution at the 4th order

    
    // ==================
    // Time stepping data
    RealType currentTime       = timeContainer.argback();        // The time propagation starts with the last time already present in the container
    RealType timeStepSize      = denseTimeStep;                  // As initial guess for the time step uses the dense step
    RealType denseTime         = currentTime + denseTimeStep;    // Initialises the dense time
    const RealType minStepSize = denseTimeStep / 10000.;
    bool     notfailed         = true;                           // Becomes false if the step size becomes too small
    bool     rejected          = false;                          // Becomes true if the previous step was rejected
    
    // ==================
    // Adaptive step parameters
    const RealType minscale = 0.333333;                     // Max reduction of step size
    const RealType maxscale = 3.0;                          // Max increment of step size
    const RealType safetyFactor = 0.98;                     // Safety factor in step size increment
    const RealType alpha    = 0.7/5.;                       // PI stepsize control parameter
    const RealType beta     = 0.4/5.;                       // PI stepsize control parameter
    RealType olderror       = 1.0;
    
    // ==================
    // DP54 Start
    // Since the method has the FSAL (First Same As Last) property, k1 will be set to be equal to k7, except for the 1st iteration. We set it now
    k1 = timeOp(currentTime, timeContainer.back());
    
    while (currentTime < finalTime && notfailed){
        
        // ==================
        // Console Output
        std::cout << "** Time: " << std::fixed << currentTime << " timestep: " << timeStepSize;

        
        // ==================
        // Calculation of DP54 stages, new function and error function
        k2 = timeOp(currentTime + c2 * timeStepSize,
                    currentY + (a21 * timeStepSize) * k1);
        k3 = timeOp(currentTime + c3 * timeStepSize,
                    currentY + (a31 * timeStepSize) * k1 + (a32 * timeStepSize) * k2 );
        k4 = timeOp(currentTime + c4 * timeStepSize,
                    currentY + (a41 * timeStepSize) * k1 + (a42 * timeStepSize) * k2 + (a43 * timeStepSize) * k3 );
        k5 = timeOp(currentTime + c5 * timeStepSize,
                    currentY + (a51 * timeStepSize) * k1 + (a52 * timeStepSize) * k2 + (a53 * timeStepSize) * k3 + (a54 * timeStepSize) * k4 );
        k6 = timeOp(currentTime + timeStepSize, // we do not multiply by b6 as it is 1.
                    currentY + (a61 * timeStepSize) * k1 + (a62 * timeStepSize) * k2 + (a63 * timeStepSize) * k3 + (a64 * timeStepSize) * k4 + (a65 * timeStepSize) * k5 );
        
        y5th = currentY + (timeStepSize * b1) * k1 + (timeStepSize * b3) * k3 + (timeStepSize * b4) * k4 + (timeStepSize * b5) * k5 + (timeStepSize * b6) * k6;  // the step size has been multiplied by each coefficient, since the multiplication of the solution vector by a scalar is a much more expensive operation than the multiplication between two scalars
        
        k7 = timeOp(currentTime + timeStepSize, y5th);
        
        err = (timeStepSize * (b1 - bh1)) * k1 + (timeStepSize * (b3 - bh3)) * k3 + (timeStepSize * (b4 - bh4)) * k4 + (timeStepSize * (b5 - bh5)) * k5 + (timeStepSize * (b6 - bh6)) * k6 - (timeStepSize * bh7) * k7;
        
        // ==================
        // Evaluation normalised error
        RealType errNorm = normErrorForm(err, y5th);             // Normalised error calculated as indicated by the user      abs(y-yh)/errtol
        std::cout << " error: " << errNorm ;

        // Checking error
        if ( errNorm <= 1. ) {  // Case when error is acceptable
            
            // ==================
            // Console Output
            std::cout << " ACCEPTED <<< ";
            outputFunct(currentTime + timeStepSize, y5th);        // user-defined output
            
            
            // ==================
            // Dense Output
            while(denseTime <= currentTime + (1.0001) * timeStepSize   && denseTime <= finalTime + .0001 * timeStepSize ){
                
                // Taken from "Solving Ordinary Differential Equations" by Hairer, Norset, Wanner
//                RealType theta = ( denseTime - currentTime ) / timeStepSize;
//                RealType  btilde1 = theta * theta * (3. - 2. * theta) * b1 + theta * (theta - 1.) * (theta - 1.) - theta * theta * (theta - 1.) * (theta - 1.) * 5. * (2558722523. - 31403016. * theta) / 11282082432.;
//                RealType  btilde3 = theta * theta * (3. - 2. * theta) * b3 + theta * theta * (theta - 1.) * (theta - 1.) * 100. *  (882725551. - 15701508. * theta) / 32700410799.;
//                RealType  btilde4 = theta * theta * (3. - 2. * theta) * b4 - theta * theta * (theta - 1.) * (theta - 1.) * 25. * (443332067. - 31403016 * theta) / 1880347072.;
//                RealType  btilde5 = theta * theta * (3. - 2. * theta) * b5 + theta * theta * (theta - 1.) * (theta - 1.) * 32805. * (23143187. - 3489224. * theta) / 199316789632.;
//                RealType  btilde6 = theta * theta * (3. - 2. * theta) * b6 - theta * theta * (theta - 1.) * (theta - 1.) * 55. * (29972135. - 7076736 * theta) / 822651844.;
//                RealType  btilde7 = theta * theta * (theta - 1.) + theta * theta * (theta - 1.) * (theta - 1.) * 10. * (7414447. - 829305. * theta) / 29380423.;
//
//                timeContainer.emplace_back(denseTime,
//                                           currentY + (timeStepSize * btilde1) * k1 + (timeStepSize * btilde3) * k3 + (timeStepSize * btilde4) * k4 + (timeStepSize * btilde5) * k5 + (timeStepSize * btilde6) * k6 + (timeStepSize * btilde7) * k7);
                
                // Dense output at the 4th order taken from J. R. Dormand "Numerical Methodds for Differential Equations"
                RealType s = (denseTime - currentTime) / timeStepSize;
                RealType btilde1 = - (435.*s*s*s - 1184.*s*s+1098.*s-384.)/(384.);
                RealType btilde3 = 500.*s*(6.*s*s -14.*s + 9.)/(1113.);
                RealType btilde4 = - 125.*s*(9.*s*s -16.*s + 6.)/(192.);
                RealType btilde5 = 729.*s*(35.*s*s -64.*s + 26.)/(6784.);
                RealType btilde6 = - 11.*s*(3.*s - 2.)*(5.*s - 6.)/(84.);
                RealType btilde7 = s*(s - 1.)*(5.*s - 3.)/(2.);
                timeContainer.emplace_back(denseTime, currentY + (s * timeStepSize * btilde1) * k1 + (s * timeStepSize * btilde3) * k3 + (s * timeStepSize * btilde4) * k4 + (s * timeStepSize * btilde5) * k5 + (s * timeStepSize * btilde6) * k6 + (s * timeStepSize * btilde7) * k7);
                
                denseTime += denseTimeStep;
            }
            std::cout << "\n";
            // ==================
            // Iteration concluding tasks (FSAL)
            currentY    = y5th;
            currentTime += timeStepSize;
            k1          = k7;
            
            // ==================
            // Time step update
            
            // The step size will then be modified by multiplying it by a scale factor
            RealType       scale;
            if(errNorm == 0.){  // If the error is 0, the scale is set to the maximum one
                scale = maxscale;
            } else {
                scale = safetyFactor * std::pow(1./errNorm, alpha) * std::pow(olderror, beta);
                scale = std::min(std::max(scale, minscale), maxscale);                // The scale is however capped between min and max allowed values
            }
            
            if (rejected){
                timeStepSize *= std::min(1.0, scale);
            } else {
                timeStepSize *= scale;
            }
            rejected = false;    // To rember that the step was accepted
            
            olderror = errNorm;
                
        } else {
            RealType scale = safetyFactor * std::pow(1./errNorm, 1/5.);
            scale = std::min(std::max(scale, minscale), maxscale);
            timeStepSize *= scale;

            rejected = true;    // To rember that the step was rejected

            // ==================
            // Console Output
            std::cout << " REJECTED !!! ";
            outputFunct(currentTime + timeStepSize, y5th);        // user-defined output
            std::cout << "\n";
            if(timeStepSize < minStepSize){
                std::cout << "The step size has become too small to be handled" << std::endl;
                notfailed = false;
            }
        }
    }
}


template <typename ContainerType, typename timeOpFunctType, typename enforceValueFunctType, typename errorFunctType, typename outputFunctType> void dormandPrince54
 (ContainerType&                                timeContainer,      // input/output
  const typename ContainerType::ArgType         finalTime,          // end time of propagation
  const typename ContainerType::ArgType         denseTimeStep,      // step for dense grid output
  timeOpFunctType                               timeOp,             // operator that gives the derivative
  // It must be >   timeOp(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> ContainerType::ValType
  enforceValueFunctType                         enforceValue,       // Function that modifies the value of the solution to enforce certain constraints
  // it must be >   enforceValue(ContainerType::ArgType currentTime, ContainerType::ValType& currentY) -> void
  errorFunctType                                normErrorForm,       //  normalised error
  // It must be >    normErrorForm(ContainerType::ValType localErrY, ContainerType::ValType currentY) -> ContainerType::ArgType
  outputFunctType                               outputFunct         // Defines some output operations
  // It must be in the form >     outputFunct(ContainerType::ArgType currentTime, ContainerType::ValType currentY)
) {
  
    // Check if the container contains at least an element to be used as an initial contition
    assert(!timeContainer.empty());                          // If the container is empty there is no usable starting point, it quits

    typedef typename ContainerType::ArgType  RealType;       // Short name for the Real type (for instance double)
    typedef typename ContainerType::ValType  vecType;        // Short name for the vector type (for instance
    
    // ==================
    // Butcher tableau
    const RealType c2 = 1./5., c3 = 3./10., c4 = 4./5., c5 = 8./9. /*, c6 = 1., c7 = 1. */ ;
    const RealType  a21 = 1./5.,
                    a31 = 3./40.,       a32 = 9./40.,
                    a41 = 44./45.,      a42 = -56./15.,         a43 = 32./9.,
                    a51 = 19372./6561., a52 = -25360./2187.,    a53 = 64448./6561., a54 = -212./729.,
                    a61 = 9017./3168.,  a62 = -355./33.,        a63 = 46732./5247., a64 = 49./176.,     a65 = -5103./18656.;
                    // a71 = b1, a72 = b2, a73 = b3, a74 = b4, a75 = b5, a76 = b6;     // Notice the FSAL property
    const RealType  b1 = 35./384.,     /* b2 = 0., */           b3 = 500./1113.,    b4 = 125./192.,     b5 = -2187./6784.,     b6 = 11./84.; // b7 = 0.;
    const RealType  bh1 = 5179./57600., /* bh2 = 0. */          bh3 = 7571./16695., bh4 = 393./640.,    bh5 = -92097./339200., bh6 = 187./2100., bh7=1./40.;

    // ==================
    // Allocation of the Y on adaptive time
    vecType currentY(timeContainer.back());
    // Allocation of the RK stages, attempted solution, and error
    // They are initialised by copy, but their initial value is irrelevant.
    // The initialisation by copy is done in order to allow for non-default constructible types,
    // since even non-default constructible types, will most likely have a copy constructor.
    // In principle one could ask the user for a user-defined constructor,
    // however the gain in efficiency is minimal and not worth the increase in technicality
    vecType k1(currentY), k2(k1), k3(k1), k4(k1), k5(k1), k6(k1), k7(k1);
    vecType y5th(k1);   // Will contain the solution at the 5th order
    vecType err(k1);    // Will contain the solution at the 5th order minus the solution at the 4th order

    
    // ==================
    // Time stepping data
    RealType currentTime       = timeContainer.argback();        // The time propagation starts with the last time already present in the container
    RealType timeStepSize      = denseTimeStep;                  // As initial guess for the time step uses the dense step
    RealType denseTime         = currentTime + denseTimeStep;    // Initialises the dense time
    const RealType minStepSize = denseTimeStep / 10000.;
    bool     notfailed         = true;                           // Becomes false if the step size becomes too small
    bool     rejected          = false;                          // Becomes true if the previous step was rejected
    
    // ==================
    // Adaptive step parameters
    const RealType minscale = 0.333333;                     // Max reduction of step size
    const RealType maxscale = 3.0;                          // Max increment of step size
    const RealType safetyFactor = 0.98;                     // Safety factor in step size increment
    const RealType alpha    = 0.7/5.;                       // PI stepsize control parameter
    const RealType beta     = 0.4/5.;                       // PI stepsize control parameter
    RealType olderror       = 1.0;
    
    // ==================
    // DP54 Start
    // Since the method has the FSAL (First Same As Last) property, k1 will be set to be equal to k7, except for the 1st iteration. We set it now
    k1 = timeOp(currentTime, enforceValue(currentTime, timeContainer.back()));
    
    int count = 0;
    
    while (currentTime < finalTime && notfailed){
        ++count;
        // ==================
        // Console Output
        std::cout << "** Time: " << std::fixed << currentTime << " timestep: " << timeStepSize;

        
        // ==================
        // Calculation of DP54 stages, new function and error function
        k2 = timeOp(currentTime + c2 * timeStepSize,
                    enforceValue(currentTime + c2 * timeStepSize, currentY + (a21 * timeStepSize) * k1));
        k3 = timeOp(currentTime + c3 * timeStepSize,
                    enforceValue(currentTime + c3 * timeStepSize, currentY + (a31 * timeStepSize) * k1 + (a32 * timeStepSize) * k2 ));
        k4 = timeOp(currentTime + c4 * timeStepSize,
                    enforceValue(currentTime + c4 * timeStepSize, currentY + (a41 * timeStepSize) * k1 + (a42 * timeStepSize) * k2 + (a43 * timeStepSize) * k3 ));
        k5 = timeOp(currentTime + c5 * timeStepSize,
                    enforceValue(currentTime + c5 * timeStepSize, currentY + (a51 * timeStepSize) * k1 + (a52 * timeStepSize) * k2 + (a53 * timeStepSize) * k3 + (a54 * timeStepSize) * k4 ));
        k6 = timeOp(currentTime + timeStepSize, // we do not multiply by b6 as it is 1.
                    enforceValue(currentTime + timeStepSize, currentY + (a61 * timeStepSize) * k1 + (a62 * timeStepSize) * k2 + (a63 * timeStepSize) * k3 + (a64 * timeStepSize) * k4 + (a65 * timeStepSize) * k5 ));
        
        y5th = enforceValue(currentTime + timeStepSize, currentY + (timeStepSize * b1) * k1 + (timeStepSize * b3) * k3 + (timeStepSize * b4) * k4 + (timeStepSize * b5) * k5 + (timeStepSize * b6) * k6);  // the step size has been multiplied by each coefficient, since the multiplication of the solution vector by a scalar is a much more expensive operation than the multiplication between two scalars
        
        k7 = timeOp(currentTime + timeStepSize, y5th);
        
        err = y5th - enforceValue(currentTime + timeStepSize, currentY + (timeStepSize * bh1) * k1 + (timeStepSize * bh3) * k3 + (timeStepSize * bh4) * k4 + (timeStepSize * bh5) * k5 + (timeStepSize * bh6) * k6 + (timeStepSize * bh7) * k7);
        
        // ==================
        // Evaluation normalised error
        RealType errNorm = normErrorForm(err, y5th);             // Normalised error calculated as indicated by the user      abs(y-yh)/errtol
//        std::cout << " error: " << errNorm ;

        // Checking error
        if ( errNorm <= 1. ) {  // Case when error is acceptable
            
            // ==================
            // Console Output
//            std::cout << " ACCEPTED <<< ";
            outputFunct(currentTime + timeStepSize, y5th);        // user-defined output
            
            
            // ==================
            // Dense Output
            while(denseTime <= currentTime + (1.0001) * timeStepSize   && denseTime <= finalTime + .0001 * timeStepSize ){
                
                // Dense output at the 4th order taken from J. R. Dormand "Numerical Methodds for Differential Equations"
                RealType s = (denseTime - currentTime) / timeStepSize;
                RealType btilde1 = - (435.*s*s*s - 1184.*s*s+1098.*s-384.)/(384.);
                RealType btilde3 = 500.*s*(6.*s*s -14.*s + 9.)/(1113.);
                RealType btilde4 = - 125.*s*(9.*s*s -16.*s + 6.)/(192.);
                RealType btilde5 = 729.*s*(35.*s*s -64.*s + 26.)/(6784.);
                RealType btilde6 = - 11.*s*(3.*s - 2.)*(5.*s - 6.)/(84.);
                RealType btilde7 = s*(s - 1.)*(5.*s - 3.)/(2.);
                timeContainer.emplace_back(denseTime, enforceValue(denseTime, currentY + (s * timeStepSize * btilde1) * k1 + (s * timeStepSize * btilde3) * k3 + (s * timeStepSize * btilde4) * k4 + (s * timeStepSize * btilde5) * k5 + (s * timeStepSize * btilde6) * k6 + (s * timeStepSize * btilde7) * k7));
                
                denseTime += denseTimeStep;
            }
            std::cout << "\n";
            // ==================
            // Iteration concluding tasks (FSAL)
            currentY    = y5th;
            currentTime += timeStepSize;
            k1          = k7;
            
            // ==================
            // Time step update
            
            // The step size will then be modified by multiplying it by a scale factor
            RealType       scale;
            if(errNorm == 0.){  // If the error is 0, the scale is set to the maximum one
                scale = maxscale;
            } else {
                scale = safetyFactor * std::pow(1./errNorm, alpha) * std::pow(olderror, beta);
                scale = std::min(std::max(scale, minscale), maxscale);                // The scale is however capped between min and max allowed values
            }
            
            if (rejected){
                timeStepSize *= std::min(1.0, scale);
            } else {
                timeStepSize *= scale;
            }
            rejected = false;    // To rember that the step was accepted
            
            olderror = errNorm;
                
        } else {
            RealType scale = safetyFactor * std::pow(1./errNorm, 1/5.);
            scale = std::min(std::max(scale, minscale), maxscale);
            timeStepSize *= scale;

            rejected = true;    // To rember that the step was rejected

            // ==================
            // Console Output
            std::cout << " REJECTED !!! ";
            outputFunct(currentTime + timeStepSize, y5th);        // user-defined output
            std::cout << "\n";
            if(timeStepSize < minStepSize){
                std::cout << "The step size has become too small to be handled" << std::endl;
                notfailed = false;
            }
        }
    }
    std::cout << "Number of adaptive steps: " << count << "\n";
}



} // namespace Tortoise



#endif /* DormandPrince54_hpp */
