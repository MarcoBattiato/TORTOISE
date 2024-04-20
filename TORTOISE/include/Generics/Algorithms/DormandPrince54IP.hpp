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
//  DormandPrince54IP.hpp
//  Tests
//
//  Created by Marco Battiato on 9/9/21.
//
//



#ifndef DormandPrince54IP_hpp
#define DormandPrince54IP_hpp

#include <stdio.h>


namespace Tortoise {

namespace Algorithms {

template <typename ContainerType, typename timeOpFunctType, typename timePropFunctType, typename errorFunctType, typename outputFunctType> void dormandPrince54IP
(ContainerType&                                timeContainer,      // input/output
 const typename ContainerType::ArgType         finalTime,          // end time of propagation
 const typename ContainerType::ArgType         denseTimeStep,      // step for dense grid output
 timeOpFunctType                               timeOpSlow,         // operator that gives the derivative due to the slow dynamics
 // It must be >   timeOpSlow(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> ContainerType::ValType
 timePropFunctType                             propagatorFast,     // >> propagator << of the fast dynamics
 // It must be >   propagatorFast(ContainerType::ArgType startTime, ContainerType::ArgType endTime, ContainerType::ValType currentY) -> ContainerType::ValType
 errorFunctType                                sqrErrorForm,       // squared normalised error
 // It must be >    normErrorForm(ContainerType::ValType localErrY, ContainerType::ValType currentY) -> ContainerType::ArgType
 outputFunctType                               outputFunct         // Defines some output operations
// It must be in the form >     outputFunct(ContainerType::ArgType currentTime, ContainerType::ValType currentY)
)  {
    
    // Check if the container contains at least an element to be used as an initial contition
    if(timeContainer.empty()) { return ;}                             // If the container is empty there is no usable starting point, it quits
    
    typedef typename ContainerType::ArgType  RealType;       // Short name for the real type (for instance double)
    typedef typename ContainerType::ValType  vecType;        // Short name for the vector type (for instance
    
    // ==================
    // Butcher tableau
    const RealType c2 = 1./5., c3 = 3./10., c4 = 4./5., c5 = 8./9., c6 = 1., c7 = 1.;
    const RealType  a21 = 1./5.,
                    a31 = 3./40.,       a32 = 9./40.,
                    a41 = 44./45.,      a42 = -56./15.,         a43 = 32./9.,
                    a51 = 19372./6561., a52 = -25360./2187.,    a53 = 64448./6561., a54 = -212./729.,
                    a61 = 9017./3168.,  a62 = -355./33.,        a63 = 46732./5247., a64 = 49./176.,     a65 = -5103./18656.;
                    // a71 = b1, a72 = b2, a73 = b3, a74 = b4, a75 = b5, a76 = b6;     // Notice the FSAL property
    const RealType  b1 = 35./384.,     /* b2 = 0., */           b3 = 500./1113.,    b4 = 125./192.,     b5 = -2187./6784.,     b6 = 11./84.; // b7 = 0.;
    const RealType  bh1 = 5179./57600., /* bh2 = 0. */          bh3 = 7571./16695., bh4 = 393./640.,    bh5 = -92097./339200., bh6 = 187./2100., bh7=1./40.;

    // ==================
    // Time stepping data
    RealType currentTime       = timeContainer.argback();        // The time propagation starts with the last time already present in the container
    RealType timeStepSize      = denseTimeStep;                  // As initial guess for the time step uses the dense step
    RealType denseTime         = currentTime + denseTimeStep;    // Initialises the dense time
    const RealType minStepSize = denseTimeStep / 10000.;
    bool     notfailed         = true;                           // Becomes false if the step size becomes too small
    bool     rejected          = false;                          // Becomes true if the previous step was rejected

    
    // ==================
    // Allocation of the Y on adaptive time
    vecType currentY(timeContainer.back());
    vecType phi0(propagatorFast(currentTime, currentTime + timeStepSize, currentY));
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
    k1 = propagatorFast(currentTime, currentTime + timeStepSize, timeOpSlow(currentTime, timeContainer.back()));
    
    while (currentTime < finalTime && notfailed){
        
        // ==================
        // Console Output
        std::cout << "** Time: " << std::fixed << currentTime << " timestep: " << timeStepSize;
        
        
        // ==================
        // Calculation of DP54 stages, new function and error function
        k2 = propagatorFast(currentTime + c2 * timeStepSize, currentTime + timeStepSize,
                            timeOpSlow(currentTime + c2 * timeStepSize,
                                   propagatorFast(currentTime + timeStepSize, currentTime + c2 * timeStepSize,
                                                  phi0 + (a21 * timeStepSize) * k1
                                                  )
                                   )
                            );
        k3 = propagatorFast(currentTime + c3 * timeStepSize, currentTime + timeStepSize,
                            timeOpSlow(currentTime + c3 * timeStepSize,
                                   propagatorFast(currentTime + timeStepSize, currentTime + c3 * timeStepSize,
                                                  phi0 + (a31 * timeStepSize) * k1 + (a32 * timeStepSize) * k2
                                                  )
                                   )
                            );
        k4 = propagatorFast(currentTime + c4 * timeStepSize, currentTime + timeStepSize,
                            timeOpSlow(currentTime + c4 * timeStepSize,
                                   propagatorFast(currentTime + timeStepSize, currentTime + c4 * timeStepSize,
                                                  phi0 + (a41 * timeStepSize) * k1 + (a42 * timeStepSize) * k2 + (a43 * timeStepSize) * k3
                                                  )
                                   )
                            );
        k5 = propagatorFast(currentTime + c5 * timeStepSize, currentTime + timeStepSize,
                            timeOpSlow(currentTime + c5 * timeStepSize,
                                   propagatorFast(currentTime + timeStepSize, currentTime + c5 * timeStepSize,
                                                  phi0 + (a51 * timeStepSize) * k1 + (a52 * timeStepSize) * k2 + (a53 * timeStepSize) * k3 + (a54 * timeStepSize) * k4
                                                  )
                                   )
                            );
        k6 = timeOpSlow(currentTime + timeStepSize, // we do not multiply by b6 as it is 1.
                    phi0 + (a61 * timeStepSize) * k1 + (a62 * timeStepSize) * k2 + (a63 * timeStepSize) * k3 + (a64 * timeStepSize) * k4 + (a65 * timeStepSize) * k5 );
        
        y5th = phi0 + (timeStepSize * b1) * k1 + (timeStepSize * b3) * k3 + (timeStepSize * b4) * k4 + (timeStepSize * b5) * k5 + (timeStepSize * b6) * k6;  // the step size has been multiplied by each coefficient, since the multiplication of the solution vector by a scalar is a much more expensive operation than the multiplication between two scalars
        
        k7 = timeOpSlow(currentTime + timeStepSize, y5th);
        
        err = (timeStepSize * (b1 - bh1)) * k1 + (timeStepSize * (b3 - bh3)) * k3 + (timeStepSize * (b4 - bh4)) * k4 + (timeStepSize * (b5 - bh5)) * k5 + (timeStepSize * (b6 - bh6)) * k6 - (timeStepSize * bh7) * k7;
        
        // ==================
        // Evaluation normalised error
        RealType errNorm = sqrErrorForm(err, y5th);             // Normalised error squared calculated as indicated by the user
        std::cout << " error: " << errNorm ;
        
        // Checking error
        if ( errNorm <= 1. ) {  // Case when error is acceptable
            
            // ==================
            // Console Output
            std::cout << " ACCEPTED <<< ";
            outputFunct(currentTime + timeStepSize, y5th);        // user-defined output
            std::cout << "\n";
            
            // ==================
            // Dense Output
            while(denseTime <= currentTime + (1.0001) * timeStepSize   && denseTime <= finalTime + .0001 * timeStepSize ){
                
                RealType theta = ( denseTime - currentTime ) / timeStepSize;
                
                // Taken from "Solving Ordinary Differential Equations" by Hairer, Norset, Wanner
                RealType  btilde1 = theta * theta * (3. - 2. * theta) * b1 + theta * (theta - 1.) * (theta - 1.) - theta * theta * (theta - 1.) * (theta - 1.) * 5. * (2558722523. - 31403016. * theta) / 11282082432.;
                RealType  btilde3 = theta * theta * (3. - 2. * theta) * b3 + theta * theta * (theta - 1.) * (theta - 1.) * 100. *  (882725551. - 15701508. * theta) / 32700410799.;
                RealType  btilde4 = theta * theta * (3. - 2. * theta) * b4 - theta * theta * (theta - 1.) * (theta - 1.) * 25. * (443332067. - 31403016 * theta) / 1880347072.;
                RealType  btilde5 = theta * theta * (3. - 2. * theta) * b5 + theta * theta * (theta - 1.) * (theta - 1.) * 32805. * (23143187. - 3489224. * theta) / 199316789632.;
                RealType  btilde6 = theta * theta * (3. - 2. * theta) * b6 - theta * theta * (theta - 1.) * (theta - 1.) * 55. * (29972135. - 7076736 * theta) / 822651844.;
                RealType  btilde7 = theta * theta * (theta - 1.) + theta * theta * (theta - 1.) * (theta - 1.) * 10. * (7414447. - 829305. * theta) / 29380423.;
                
                timeContainer.emplace_back(denseTime,
                                           propagatorFast(currentTime + timeStepSize, denseTime,
                                                          phi0 + (timeStepSize * btilde1) * k1 + (timeStepSize * btilde3) * k3 + (timeStepSize * btilde4) * k4 + (timeStepSize * btilde5) * k5 + (timeStepSize * btilde6) * k6 + (timeStepSize * btilde7) * k7
                                                          ));
                denseTime += denseTimeStep;
            }
            
            currentY    = y5th;
            currentTime += timeStepSize;
                        
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
            
            
            // ==================
            // Iteration concluding tasks (FSAL)
            k1 = propagatorFast(currentTime , currentTime + timeStepSize, k7);
            phi0 = propagatorFast(currentTime, currentTime + timeStepSize, currentY);
                                
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


} // namespace Algorithms 

} // namespace Tortoise
#endif /* DormandPrince54IP_hpp */
