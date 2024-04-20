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
//  IntegratorCore.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 18/6/21.
//
// The class builds an object capable of integrating
// (\int_T)^NLegs \delta(D K - N) * f1(K) * f2(K) * ... dK
// where:
//      K is the vector of all the coordinates per leg. Its dimensionality is NDim*NLegs
//      The coordinates are expected to be grouped by leg
//          [ k1 ]
//      K = [ k2 ]
//          [ k3 ]
//      T is the reference tetrahedron in NDim dimensions
//             i.e. is interval (0,1) in 1D; triangle ([0,0],[1,0],[0,1]) in 2D; tetrahedron ([0,0,0],[1,0,0],[0,1,0], [0,0,1]) in 3D etc.
//      NLegs is the number of NDim-D spaces which are multiplied in cartesian sense and integrated over
//
//      The Dirac delta is specified as a linear Dirac delta, where
//          D is a [NDim+1,NLegs*NDim]-dimensional matrix
//          N is a NDim+1 dimensional vector
//      Notice that this does not allow for the construction of a generic Dirac delta as the number of contraints is always NDim+1
//
//      f1(K), f2(K), ... are polynomial functions of the full coordinates K
//
//
// We realise that NDim+1 of the componets of the K vector are fixed by the Dirac delta.
// We implement here a quick version where the constrained variables are always the last ones.
// This might deliver bad numerical results if the energy dispersion of the last two legs in the last direction are similar in values
// (for instance if the last two legs are the same band and we calcute the scattering to the same element for both legs).
// It is up to the user to send an appropriate integral to have numerically stable results
//
// We write the K vector as composed by two parts: the free coordinates and the constrained coordinates
//          [ kF ]
//      K = [    ]
//          [ kC ]
//
// kF has dimension (NLegs-1)*NDim-1
//          [       k0      ]
//     KF = [       k1      ]
//          [      ...      ]
//          [  k(NLegs-2)P  ]
// Where k(NLegs-2)P is a vector of k coordinates corresponding to the (NLegs-1)-th leg, but where the last coordinate is absent (because it is constrained)
// Notice that in 1D this vector has 0 components
//
// kC has dimension NDim+1
//          [   k(NLegs-2)(Ndim-1)  ]
//     KC = [       k(NLegs-1)      ]
// Where the k(NLegs-2)(Ndim-1) is the last coordinate of the (NLegs-1)-th leg. Notice that this is a single number.
// The second part k(NLegs-1) is the k coordinates of the last leg.
//
// Similarly we split the D matrix in two parts
// D = [DF DC]
// splitting the parts that multiply the different parts of the K vector.
//
// So we can analytically integrate
// (\int_T)^NLegs \delta( D K - N) * f1(K) * f2(K) * ... dK =
//   = (\int_T)^NLegs \delta( DF KF + DC KC - N) * f1(K) * f2(K) * ... dK =
//   = abs(1/Det(DC))  (\int_T)^(NLegs-2) (\int_T')  \Theta_T(k(NLegs-2)) * \Theta_T(k(NLegs-1)) * f1(K) * f2(K) * ... dKF
// where:
//      Det(DC) is the determinant of DC
//      \Theta_T(k(NLegs-2)) is a Theta function that tests if the coordinates of the (NLegs-1)-th leg is within the reference tetrahedron
//      \Theta_T(k(NLegs-1)) is a Theta function that tests if the coordinates of the NLegs-th leg is within the reference tetrahedron
//
// We can do a Monte Carlo integration of this formula
//
// First we need to find the volume of the integration.
// The volume of an n-Dim tetrahedron is 1/n!.
// So the volume of integration is
//   VolIntegr = (1/NDim!)^(NLegs-2) * (1/(NDim-1)!)
//
// Now we need to contruct the a MC point (i.e. a KF)
//          [       k0      ]
//     KF = [       k1      ]
//          [      ...      ]
//          [  k(NLegs-2)P  ]
// We need to extract k0, k1, ..., k(NLegs-3) each one uniformly in the reference tetrahedron in NDim dimansions. I do not explain here how to do it, but the method is known.
// We need to extract k(NLegs-2)P uniformly in the reference tetrahedron in NDim-1 dimansions. Same method as above but in 1 dimension less, or not constructed at all if NDim = 1.
//
// We will actually construct several MC points and we will put them in a matrix
//           [       k0             k0             k0            ...             k0      ]
//     KFM = [       k1             k1             k1            ...             k1      ]
//           [      ...            ...            ...            ...            ...      ]
//           [  k(NLegs-2)P    k(NLegs-2)P    k(NLegs-2)P        ...        k(NLegs-2)P  ]
// with dimensions (NLegs-1*NDim)-1 by NMC (where NMC is the number of Monte Carlo points)
//
// We compute KC
//    KC = - DC^-1 ( DF KF - N )
// We can construct the matrix
//    KCM = - DC^-1 ( DF KFM - NM )
//  where NM is simply the repetion NMC times on N
//
// We construct K0M, K1M,  KM(NLegs-1) and KM
//           [ kFM ]
//      KM = [     ]
//           [ kCM ]
//
// We now need to construct all the functions to integrate which will be vectors 1 by NMC
//  \Theta_T(k(NLegs-1))
//  \Theta_T(k(NLegs-2))
//  f1(K)
//
// We finally multiply all these vectors and take the total
//


#ifndef IntegratorCore_hpp
#define IntegratorCore_hpp

#include <Geometry/GeometryCore/Geometry.hpp>

#include <stdio.h>
#include <cassert>
#include <iostream>

using namespace Tortoise::GeometryCore;

namespace Tortoise {

namespace PhysicsCore {

template <int NDim, int NLegs, int NForms> using ArrayMultiLegLinearForm       = Eigen::Matrix<Real, NForms, NLegs*(NDim+1)>;

template<int NDIMS, int NLEGS> Real integrationInversionDeterminant(const Eigen::Matrix<Real, NDIMS+1,NLEGS*NDIMS> D){
    return D.template block<NDIMS+1,NDIMS+1>(0,(NLEGS-1)*NDIMS-1).determinant();
};

// Note for self: All the results are returned by value. Also they are not returned using "auto", which in principle would have allowed to return an expression.
// However when I tried to return auto (to be able to return an expression, rather than the matrix values), I kept getting compilation errors. I suspect that somewhere
// during the lazy evaluation, one of the parts of the expression that are needed for the evaluation have gone out of scope.


//template<int NDim, int NLEGS, typename DerivedA, typename DerivedB, typename DerivedC, typename DerivedD, typename DerivedE>
//Eigen::Matrix<Real, NLEGS*(NDim+1), 1>                     // Expected -> Eigen::Matrix<Real, NLEGS*(NDim+1), 1>
//scatteringIntegrationTypeB
//(const Eigen::MatrixBase<DerivedA>& D,                              // Expected -> Eigen::Matrix<Real, NDim+1,NLEGS*NDim>
// const Eigen::MatrixBase<DerivedB>& N,                              // Expected -> Eigen::Matrix<Real,NDim+1, 1>
// const Eigen::MatrixBase<DerivedC>& commonLinFormA,                 // Expected -> ArrayMultiLegLinearForm<NDim,NLEGS,1>
// const Eigen::MatrixBase<DerivedD>& commonLinFormB,                 // Expected -> ArrayMultiLegLinearForm<NDim,NLEGS,1>
// const Eigen::MatrixBase<DerivedE>& outputformsLinForm,             // Expected -> ArrayMultiLegLinearForm<NDIMS,NLEGS,NLEGS*(NDIMS+1)>
// const int NMCPoints) {
//    
//    assert( (D.rows() == NDim+1) && (D.cols() == NLEGS*NDim) );
//    assert( (N.rows() == NDim+1) && (N.cols() == 1) );
//    assert( (commonLinFormA.rows() == 1) && (commonLinFormA.cols() == NLEGS*(NDim+1)) );
//    assert( (commonLinFormB.rows() == 1) && (commonLinFormB.cols() == NLEGS*(NDim+1)) );
//    assert( (outputformsLinForm.rows() == NLEGS*(NDim+1)) && (outputformsLinForm.cols() == NLEGS*(NDim+1)) );
//
//    typedef Eigen::Matrix<Real, NDim, Eigen::Dynamic> LegMCPoints;
//    
//    // Construct DF and DC
//    auto D0F  (D.template block<NDim+1,NDim>(0,0));                      // Free part of D: Used only if NLEGS >= 3
//    auto D1F  (D.template block<NDim+1,NDim>(0,NDim));                   // Free part of D: Used only if NLEGS >= 4
//    auto DLpF (D.template block<NDim+1,NDim-1>(0,(NLEGS-2)*NDim));       // Free part of D: Partially free leg: corresponds to leg NLEGS-2
//    auto DC   (D.template block<NDim+1,NDim+1>(0,(NLEGS-1)*NDim-1));     // Constrained part of D present for any NLEGS
//    
//    // Free Coordinates
//    LegMCPoints K0(randomPointsReference<NDim>(NMCPoints));              // Unconstrained leg: Used only if NLEGS >= 3
//    LegMCPoints K1(randomPointsReference<NDim>(NMCPoints));              // Unconstrained leg: Used only if NLEGS >= 4
//    LegMCPoints KpF;                                     // Partially constrained leg: corresponds to leg NLEGS-2
//    if constexpr((NLEGS == 2) && (NDim == 1)){ KpF.resize(NDim, 1);
//    } else { KpF.resize(NDim, NMCPoints); }
//    if constexpr(NDim>1){
//        KpF.block(0, 0, NDim-1, NMCPoints) = randomPointsReference<NDim-1>(NMCPoints);
//    }
//    
//    // Constrained Coordinates
//    Eigen::Matrix<Real, NDim+1, Eigen::Dynamic> alpha;
//    if constexpr(NLEGS == 2){
//        if constexpr (NDim == 1){ alpha = DC.inverse() * N; }
//        else                    { alpha = - DC.inverse() * ( (DLpF * KpF.block(0, 0, NDim-1, NMCPoints) ).colwise() - N) ; }}
//    if constexpr(NLEGS == 3){
//        if constexpr (NDim == 1){ alpha = - DC.inverse() * ( (D0F * K0 ).colwise() - N) ; }
//        else                    { alpha = - DC.inverse() * ( (D0F * K0 + DLpF * KpF.block(0, 0, NDim-1, NMCPoints) ).colwise() - N) ; }}
//    if constexpr(NLEGS == 4){
//        if constexpr (NDim == 1){ alpha = - DC.inverse() * ( (D0F * K0 + D1F * K1 ).colwise() - N) ; }
//        else                    { alpha = - DC.inverse() * ( (D0F * K0 + D1F * K1 +  DLpF * KpF.block(0, 0, NDim-1, NMCPoints) ).colwise() - N) ; }}
//        
//    LegMCPoints KC;
//    if constexpr( (NLEGS == 2) && (NDim == 1)){
//        KpF.block(NDim-1, 0, 1, 1) = alpha.block(0, 0, 1, 1);
//        KC = alpha.block(1, 0, NDim, 1);                  // Fully constrained leg: corresponds to leg NLEGS-1
//    } else {
//        KpF.block(NDim-1, 0, 1, NMCPoints) = alpha.block(0, 0, 1, NMCPoints);
//        KC = alpha.block(1, 0, NDim, NMCPoints);                  // Fully constrained leg: corresponds to leg NLEGS-1
//    }
//
//    
//    
//    // Theta functions
//    auto ThetapC ( (((KpF.array()>0.0).colwise().all() && (KpF.colwise().sum().array()<1.0)).template cast<Real>()) );
//    auto ThetaFC ( (((KC.array()>0.0).colwise().all() && (KC.colwise().sum().array()<1.0)).template cast<Real>()) );
//    
//    auto evalCommonA0  (apply(commonLinFormA.template block<1,NDim+1>(0,0),                  K0).array());    // Used only if NLEGS >= 3
//    auto evalCommonA1  (apply(commonLinFormA.template block<1,NDim+1>(0,NDim+1),             K1).array());    // Used only if NLEGS >= 4
//    auto evalCommonApF (apply(commonLinFormA.template block<1,NDim+1>(0,(NLEGS-2)*(NDim+1)), KpF).array());   //
//    auto evalCommonAC  (apply(commonLinFormA.template block<1,NDim+1>(0,(NLEGS-1)*(NDim+1)), KC).array());    //
//    auto evalCommonB0  (apply(commonLinFormB.template block<1,NDim+1>(0,0),                  K0).array());    // Used only if NLEGS >= 3
//    auto evalCommonB1  (apply(commonLinFormB.template block<1,NDim+1>(0,NDim+1),             K1).array());    // Used only if NLEGS >= 4
//    auto evalCommonBpF (apply(commonLinFormB.template block<1,NDim+1>(0,(NLEGS-2)*(NDim+1)), KpF).array());   //
//    auto evalCommonBC  (apply(commonLinFormB.template block<1,NDim+1>(0,(NLEGS-1)*(NDim+1)), KC).array());    //
//
//    auto evalOutput0  (apply(outputformsLinForm.template block<NLEGS*(NDim+1),NDim+1>(0,0),                  K0).array()); // Used only if NLEGS >= 3
//    auto evalOutput1  (apply(outputformsLinForm.template block<NLEGS*(NDim+1),NDim+1>(0,NDim+1),             K1).array()); // Used only if NLEGS >= 4
//    auto evalOutputpF (apply(outputformsLinForm.template block<NLEGS*(NDim+1),NDim+1>(0,(NLEGS-2)*(NDim+1)), KpF).array());
//    auto evalOutputC  (apply(outputformsLinForm.template block<NLEGS*(NDim+1),NDim+1>(0,(NLEGS-1)*(NDim+1)), KC).array());
//    
//    if constexpr(NLEGS == 2){
//        if constexpr(NDim == 1){
//            return (1./(std::fabs(DC.determinant()) ) ) *
//                     ( (evalOutputpF * evalOutputC).rowwise() *
//                      ( (evalCommonApF * evalCommonAC + evalCommonBpF * evalCommonBC) * ThetapC.array() * ThetaFC.array() )).rowwise().sum();
//        } else {
//            return (1./(static_cast<double>(Utilities::factorial(NDim-1)) * static_cast<double>(NMCPoints) * std::fabs(DC.determinant()) ) ) *
//                     ( (evalOutputpF * evalOutputC).rowwise() *
//                      ( (evalCommonApF * evalCommonAC + evalCommonBpF * evalCommonBC) * ThetapC.array() * ThetaFC.array() )).rowwise().sum();
//        }
//    }
//    if constexpr(NLEGS == 3){
//        return (1./( static_cast<double>(Utilities::factorial(NDim)) * static_cast<double>(Utilities::factorial(NDim-1)) * static_cast<double>(NMCPoints) * std::fabs(DC.determinant()) ) ) *
//                 ( (evalOutput0 * evalOutputpF * evalOutputC).rowwise() *
//                  ( (evalCommonA0 *  evalCommonApF * evalCommonAC + evalCommonB0 *  evalCommonBpF * evalCommonBC) * ThetapC.array() * ThetaFC.array() )).rowwise().sum();
//    }
//    if constexpr(NLEGS == 4){
//        return (1./( static_cast<double>(Utilities::factorial(NDim)) * static_cast<double>(Utilities::factorial(NDim)) * static_cast<double>(Utilities::factorial(NDim-1)) * static_cast<double>(NMCPoints) * std::fabs(DC.determinant()) ) ) *
//                 ( (evalOutput0 * evalOutput1 * evalOutputpF * evalOutputC).rowwise() *
//                  ( (evalCommonA0 * evalCommonA1 * evalCommonApF * evalCommonAC + evalCommonB0 * evalCommonB1 * evalCommonBpF * evalCommonBC) * ThetapC.array() * ThetaFC.array() )).rowwise().sum();
//    }
//}

template<int NDim, int NLEGS, typename DerivedA, typename DerivedB, typename DerivedC, typename DerivedD, typename DerivedE, typename DerivedF>
Eigen::Matrix<Real, NLEGS*(NDim+1), 1>                     // Expected -> Eigen::Matrix<Real, NLEGS*(NDim+1), 1>
scatteringIntegrationTypeB
(const Eigen::MatrixBase<DerivedA>& D,                              // Expected -> Eigen::Matrix<Real, NDim+1,NLEGS*NDim>
 const Eigen::MatrixBase<DerivedB>& N,                              // Expected -> Eigen::Matrix<Real,NDim+1, 1>
 const Eigen::MatrixBase<DerivedC>& commonLinFormA,                 // Expected -> ArrayMultiLegLinearForm<NDim,NLEGS,1>
 const Eigen::MatrixBase<DerivedD>& commonLinFormB,                 // Expected -> ArrayMultiLegLinearForm<NDim,NLEGS,1>
 const Eigen::MatrixBase<DerivedE>& outputformsLinForm,             // Expected -> ArrayMultiLegLinearForm<NDIMS,NLEGS,NLEGS*(NDIMS+1)>
 const Eigen::MatrixBase<DerivedF>& MCPoints) {
    
    assert( (D.rows() == NDim+1) && (D.cols() == NLEGS*NDim) );
    assert( (N.rows() == NDim+1) && (N.cols() == 1) );
    assert( (commonLinFormA.rows() == 1) && (commonLinFormA.cols() == NLEGS*(NDim+1)) );
    assert( (commonLinFormB.rows() == 1) && (commonLinFormB.cols() == NLEGS*(NDim+1)) );
    assert( (outputformsLinForm.rows() == NLEGS*(NDim+1)) && (outputformsLinForm.cols() == NLEGS*(NDim+1)) );

    // Construct DF and DC
    auto D0F  (D.template block<NDim+1,NDim>(0,0));                      // Free part of D: Used only if NLEGS >= 3
    auto D1F  (D.template block<NDim+1,NDim>(0,NDim));                   // Free part of D: Used only if NLEGS >= 4
    auto DLpF (D.template block<NDim+1,NDim-1>(0,(NLEGS-2)*NDim));       // Free part of D: Partially free leg: corresponds to leg NLEGS-2
    auto DC   (D.template block<NDim+1,NDim+1>(0,(NLEGS-1)*NDim-1));     // Constrained part of D present for any NLEGS
    
    typedef Eigen::Matrix<Real, NDim, Eigen::Dynamic> LegMCPoints;
    
    // Free Coordinates
    auto K0 = MCPoints.template topRows<NDim>();                   // Unconstrained leg: Used only if NLEGS >= 3
    auto K1 = MCPoints.template middleRows<NDim>(NDim);            // Unconstrained leg: Used only if NLEGS >= 4
    LegMCPoints KpF;                                     // Partially constrained leg: corresponds to leg NLEGS-2
    if constexpr((NLEGS == 2) && (NDim == 1)){
        KpF.resize(NDim, 1);
    } else {
        KpF.resize(NDim, MCPoints.cols());
    }
    if constexpr(NDim>1){
        KpF.template topRows<NDim-1>() = MCPoints.template bottomRows<NDim-1>();
    }
    
    // Constrained Coordinates
    Eigen::Matrix<Real, NDim+1, Eigen::Dynamic> alpha;
    if constexpr(NLEGS == 2){
        if constexpr (NDim == 1){ alpha = DC.inverse() * N; }
        else                    { alpha = - DC.inverse() * ( (DLpF * KpF.template topRows<NDim-1>() ).colwise() - N) ; }}
    if constexpr(NLEGS == 3){
        if constexpr (NDim == 1){ alpha = - DC.inverse() * ( (D0F * K0 ).colwise() - N) ; }
        else                    { alpha = - DC.inverse() * ( (D0F * K0 + DLpF * KpF.template topRows<NDim-1>() ).colwise() - N) ; }}
    if constexpr(NLEGS == 4){
        if constexpr (NDim == 1){ alpha = - DC.inverse() * ( (D0F * K0 + D1F * K1 ).colwise() - N) ; }
        else                    { alpha = - DC.inverse() * ( (D0F * K0 + D1F * K1 +  DLpF * KpF.template topRows<NDim-1>() ).colwise() - N) ; }}
        
    LegMCPoints KC;
    if constexpr( (NLEGS == 2) && (NDim == 1)){
        KpF.block(NDim-1, 0, 1, 1) = alpha.block(0, 0, 1, 1);
        KC = alpha.block(1, 0, NDim, 1);                  // Fully constrained leg: corresponds to leg NLEGS-1
    } else {
        KpF.template bottomRows<1>() = alpha.template topRows<1>();
        KC = alpha.template bottomRows<NDim>();                  // Fully constrained leg: corresponds to leg NLEGS-1
    }

    // Theta functions
    auto ThetapC ( (((KpF.array()>0.0).colwise().all() && (KpF.colwise().sum().array()<1.0)).template cast<Real>()) );
    auto ThetaFC ( (((KC.array()>0.0).colwise().all() && (KC.colwise().sum().array()<1.0)).template cast<Real>()) );
    
    auto evalCommonA0  (apply(commonLinFormA.template block<1,NDim+1>(0,0),                  K0).array());    // Used only if NLEGS >= 3
    auto evalCommonA1  (apply(commonLinFormA.template block<1,NDim+1>(0,NDim+1),             K1).array());    // Used only if NLEGS >= 4
    auto evalCommonApF (apply(commonLinFormA.template block<1,NDim+1>(0,(NLEGS-2)*(NDim+1)), KpF).array());   //
    auto evalCommonAC  (apply(commonLinFormA.template block<1,NDim+1>(0,(NLEGS-1)*(NDim+1)), KC).array());    //
    auto evalCommonB0  (apply(commonLinFormB.template block<1,NDim+1>(0,0),                  K0).array());    // Used only if NLEGS >= 3
    auto evalCommonB1  (apply(commonLinFormB.template block<1,NDim+1>(0,NDim+1),             K1).array());    // Used only if NLEGS >= 4
    auto evalCommonBpF (apply(commonLinFormB.template block<1,NDim+1>(0,(NLEGS-2)*(NDim+1)), KpF).array());   //
    auto evalCommonBC  (apply(commonLinFormB.template block<1,NDim+1>(0,(NLEGS-1)*(NDim+1)), KC).array());    //

    auto evalOutput0  (apply(outputformsLinForm.template block<NLEGS*(NDim+1),NDim+1>(0,0),                  K0).array()); // Used only if NLEGS >= 3
    auto evalOutput1  (apply(outputformsLinForm.template block<NLEGS*(NDim+1),NDim+1>(0,NDim+1),             K1).array()); // Used only if NLEGS >= 4
    auto evalOutputpF (apply(outputformsLinForm.template block<NLEGS*(NDim+1),NDim+1>(0,(NLEGS-2)*(NDim+1)), KpF).array());
    auto evalOutputC  (apply(outputformsLinForm.template block<NLEGS*(NDim+1),NDim+1>(0,(NLEGS-1)*(NDim+1)), KC).array());
    
    if constexpr(NLEGS == 2){
        if constexpr(NDim == 1){
            return (1./(std::fabs(DC.determinant()) ) ) *
                     ( (evalOutputpF * evalOutputC).rowwise() *
                      ( (evalCommonApF * evalCommonAC + evalCommonBpF * evalCommonBC) * ThetapC.array() * ThetaFC.array() )).rowwise().sum();
        } else {
            return (1./(static_cast<double>(Utilities::factorial(NDim-1)) * static_cast<double>(MCPoints.cols()) * std::fabs(DC.determinant()) ) ) *
                     ( (evalOutputpF * evalOutputC).rowwise() *
                      ( (evalCommonApF * evalCommonAC + evalCommonBpF * evalCommonBC) * ThetapC.array() * ThetaFC.array() )).rowwise().sum();
        }
    }
    if constexpr(NLEGS == 3){
        return (1./( static_cast<double>(Utilities::factorial(NDim)) * static_cast<double>(Utilities::factorial(NDim-1)) * static_cast<double>(MCPoints.cols()) * std::fabs(DC.determinant()) ) ) *
                 ( (evalOutput0 * evalOutputpF * evalOutputC).rowwise() *
                  ( (evalCommonA0 *  evalCommonApF * evalCommonAC + evalCommonB0 *  evalCommonBpF * evalCommonBC) * ThetapC.array() * ThetaFC.array() )).rowwise().sum();
    }
    if constexpr(NLEGS == 4){
        return (1./( static_cast<double>(Utilities::factorial(NDim)) * static_cast<double>(Utilities::factorial(NDim)) * static_cast<double>(Utilities::factorial(NDim-1)) * static_cast<double>(MCPoints.cols()) * std::fabs(DC.determinant()) ) ) *
                 ( (evalOutput0 * evalOutput1 * evalOutputpF * evalOutputC).rowwise() *
                  ( (evalCommonA0 * evalCommonA1 * evalCommonApF * evalCommonAC + evalCommonB0 * evalCommonB1 * evalCommonBpF * evalCommonBC) * ThetapC.array() * ThetaFC.array() )).rowwise().sum();
    }
}

template<int NDim, int NLEGS, typename DerivedA, typename DerivedB, typename DerivedC, typename DerivedD>
Eigen::Matrix<Real, NLEGS*(NDim+1), Eigen::Dynamic>                 // Expected -> Eigen::Matrix<Real, NLEGS*(NDim+1), 1>
scatteringIntegrationTypeC
(const Eigen::MatrixBase<DerivedA>& D,                              // Expected -> Eigen::Matrix<Real, NDim+1,NLEGS*NDim>
 const Eigen::MatrixBase<DerivedB>& N,                              // Expected -> Eigen::Matrix<Real,NDim+1, 1>
 const Eigen::MatrixBase<DerivedC>& commonLinForms,                 // Expected -> ArrayMultiLegLinearForm<NDim,NLEGS, >>Variable<<  >
 const Eigen::MatrixBase<DerivedD>& outputformsLinForm,             // Expected -> ArrayMultiLegLinearForm<NDIMS,NLEGS,NLEGS*(NDIMS+1)>
 const int NMCPoints) {
    
    assert( (D.rows() == NDim+1) && (D.cols() == NLEGS*NDim) );
    assert( (N.rows() == NDim+1) && (N.cols() == 1) );
    assert( commonLinForms.cols() == NLEGS*(NDim+1) );
    assert( (outputformsLinForm.rows() == NLEGS*(NDim+1)) && (outputformsLinForm.cols() == NLEGS*(NDim+1)) );

    typedef Eigen::Matrix<Real, NDim, Eigen::Dynamic> LegMCPoints;
    
    // Construct DF and DC
    auto D0F  (D.template block<NDim+1,NDim>(0,0));                      // Free part of D: Used only if NLEGS >= 3
    auto D1F  (D.template block<NDim+1,NDim>(0,NDim));                   // Free part of D: Used only if NLEGS >= 4
    auto DLpF (D.template block<NDim+1,NDim-1>(0,(NLEGS-2)*NDim));       // Free part of D: Partially free leg: corresponds to leg NLEGS-2
    auto DC   (D.template block<NDim+1,NDim+1>(0,(NLEGS-1)*NDim-1));     // Constrained part of D present for any NLEGS
    
    // Free Coordinates
    LegMCPoints K0(randomPointsReference<NDim>(NMCPoints));              // Unconstrained leg: Used only if NLEGS >= 3
    LegMCPoints K1(randomPointsReference<NDim>(NMCPoints));              // Unconstrained leg: Used only if NLEGS >= 4
    LegMCPoints KpF;                                                     // Partially constrained leg: corresponds to leg NLEGS-2
    if constexpr((NLEGS == 2) && (NDim == 1)){ KpF.resize(NDim, 1);
    } else { KpF.resize(NDim, NMCPoints); }
    if constexpr(NDim>1){
        KpF.block(0, 0, NDim-1, NMCPoints) = randomPointsReference<NDim-1>(NMCPoints);
    }
    
    // Constrained Coordinates
    Eigen::Matrix<Real, NDim+1, Eigen::Dynamic> alpha;
    if constexpr(NLEGS == 2){
        if constexpr (NDim == 1){ alpha = DC.inverse() * N; }
        else                    { alpha = - DC.inverse() * ( (DLpF * KpF.block(0, 0, NDim-1, NMCPoints) ).colwise() - N) ; }}
    if constexpr(NLEGS == 3){
        if constexpr (NDim == 1){ alpha = - DC.inverse() * ( (D0F * K0 ).colwise() - N) ; }
        else                    { alpha = - DC.inverse() * ( (D0F * K0 + DLpF * KpF.block(0, 0, NDim-1, NMCPoints) ).colwise() - N) ; }}
    if constexpr(NLEGS == 4){
        if constexpr (NDim == 1){ alpha = - DC.inverse() * ( (D0F * K0 + D1F * K1 ).colwise() - N) ; }
        else                    { alpha = - DC.inverse() * ( (D0F * K0 + D1F * K1 +  DLpF * KpF.block(0, 0, NDim-1, NMCPoints) ).colwise() - N) ; }}
        
    LegMCPoints KC;
    if constexpr( (NLEGS == 2) && (NDim == 1)){
        KpF.block(NDim-1, 0, 1, 1) = alpha.block(0, 0, 1, 1);
        KC = alpha.block(1, 0, NDim, 1);                  // Fully constrained leg: corresponds to leg NLEGS-1
    } else {
        KpF.block(NDim-1, 0, 1, NMCPoints) = alpha.block(0, 0, 1, NMCPoints);
        KC = alpha.block(1, 0, NDim, NMCPoints);                  // Fully constrained leg: corresponds to leg NLEGS-1
    }

    Eigen::Matrix<Real, NLEGS*(NDim+1), Eigen::Dynamic> toreturn(NLEGS*(NDim+1),commonLinForms.rows());
    
    // Theta functions
    auto ThetapC ( (((KpF.array()>0.0).colwise().all() && (KpF.colwise().sum().array()<1.0)).template cast<Real>()) );
    auto ThetaFC ( (((KC.array()>0.0).colwise().all() && (KC.colwise().sum().array()<1.0)).template cast<Real>()) );
    
    auto evalOutput0  (apply(outputformsLinForm.template block<NLEGS*(NDim+1),NDim+1>(0,0),                  K0).array()); // Used only if NLEGS >= 3
    auto evalOutput1  (apply(outputformsLinForm.template block<NLEGS*(NDim+1),NDim+1>(0,NDim+1),             K1).array()); // Used only if NLEGS >= 4
    auto evalOutputpF (apply(outputformsLinForm.template block<NLEGS*(NDim+1),NDim+1>(0,(NLEGS-2)*(NDim+1)), KpF).array());
    auto evalOutputC  (apply(outputformsLinForm.template block<NLEGS*(NDim+1),NDim+1>(0,(NLEGS-1)*(NDim+1)), KC).array());

    for (int i=0; i<commonLinForms.rows(); ++i){
        auto evalCommonA0  (apply(commonLinForms.template block<1,NDim+1>(i,0),                  K0).array());    // Used only if NLEGS >= 3
        auto evalCommonA1  (apply(commonLinForms.template block<1,NDim+1>(i,NDim+1),             K1).array());    // Used only if NLEGS >= 4
        auto evalCommonApF (apply(commonLinForms.template block<1,NDim+1>(i,(NLEGS-2)*(NDim+1)), KpF).array());   //
        auto evalCommonAC  (apply(commonLinForms.template block<1,NDim+1>(i,(NLEGS-1)*(NDim+1)), KC).array());    //
        
        if constexpr(NLEGS == 2){
            if constexpr(NDim == 1){
                toreturn.col(i) = (1./(std::fabs(DC.determinant()) ) ) *
                         ( (evalOutputpF * evalOutputC).rowwise() *
                          ( (evalCommonApF * evalCommonAC) * ThetapC.array() * ThetaFC.array() )).rowwise().sum();
            } else {
                toreturn.col(i) =  (1./(static_cast<double>(Utilities::factorial(NDim-1)) * static_cast<double>(NMCPoints) * std::fabs(DC.determinant()) ) ) *
                         ( (evalOutputpF * evalOutputC).rowwise() *
                          ( (evalCommonApF * evalCommonAC) * ThetapC.array() * ThetaFC.array() )).rowwise().sum();
            }
        }
        if constexpr(NLEGS == 3){
            toreturn.col(i) =  (1./( static_cast<double>(Utilities::factorial(NDim)) * static_cast<double>(Utilities::factorial(NDim-1)) * static_cast<double>(NMCPoints) * std::fabs(DC.determinant()) ) ) *
                     ( (evalOutput0 * evalOutputpF * evalOutputC).rowwise() *
                      ( (evalCommonA0 *  evalCommonApF * evalCommonAC) * ThetapC.array() * ThetaFC.array() )).rowwise().sum();
        }
        if constexpr(NLEGS == 4){
            toreturn.col(i) =  (1./( static_cast<double>(Utilities::factorial(NDim)) * static_cast<double>(Utilities::factorial(NDim)) * static_cast<double>(Utilities::factorial(NDim-1)) * static_cast<double>(NMCPoints) * std::fabs(DC.determinant()) ) ) *
                     ( (evalOutput0 * evalOutput1 * evalOutputpF * evalOutputC).rowwise() *
                      ( (evalCommonA0 * evalCommonA1 * evalCommonApF * evalCommonAC) * ThetapC.array() * ThetaFC.array() )).rowwise().sum();
        }
    }
    return toreturn;
}

} // namespace PhysicsCore

} // namespace Tortoise

#endif /* IntegratorCore_hpp */
