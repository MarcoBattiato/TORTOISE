//
//  GaussQuadratureConst.hpp
//  TORTOISE_Development
//
//  Created by Marco Battiato on 2/11/22.
//

#ifndef GaussQuadratureConst_hpp
#define GaussQuadratureConst_hpp

#include <stdio.h>

namespace Tortoise {

namespace GeometryCore {

//*******************************
// Gauss quadrature
//*******************************
// The gauss points and weights are for reference elements
// The reference element in 1D is the [0,1] interval
// The reference element in 2D is the {[0,0],[1,0],[0,1]} triangle
// The reference element in 3D is the {[0,0,0],[1,0,0],[0,1,0],[0,0,1]} tetrahedron
constexpr int ngausspoint(const int dim, const int ord) {
    switch (dim) {
        case 1: return ord; break;
        case 2: switch (ord) {
                    case 1: return 1; break;
                    case 2: return 3; break;
                    case 3: return 4; break;
                    case 4: return 6; break;
                    case 5: return 7; break;
                    case 6: return 12; break;
                    case 7: return 13; break;
                    default: return -1; break;
                }
        case 3: switch (ord) {
                    case 1: return 1; break;
                    case 2: return 4; break;
                    case 3: return 5; break;
                    case 4: return 10; break;
                    case 5: return 11; break;
                    case 6: return 16; break;
                    case 7: return 13; break;
                    default: return -1; break;
            }
        default: return -1; break;
    }
}

template <int NDim, int order> const ArrayPoint<NDim,ngausspoint(NDim,order)> gaussPoints;
template <int NDim, int order> const ArrayPoint<1,ngausspoint(NDim,order)> gaussWeigths;

// 1D
// ==============
template<> const inline ArrayPoint<1,1> gaussPoints<1,1> {0.5};
template<> const inline ArrayPoint<1,1> gaussWeigths<1,1> {1.0};
// -----
template<> const inline ArrayPoint<1,2> gaussPoints<1,2> {0.2113248654051871,0.7886751345948129};
template<> const inline ArrayPoint<1,2> gaussWeigths<1,2> {0.5,0.5};
// -----
template<> const inline ArrayPoint<1,3> gaussPoints<1,3> {0.1127016653792583,0.5,0.8872983346207417};
template<> const inline ArrayPoint<1,3> gaussWeigths<1,3> {0.2777777777777778,0.4444444444444444,0.2777777777777778};
// -----
template<> const inline ArrayPoint<1,4> gaussPoints<1,4> {0.0694318442029737,0.3300094782075719,0.6699905217924281,0.9305681557970263};
template<> const inline ArrayPoint<1,4> gaussWeigths<1,4> {0.1739274225687269,0.3260725774312731,0.3260725774312731,0.1739274225687269};
// -----
template<> const inline ArrayPoint<1,5> gaussPoints<1,5> = [] {ArrayPoint<1,5> tmp;
    tmp << 0.0469100770306680,0.0230765344947158,0.5,0.7692346550528415,0.9530899229693319;
    return tmp; }();
template<> const inline ArrayPoint<1,5> gaussWeigths<1,5> = [] {ArrayPoint<1,5> tmp;
    tmp << 0.1184634425280945,0.2393143352496832,0.2844444444444444,0.2393143352496832,0.1184634425280945;
    return tmp;}();
// 2D
// ==============
template<> const inline ArrayPoint<2,ngausspoint(2,1)> gaussPoints<2,1> = [] {ArrayPoint<2,ngausspoint(2,1)> tmp;
tmp << 0.33333333333333, 0.33333333333333;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(2,1)> gaussWeigths<2,1> = [] {ArrayPoint<1,ngausspoint(2,1)> tmp;
tmp << 1.00000000000000;
return tmp; }();
// -----
template<> const inline ArrayPoint<2,ngausspoint(2,2)> gaussPoints<2,2> = [] {ArrayPoint<2,ngausspoint(2,2)> tmp;
tmp << 0.16666666666667, 0.66666666666667, 0.16666666666667,
    0.16666666666667, 0.16666666666667, 0.66666666666667;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(2,2)> gaussWeigths<2,2> = [] {ArrayPoint<1,ngausspoint(2,2)> tmp;
tmp << 0.33333333333333, 0.33333333333333, 0.33333333333333;
return tmp; }();
// -----
template<> const inline ArrayPoint<2,ngausspoint(2,3)> gaussPoints<2,3> = [] {ArrayPoint<2,ngausspoint(2,3)> tmp;
tmp << 0.20000000000000, 0.60000000000000, 0.20000000000000, 0.33333333333333,
    0.20000000000000, 0.20000000000000, 0.60000000000000, 0.33333333333333;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(2,3)> gaussWeigths<2,3> = [] {ArrayPoint<1,ngausspoint(2,3)> tmp;
tmp << 0.52083333333333, 0.52083333333333, 0.52083333333333, -0.5625000000000;
return tmp; }();
// -----
template<> const inline ArrayPoint<2,ngausspoint(2,4)> gaussPoints<2,4> = [] {ArrayPoint<2,ngausspoint(2,4)> tmp;
tmp << 0.09157621350977, 0.44594849091597, 0.81684757298046,  0.44594849091597, 0.09157621350977, 0.10810301816807,
    0.09157621350977, 0.10810301816807, 0.09157621350977, 0.44594849091597, 0.81684757298046, 0.44594849091597;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(2,4)> gaussWeigths<2,4> = [] {ArrayPoint<1,ngausspoint(2,4)> tmp;
    tmp << 0.10995174365532, 0.22338158967801, 0.10995174365532, 0.22338158967801, 0.10995174365532, 0.22338158967801;
return tmp; }();
// -----
template<> const inline ArrayPoint<2,ngausspoint(2,5)> gaussPoints<2,5> = [] {ArrayPoint<2,ngausspoint(2,5)> tmp;
tmp << 0.10128650732346, 0.47014206410511, 0.79742698535309, 0.47014206410511, 0.10128650732346, 0.05971587178977, 0.33333333333333,
    0.10128650732346, 0.05971587178977, 0.10128650732346, 0.47014206410511, 0.79742698535309, 0.47014206410511, 0.33333333333333;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(2,5)> gaussWeigths<2,5> = [] {ArrayPoint<1,ngausspoint(2,5)> tmp;
    tmp << 0.12593918054483, 0.13239415278851, 0.12593918054483, 0.13239415278851, 0.12593918054483, 0.13239415278851, 0.22500000000000;
return tmp; }();
// 3D
// ==============
template<> const inline ArrayPoint<3,ngausspoint(3,1)> gaussPoints<3,1> = [] {ArrayPoint<3,ngausspoint(3,1)> tmp;
tmp << 0.250000000000,
    0.2500000000000,
    0.2500000000000;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(3,1)> gaussWeigths<3,1> = [] {ArrayPoint<1,ngausspoint(3,1)> tmp;
tmp << 1.00000000000000;
return tmp; }();
// -----
template<> const inline ArrayPoint<3,ngausspoint(3,2)> gaussPoints<3,2> = [] {ArrayPoint<3,ngausspoint(3,2)> tmp;
tmp << 0.585410196624, 0.138196601125, 0.138196601125, 0.138196601125,
    0.1381966011250, 0.5854101966249, 0.1381966011250, 0.1381966011250,
    0.1381966011250, 0.1381966011250, 0.5854101966249, 0.1381966011250;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(3,2)> gaussWeigths<3,2> = [] {ArrayPoint<1,ngausspoint(3,2)> tmp;
    tmp << 0.25000000000000, 0.25000000000000, 0.25000000000000, 0.25000000000000;
return tmp; }();
// -----
template<> const inline ArrayPoint<3,ngausspoint(3,3)> gaussPoints<3,3> = [] {ArrayPoint<3,ngausspoint(3,3)> tmp;
    tmp << 0.250000000000, 0.500000000000, 0.166666666666, 0.166666666666, 0.166666666666,
    0.2500000000000, 0.1666666666666, 0.1666666666666, 0.1666666666666, 0.5000000000000,
    0.2500000000000, 0.1666666666666, 0.1666666666666, 0.5000000000000, 0.1666666666666;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(3,3)> gaussWeigths<3,3> = [] {ArrayPoint<1,ngausspoint(3,3)> tmp;
    tmp << -0.8000000000000, 0.45000000000000, 0.45000000000000, 0.45000000000000, 0.45000000000000;
return tmp; }();
// -----
template<> const inline ArrayPoint<3,ngausspoint(3,4)> gaussPoints<3,4> = [] {ArrayPoint<3,ngausspoint(3,4)> tmp;
    tmp << 0.568430584196, 0.143856471934, 0.143856471934, 0.143856471934, 0.000000000000, 0.500000000000, 0.500000000000, 0.500000000000, 0.000000000000, 0.000000000000,
    0.1438564719343, 0.1438564719343, 0.1438564719343, 0.5684305841968, 0.5000000000000, 0.0000000000000, 0.5000000000000, 0.0000000000000, 0.5000000000000, 0.0000000000000,
    0.1438564719343, 0.1438564719343, 0.5684305841968, 0.1438564719343, 0.5000000000000, 0.5000000000000, 0.0000000000000, 0.0000000000000, 0.0000000000000, 0.5000000000000;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(3,4)> gaussWeigths<3,4> = [] {ArrayPoint<1,ngausspoint(3,4)> tmp;
    tmp << 0.21776506988041, 0.21776506988041, 0.21776506988041, 0.21776506988041, 0.02148995341306, 0.02148995341306, 0.02148995341306, 0.02148995341306, 0.02148995341306, 0.02148995341306;
return tmp; }();
// -----
template<> const inline ArrayPoint<3,ngausspoint(3,5)> gaussPoints<3,5> = [] {ArrayPoint<3,ngausspoint(3,5)> tmp;
    tmp << 0.250000000000, 0.785714285714, 0.071428571428, 0.071428571428, 0.071428571428, 0.100596423833, 0.399403576166, 0.399403576166, 0.399403576166, 0.100596423833, 0.100596423833,
    0.2500000000000, 0.0714285714285, 0.0714285714285, 0.0714285714285, 0.7857142857142, 0.3994035761668, 0.1005964238332, 0.3994035761668, 0.1005964238332, 0.3994035761668, 0.1005964238330,
    0.2500000000000, 0.0714285714285, 0.0714285714285, 0.7857142857142, 0.0714285714285, 0.3994035761668, 0.3994035761668, 0.1005964238332, 0.1005964238332, 0.1005964238332, 0.3994035761668;
return tmp; }();
template<> const inline ArrayPoint<1,ngausspoint(3,5)> gaussWeigths<3,5> = [] {ArrayPoint<1,ngausspoint(3,5)> tmp;
    tmp << -0.0789333333333, 0.04573333333333, 0.04573333333333, 0.04573333333333, 0.04573333333333, 0.14933333333333, 0.14933333333333, 0.14933333333333, 0.14933333333333, 0.14933333333333, 0.14933333333333;
return tmp; }();

} // namespace GeometryCore

} // namespace Tortoise 

#endif /* GaussQuadratureConst_hpp */
