//
//  main.cpp
//  Workshop
//
//  Created by Marco Battiato on 15/2/23.
//

// #define NDEBUG            // Uncomment to deactivate debug mode

#include <TORTOISE>
#include <iostream>

using namespace Tortoise;
using std::cout;
using Tortoise::AnalyticExpression::x, Tortoise::AnalyticExpression::y, Tortoise::AnalyticExpression::z;
using Tortoise::AnalyticExpression::k;

int main(int argc, const char * argv[]) {
    
    Real                pi = 3.1415;
    Point<2>            originBZ({0.0, 0.}),
                        sideBZ0({1.0, 0.0}),
                        sideBZ1({-.5, std::sqrt(3.0)/2});
    Region<2>           brillouinZone(originBZ,{sideBZ0,sideBZ1});

    Material<2>         exMat("ExampleMaterial", brillouinZone);
    
    int                 resolutionFull = 15, resolutionPart = 7;
    Mesh<2>             meshFull(brillouinZone, {0.,0.},{1.,1.},{resolutionFull,resolutionFull});
    Mesh<2>             meshPart(brillouinZone, {0.3,0.3},{.4,.4},{resolutionPart,resolutionPart});
    
    Real                phononOpt0Freq = 0.1;
    exMat.addBand("pn:opt0", 0, boson, meshFull, phononOpt0Freq);
    Point<2>            ortsideBZ0({.0, 1.0}), ortsideBZ1({std::sqrt(3.0)/2, .5});
    exMat.addBand("el:V_A", -1, fermion, meshFull, 1.1+sin(0.5*pi*dot(ortsideBZ0,k - 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2))*sin(0.5*pi*dot(ortsideBZ1,k- 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2)));
    exMat.addBand("el:V_B", -1, fermion, meshPart, 2.5 + 9*(x-0.25)*(x-0.25) + 9*(y-0.45)*(y-0.45));
    
    exMat.plot();
    
    exMat["el:V_B"].dispersion.plot("Band's dispersion plot");
    
    int                 resolutionPhoton = 20;
    Mesh<2>             photonMesh = exMat.photonMesh(resolutionPhoton);
    Real                photonLowestEnergy = 0, photonHighestEnergy = 3.5;
    exMat.addPhotonBand("pt:ppol", photonMesh, photonLowestEnergy, photonHighestEnergy);
    
    exMat.plot();
    
    cout << exMat;
    
    return 0;
}
