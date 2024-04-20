//
//  main.cpp
//  TORTOISE
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
    Region<2>           brillouinZone (originBZ,{sideBZ0,sideBZ1});

    Material<2>         exMat ("ExampleMaterial", brillouinZone);
    
    int                 resolutionFull = 15, resolutionPart = 7;
    Mesh<2>             meshFull(brillouinZone, {0.,0.},{1.,1.},{resolutionFull,resolutionFull});
    Mesh<2>             meshPart(brillouinZone, {0.3,0.3},{.4,.4},{resolutionPart,resolutionPart});
    
    Real                phononOpt0Freq = 0.1;
    exMat.addBand("pn:opt0", 0, boson, meshFull, phononOpt0Freq);
    Point<2>            ortsideBZ0({.0, 1.0}), ortsideBZ1({std::sqrt(3.0)/2, .5});
    exMat.addBand("el:V_A", -1, fermion, meshFull, 1.1+sin(0.5*pi*dot(ortsideBZ0,k - 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2))*sin(0.5*pi*dot(ortsideBZ1,k- 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2)));
    exMat.addBand("el:V_B", -1, fermion, meshPart, 2.5 + 9*(x-0.25)*(x-0.25) + 9*(y-0.45)*(y-0.45));
    
    int                 resolutionPhoton = 20;
    Mesh<2>             photonMesh = exMat.photonMesh(resolutionPhoton);
    Real                photonLowestEnergy = 0, photonHighestEnergy = 3.5;
    exMat.addPhotonBand("pt:ppol", photonMesh, photonLowestEnergy, photonHighestEnergy);
    
    exMat.plot();
    
    MaterialStatus<2>   equilibrium(exMat);
    Real                kBT = .2, electronChemPot = 1, photonChemPot = 1;
    equilibrium["pn:opt0"] = 1/(exp(exMat["pn:opt0"].dispersion/kBT) - 1 );
    (equilibrium["el:V_A"] = 1/(exp((exMat["el:V_A"].dispersion-electronChemPot)/kBT) + 1 )).plot();
    
    equilibrium.setToEquilibrium(kBT, electronChemPot, photonChemPot);
    equilibrium["pt:ppol"] = 0;
   
    
    cout << " Electron density in band el:V_A  : " << equilibrium["el:V_A"].integrate(exMat["el:V_A"].densityOfStates) << "\n";
    cout << " Energy density in band el:V_A    : " << equilibrium["el:V_A"].integrate(exMat["el:V_A"].dispersion) << "\n";
    cout << " Momentum density in band el:V_A  : [" << equilibrium["el:V_A"].integrate(exMat["el:V_A"].crystalMomentum[0]) << "," << equilibrium["el:V_A"].integrate(exMat["el:V_A"].crystalMomentum[1]) << "]\n\n";
    
    cout << " Average group vel in band el:V_A : [" << equilibrium["el:V_A"].integrate(exMat["el:V_A"].dispersion.derivative(0)) << "," << equilibrium["el:V_A"].integrate(exMat["el:V_A"].dispersion.derivative(1)) << "] (it should be zero, but it is not because the band is not periodic, eh eh sorry)\n";
    
    
    // Advanced
    Region<1>           energyRange(photonLowestEnergy, photonHighestEnergy-photonLowestEnergy);
    Mesh<1>             energyMesh (energyRange, 0., 1., 50);
    Function<1>         densityOfStates(energyMesh);
    densityOfStates     = [&exMat](Point<1> en) {return exMat["el:V_A"].densityOfStates.integrateDiracDelta(exMat["el:V_A"].dispersion - en(0)); };
    
    Function<1>         populationInEnergy(energyMesh, [&exMat,&equilibrium](Point<1> en) {return equilibrium["el:V_A"].integrateDiracDelta(exMat["el:V_A"].dispersion - en(0)); });
    
    plotter2d << densityOfStates << LABEL("Density of states") << COLOR("black") << populationInEnergy << LABEL("Population") << COLOR("red") << PLOT;
    
    return 0;
}
