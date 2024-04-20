//
//  main.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 22/2/23.
//

// #define NDEBUG            // Uncomment to deactivate debug mode

#include <TORTOISE>
#include <iostream>
using namespace Tortoise;
using std::cout;
using Tortoise::AnalyticExpression::x, Tortoise::AnalyticExpression::y, Tortoise::AnalyticExpression::z;
using Tortoise::AnalyticExpression::k;

int main(int argc, const char * argv[]) {
    
    const Real          pi = 3.1415;
    const Region<1>     brillouinZone (-1,2);

    Material<1>         exMat ("1DMaterial", brillouinZone);
    
    int                 resolutionFull = 500, resolutionPart = 200;
    const Mesh<1>       meshFull(brillouinZone, 0, 1, resolutionFull);
    const Mesh<1>       meshPart(brillouinZone, 0.3, 0.4, resolutionPart);
    
    exMat.addBand("pn:ac", 0, boson, meshFull, 0.1*sin(abs(pi*x/2.))+0.05);
    exMat.addBand("el:C0", 0, fermion, meshPart, 15*x*x-0.3);
    
    for(auto scatt : exMat.addScatteringChannelsWildcards<3>({"el:*", "pn:*", "el:*"}, {in, in, out}, 2)){
        scatt->amplitude = 1.e-0;
        scatt->numMCPoints = 500;
    }
    
    exMat.plot();
    cout << exMat;
    
    MaterialStatus<1>   startPopul(exMat);
    Real temp = 0.05, chemPotFerm = 0., chemPotBos = 0;
    startPopul.setToEquilibrium(temp, chemPotFerm, chemPotBos);
    
    (1./exMat.scatteringRates<3>({"el:C0", "pn:ac", "el:C0"}, {in, in, out}, startPopul, "el:C0")).plot();
    
    startPopul["el:C0"] += 0.1*exp(-200.*(x+0.3)*(x+0.3));
    startPopul["el:C0"].plot();
    
    MaterialTimeStatus<1>   thermalDynamics(startPopul);
    thermalDynamics.propagateDeterministic(50, 2.);
    
    for(int i=0; i < thermalDynamics.size(); ++i){
        plotter2d << thermalDynamics[i]["el:C0"] << LABEL(std::to_string(thermalDynamics.arg(i)));
    }
    plotter2d << PLOTTITLE("Electron dynamics") << PLOT;
    
    for(int i=0; i < thermalDynamics.size(); ++i){
        plotter2d << thermalDynamics[i]["pn:ac"] << LABEL(std::to_string(thermalDynamics.arg(i)));
    }
    plotter2d << PLOTTITLE("Phonon dynamics") << PLOT;
    
    
    cout << "Tot Electron Density | Electron Energy Density | Phonon Energy Density | Total Energy Density | \n";
    for(int i=0; i < thermalDynamics.size(); ++i){
        cout << thermalDynamics[i]["el:C0"].integrate() << " | ";
        cout << thermalDynamics[i]["el:C0"].integrate(exMat["el:C0"].dispersion) << " | ";
        cout << thermalDynamics[i]["pn:ac"].integrate(exMat["pn:ac"].dispersion) << " | ";
        cout << thermalDynamics[i]["el:C0"].integrate(exMat["el:C0"].dispersion) + thermalDynamics[i]["pn:ac"].integrate(exMat["pn:ac"].dispersion) << " | ";
        cout << "\n";
    }
    
    
    
    return 0;
}

