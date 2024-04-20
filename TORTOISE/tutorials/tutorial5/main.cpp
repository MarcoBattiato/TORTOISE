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
//  main.cpp
//  Tutorial5
//
//  Created by Marco Battiato on 31/10/22.
//

// Please remember that the graphical output will work only if you have written
// the GNUPlot running directory in the Configuration file

#include <TORTOISE>
#include <iostream>

using namespace Tortoise;
using std::cout;

int main(int argc, const char * argv[]) {
    
    // Let us build the material
    Real                pi = 3.1415;
    Point<2>            originBZ({0.0, 0.}), sideBZ0({1.0, 0.0}), sideBZ1({-.5, std::sqrt(3.0)/2});
    Point<2>            ortsideBZ0({.0, 1.0}), ortsideBZ1({std::sqrt(3.0)/2, .5});
    Region<2>           brillouinZone (originBZ,{sideBZ0,sideBZ1});
    Material<2>         exMat ("ExampleMaterial", brillouinZone);
    
    int                 resolution = 15;
    Mesh<2>             meshFull(brillouinZone, {0.,0.},{1.,1.},{resolution,resolution});
    Mesh<2>             meshA(brillouinZone, {0.1,0.1},{0.3,0.3},{resolution,resolution});
    Mesh<2>             meshB(brillouinZone, {0.6,0.6},{0.3,0.3},{resolution,resolution});
    Mesh<2>             meshPhot(brillouinZone, {0.,0.},{0.001,0.0001},{5*resolution,1});
    
    Real                photMaxEnergy = 2.5;
    
    exMat.addBand("pn:ac0", 0, boson, meshFull, 0.05);
    exMat.addBand("el:V_A", -1, fermion, meshA, -1.1+sin(pi*dot(ortsideBZ0,k - 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2))*sin(pi*dot(ortsideBZ1,k- 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2)));
    exMat.addBand("el:C_A", -1, fermion, meshA, 1.1-sin(pi*dot(ortsideBZ0,k - 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2))*sin(pi*dot(ortsideBZ1,k- 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2)));
    exMat.addBand("pt:Lin", 0, boson, meshPhot, photMaxEnergy*dot(ortsideBZ1, k)/(0.001*std::sqrt(3.0)/2));
    
    for (auto scatt : exMat.addScatteringChannelsWildcards<3>({"el:*", "pn:*", "el:*"}, {in, in, out}, 2)){
        scatt->numMCPoints = 250; scatt->minimumDeterminant = 1.e-6; scatt->amplitude = 1;
    }
    for (auto scatt : exMat.addScatteringChannelsWildcards<3>({"el:*", "pt:*", "el:*"}, {in, in, out}, 1)) {
        scatt->numMCPoints = 250; scatt->minimumDeterminant = 1.e-6; scatt->amplitude = 1;
    }

    
    MaterialStatus<2>   equil(exMat);
    const Real chemPot = 0.0, temp = 0.1;
    for (auto band: exMat.matchListId("*")){
        if(exMat[band].statistics == fermion)
            equil[band] = 1/(exp((exMat[band].dispersion-chemPot)/temp)+1);
        else
            equil[band] = 1/(exp((exMat[band].dispersion-chemPot)/temp)-1);
    }
    
    // We now compute the time propagation under a laser pulse
    // The starting population will be in this case equilibrium, as the laser itself will trigger the non-equilibrium distribution
    MaterialTimeStatus<2>  evolution(equil);
    
    // We need to tell the time propagator that the population of the photon band should not be allowed to evolve freely.
    // That population is controlled by the femtosecond laser
    // We do this by specifying within the band object that the population is contrained
    
    exMat["pt:Lin"].constrained = true;
    
    // And specifying a lamda function that for any given time gives the desired population
    // In this case the population will be a gaussian in momentum around the desired central laser frequency and with the appropriate
    // laser energy width
    // The dependence on time will be a gaussian as well
    exMat["pt:Lin"].constrainedPopulation = [&exMat](Real t){return 1.e7*std::exp(-std::pow(t-0.3,2.)/0.1)*exp(-pow(exMat["pt:Lin"].dispersion-1.4,2)/0.04);};
    
    // The time evolution is calculated in the same way. TORTOISE will see that one of the bands is constrained and act accordingly
    stopWatch.tic();
    evolution.propagateDeterministic(5., 1.);
    cout << "Time taken for time propagation: " << stopWatch.toc() << "s\n";
    
    for (auto elBand : exMat.matchListId("el:*")){
        evolution.back()[elBand].plot(elBand);
    }

    for (int i =0; i<evolution.size(); ++i){
        plotter3d << evolution[i]["el:V_A"] << LABEL("time = "+toString(evolution.arg(i), 1));
    }
    plotter3d << PLOT;
    
    return 0;
}
