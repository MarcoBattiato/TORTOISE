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
//  Tutorial4
//
//  Created by Marco Battiato on 31/10/22.
//

#include <iostream>

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
        scatt->numMCPoints = 200; scatt->minimumDeterminant = 1.e-6; scatt->amplitude = 1;
    }
    
    MaterialStatus<2>   equil(exMat);
    const Real chemPot = 0.0, temp = 0.1;
    for (auto band: exMat.matchListId("*")){
        if(exMat[band].statistics == fermion)
            equil[band] = 1/(exp((exMat[band].dispersion-chemPot)/temp)+1);
        else
            equil[band] = 1/(exp((exMat[band].dispersion-chemPot)/temp)-1);
    }
    
    // We now want to compute the time propagation. This is computed internalli using Dormand-Prince 5(4)
    // DP54 is an adaptive step algorithm, which automatially adapts the time propagation step size in order
    // to maintain the error below a certain specified tolerance.
    
    // It is possible to specify two different types of tolerances: absolute and relative
    // Absolute tolerance is the absolute value of tolerated error (i.e. the numerical solution
    // is expected to be within the exact solution +- absolute tolerance)
    // The relative tolerance specifies the error relatively to the actual value of the solution
    // (i.e. the numerical solution is expected to be within the exact solution +- relative tolerance * solution)
    // In TORTOISE the above condition has been extended, since for fermions the difference of the population from
    // 1 is the number of holes. Therefore The time stepper requires the numerical solution to be within the
    // exact solution +- relative tolerance * solution or +- relative tolerance * (1-solution) in the case of fermions
    // which ever is the most restrictive.
    
    // It possible to define absolute and relative tolerances for each band independently. TORTOISE uses some default
    // values, unless specifically specified
    
    exMat["el:V_A"].absoluteErrorTolerance = 0.01;
    exMat["el:V_A"].relativeErrorTolerance = 0.01;
    
    // Whenever TORTOISE will be asked to run the time propagation, it will retrieve and use those values
    // for the tolerances
    
    // The time propagation needs an object of type MaterialTimeStatus. This is a container of populations at
    // different times. We need to initialise the first time

    MaterialTimeStatus<2>  evolution(equil);
    // The constructor above add the same population as the one in the equil object in the evolution object
    // Since the initial time is not explicitly specified, it assumes that that is the population at time 0.
    
    // Let us modify the population to make it a non equilibrium one.
    evolution[0]["el:V_A"] -= 0.2*exp(-pow(x-0.12,2.)/0.001-pow(y-0.12,2.)/0.001);
    // We added some holes around a region in the brillouin zone.
    evolution[0]["el:V_A"].plot();
    
    stopWatch.tic(); // We start the timer
    // We can now propagate simply with
    evolution.propagateDeterministic(5., 1.);
    // The first value is the final time and the second value is the time at which we want TORTOISE to
    // add a population within evolution. In other words, the second value is the dense time step.
    // Please notice that the dense time step does not control the precision, it only specifies how
    // frequently the output should be saved.

    cout << "Time taken for time propagation: " << stopWatch.toc() << "s\n";
    
    // We can plot the populations at the final time with
    for (auto elBand : exMat.matchListId("el:*")){
        evolution.back()[elBand].plot(elBand);
    }
    
    for (int i =0; i<evolution.size(); ++i){
        plotter3d << evolution[i]["el:V_A"] << LABEL("time = "+toString(evolution.arg(i), 1));
    }
    plotter3d << PLOT;
    
    return 0;
}
