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
//  main.cpp
//  Tutorial3
//
//  Created by Marco Battiato on 31/10/22.
//
// TUTORIAL 3: Scattering rates


// Please remember that the graphical output will work only if you have written
// the GNUPlot running directory in the Configuration file

#include <TORTOISE>
#include <iostream>

using namespace Tortoise;
using std::cout;

int main(int argc, const char * argv[]) {
    
    // Let us build the material
    real                pi = 3.1415;
    Point<2>            originBZ({0.0, 0.}), sideBZ0({1.0, 0.0}), sideBZ1({-.5, std::sqrt(3.0)/2});
    Point<2>            ortsideBZ0({.0, 1.0}), ortsideBZ1({std::sqrt(3.0)/2, .5});
    Region<2>           brillouinZone (originBZ,{sideBZ0,sideBZ1});
    Material<2>         exMat ("ExampleMaterial", brillouinZone);
    
    int                 resolution = 13;
    Mesh<2>             meshFull(brillouinZone, {0.,0.},{1.,1.},{resolution,resolution});
    Mesh<2>             meshA(brillouinZone, {0.1,0.1},{0.3,0.3},{resolution,resolution});
    Mesh<2>             meshB(brillouinZone, {0.6,0.6},{0.3,0.3},{resolution,resolution});
    Mesh<2>             meshPhot(brillouinZone, {0.,0.},{0.001,0.0001},{5*resolution,1});

    real                photMaxEnergy = 2.5;

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
    
 
    // We now construct the population of the material.
    // Population can be equilibrium or out of equilibrium. Let us start with the population at equilibrium

    // We first create an object that can contain the the population
    MaterialStatus<2>   equil(exMat);
    // Notice how the constructor needs the material, as it needs to know the bands for which it needs to accommodate populations
    
    const real chemPot = 0.0, temp = 0.1;
    // We now assign the populations. The population of each band in the MaterialStatus object is accessed using the same string
    // id as the original band. For instance we can
    equil["el:V_A"] = 1/(exp((exMat["el:V_A"].dispersion-chemPot)/temp)+1);
    
    // The above method however can be rather tedious for many bands. We can iterate with
    for (auto band: exMat.matchListId("*")){
        if(exMat[band].statistics == fermion)
            equil[band] = 1/(exp((exMat[band].dispersion-chemPot)/temp)+1);
        else
            equil[band] = 1/(exp((exMat[band].dispersion-chemPot)/temp)-1);
    }
    // For later convenience let us empty the photon band, since we will use it to contain photons from a fs laser
    equil["pt:Lin"] = 0.;
    
    equil["el:V_A"].plot("Equil pop in el:V_A");

    
    // We can now calculate the scattering lifetimes of scatterings. In reality TORTOISE calculates the scattering rates (the inverse
    // of the scattering lifetimes)
    // This is done simply by
    auto phoScattRate = exMat.scatteringRates<3>({"el:V_A", "pn:ac0", "el:V_A"}, {in, in, out}, equil, "el:V_A");
    // This requires us to specify the scattering we want to addres by the names of the involved bands and their directions.
    // Notice that if the scattering does not exist the scattering rate will be a constant 0.
    // To calculate the scatterin rate we also need to provide the population over which to calculate it.
    // Usually scattering rates are calcualted at equilibrium, yet we will see that there are other interesting applications of
    // scatterin rates
    // Finally we need to specify the band we want to calculate the scattering rates for. If the mentioned band appears more than
    // once in the list of bands involved in the scattering, all contributions will be included.
    
    // We could obtain the lifetimes by inverting, but we have to be careful as the inversion does not preserve numerical precision
    // Even more importantly, TORTOISE's plotter does not allow yet for plotting over a certain range. This means that when thes cattering
    // rates become very small the lifetimes become very large, making the figure for lifetimes often difficult to interpret (keep
    // in mind that the interesting regions are the ones with the smallest scattering lifetimes, which are often hard to see in a plot where
    // we cannot control the plot range)
    (1./phoScattRate).plot("el:V_A + pn:ac0 <-> el:V_A lifetimes for el:V_A");
    // It is better to visualise the scattering rates themselves
    phoScattRate.plot("el:V_A + pn:ac0 <-> el:V_A scatt rates for el:V_A");
    
    // Often however scattering rates of single scattering channels are not directly comparable with experiments. It is often needed to
    // sum over a number of scattering channels. this can be obtained with the same method, using wildcards

    exMat.scatteringRates<3>({"el:*", "pn:*", "el:*"}, {in, in, out}, equil, "el:V_A").plot("electron-phonon scatt rates for band el:V_A");
    
    
    
    // Interestingly the absorption spectrum of a material can be constructed as a scattering rate
    auto absSpectK = exMat.scatteringRates<3>({"el:*", "pt:*", "el:*"}, {in, in, out}, equil, "pt:Lin");

    // The above expression returns a photon momentum resolved spectrum. That is not very comfortable to look at, as we would like to know the
    // photon energy resolved absorption spectrum.
    // Since we want to use TORTOISE's plotting routines, we would like to get a Function<1> in the end.
    // The domain of that function should be the range of possible photon energies
    Region<1>    energyRange (exMat["pt:Lin"].dispersion.min(), exMat["pt:Lin"].dispersion.max()-exMat["pt:Lin"].dispersion.min());
    Mesh<1>      energyMesh  (energyRange, 0., 1., 5*resolution);
    // We now use the interpolator defined within TORTOISE
    Interpolator functInterp;
    for (auto elem = exMat["pt:Lin"].dispersion.elementIterator(); elem.unfinished; ++elem) {
        functInterp.addData( exMat["pt:Lin"].dispersion(elem, exMat["pt:Lin"].dispersion.mesh->elemCentre(elem)), absSpectK(elem, exMat["pt:Lin"].mesh->elemCentre(elem)) );
        for (auto vert = exMat["pt:Lin"].dispersion.vertexInElemIterator(); vert.unfinished; ++vert){
            Point<2> point = 0.98 * (*exMat["pt:Lin"].mesh)(elem, vert) + 0.02 * exMat["pt:Lin"].dispersion.mesh->elemCentre(elem);
            functInterp.addData( exMat["pt:Lin"].dispersion(elem, point), absSpectK(elem, point) );
            point = 0.5 * (*exMat["pt:Lin"].mesh)(elem, vert) + 0.5 * exMat["pt:Lin"].mesh->elemCentre(elem);
            functInterp.addData( exMat["pt:Lin"].dispersion(elem, point), absSpectK(elem, point) );
        }
    }
    functInterp.prepareForInterpolation();
    auto absSpect = Function<1>(energyMesh, functInterp);
    
    absSpect.plot();

    return 0;
}
