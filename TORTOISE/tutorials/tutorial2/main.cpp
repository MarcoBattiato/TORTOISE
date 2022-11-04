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
//  Tutorial2
//
//  Created by Marco Battiato on 31/10/22.
//
// TUTORIAL 2: Constructing a material


// Please remember that the graphical output will work only if you have written
// the GNUPlot running directory in the Configuration file

#include <TORTOISE>
#include <iostream>

using namespace Tortoise;
using std::cout;

int main(int argc, const char * argv[]) {
    
    // We will build a very simple paterial with bands with analytical expressions
    
    // We first must build the brillouin zone.
    
    real                pi = 3.1415;
    Point<2>            originBZ({0.0, 0.}),
                        sideBZ0({1.0, 0.0}),
                        sideBZ1({-.5, std::sqrt(3.0)/2});
    Region<2>           brillouinZone (originBZ,{sideBZ0,sideBZ1});
    
    // Notice how what we (improperly) call the Brillouin zone is actually the primiteve cell of the
    // reciprocal lattice.
    // TORTOISE can only handle regions that are lines in 1D, parallelograms in 2D and parallelepypeds in 3D
    
    // We build the material by passing a name, and its BZ
    Material<2>         exMat ("ExampleMaterial", brillouinZone);
    
    // The material can be thought as an, at the momen empty, container of bands
    
    // We now add bands
    // Each band needs to be defined over a mesh.
    // VERY IMPORTANTLY the used meshes must live through the whole program since the bands will use the mesh
    // object they were created on. Do not create temporary Mesh objects, or Mesh objects with limited scope!!!
    
    
    int                 resolution = 10;
    // Let us first add a phonon band over the full BZ. We need to construct the mesh
    Mesh<2>             meshFull(brillouinZone, {0.,0.},{1.,1.},{resolution,resolution});
    
    // The constructor of a Mesh object requires:
    // - region
    // - origin of the mesh
    // - size of the mesh
    // - resolution per side
    //
    // It is important to know that the origin must be specified as a fraction of the region.
    // This means that {0,0} means that the origin of the mesh is the same as the one of the region
    // {0.5,0.5} means that the origin on the mesh is at the centre of the region
    // {1,1} means that the origin on the mesh is at the top right corner of the region
    //
    // Furthermore the size of the mesh is given in fractions of the region sizes.
    // {0.1,0.5} means that the mesh has the first side which is 10% of the first side of the mesh
    // and the second side which is half that of the mesh.
    // Notice that the mesh does not need to stay within the region.
    // However the user is strongly discouraged to use meshes that are larger than the region itself.
    // That will impact how umklapp is calculated during the evaluation of the scattering.
    
    // We can now add the phonon band. For simplicity we assume it to be a very simplified optical branch
    // with constant energy
    
    exMat.addBand("pn:ac0", 0, boson, meshFull, 0.05);
    
    // Notice how we pass to the contructor:
    // - name
    // - charge
    // - statistics
    // - mesh
    // - dispersion
    //
    // The name is extremely important since it will be the primary way to address the bands.
    // Since TORTOISE provides way to iterate through bands with wildcards, the user should use names
    // that allow for an efficient use of wildcards.
    // The statistic can be fermion or boson (notice that these names are not strings)

    // Let us now add electronic bands
    // A material does not group different quasiparticles. They all go into the same container of bands.
    // We build two new meshes for the electronic bands
    
    Mesh<2>             meshA(brillouinZone, {0.1,0.1},{0.3,0.3},{resolution,resolution});
    Mesh<2>             meshB(brillouinZone, {0.6,0.6},{0.3,0.3},{resolution,resolution});
    
    // and add four electronic bands
    Point<2>            ortsideBZ0({.0, 1.0}), ortsideBZ1({std::sqrt(3.0)/2, .5});
    
    exMat.addBand("el:V_A", -1, fermion, meshA, -1.1+sin(pi*dot(ortsideBZ0,k - 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2))*sin(pi*dot(ortsideBZ1,k- 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2)));
    exMat.addBand("el:C_A", -1, fermion, meshA, 1.1-sin(pi*dot(ortsideBZ0,k - 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2))*sin(pi*dot(ortsideBZ1,k- 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2)));
    exMat.addBand("el:V_B", -1, fermion, meshB, -1.1+sin(pi*dot(ortsideBZ0,k - 0.601*sideBZ0 - 0.6*sideBZ1)/(0.3*std::sqrt(3.0)/2))*sin(pi*dot(ortsideBZ1,k- 0.601*sideBZ0 - 0.6*sideBZ1)/(0.3*std::sqrt(3.0)/2)));
    exMat.addBand("el:C_B", -1, fermion, meshB, 1.1-sin(pi*dot(ortsideBZ0,k - 0.601*sideBZ0 - 0.6*sideBZ1)/(0.3*std::sqrt(3.0)/2))*sin(pi*dot(ortsideBZ1,k- 0.601*sideBZ0 - 0.6*sideBZ1)/(0.3*std::sqrt(3.0)/2)));

    // Notice how we used the name to contain the name of the quasiparticle
    // We construct the bands by passing their analytic dispersion.
    
    // We finally add a photonic band. This is because we can treat the laser excitation as a scattering between a photon
    // and quasiparticles within the crystal.
    // To save some computational cost, we will treat optical photons in a simplified way.
    // First of all, we realise that while photons have momentum, momentum at optical energies is so small that can be completely neglected.
    // What will really matter is the energy. From this point of view we will build a photon band that does not accurately describe the
    // momentum of the phonon, but describes well its energy.
    // We notice that the direction of the photon (not its polarisation) is quite irrelevant, given the very small size of the momentum.
    // Therefore we might as well describe only photons in any given direction. Furthermore we want to create a photon's band domain small
    // enough in k space to have negligible impact on the momentum conservation, but not too small to create problems with numerical precision
    // and rounding during plotting
    // For that reason we do not really use the realistic maximum momentum that a photon at the energies we are interested in will have (as it
    // might be so small to create problems with rounding at machine precision.
    
    Mesh<2>             meshPhot(brillouinZone, {0.,0.},{0.001,0.0001},{5*resolution,1});

    // Notice that we have created a mesh that is small (but not too small). It is 10 times largen in the first dimesion compared to the second.
    // We have then assigned the desired resolution on the longest direction (that will reflect in the obtained energy resolution) and only one in
    // the other direction.
    
    real                photMaxEnergy = 2.5;
    exMat.addBand("pt:Lin", 0, boson, meshPhot, photMaxEnergy*dot(ortsideBZ1, k)/(0.001*std::sqrt(3.0)/2));
    
    // We define the photon dispersion as a plane. Notice that when constructing the dispersion, we need to be careful since using a more straightforward
    // dispersion (2.5*kx) would create photons with negative energy (plot it to see what I mean).
    
    // We will see later how to add bands with dispersions taken from files or constructed via other codes.
    
    
    // We can now plot the material
    
    exMat.plot();
    // Notice how the photon band is visible
    // Or plot only specific bands
    exMat.plot({"pt:Lin"});
    exMat.plot(exMat.matchListId("el:*"));
    
    // Each band within the material is accessed by its name
    // For instance exMat["el:V_A"] is an object of type band.
    // A Band object contains information about the band itself.
    // The most important characteristic stored is the the dispersion of the band. The dispersion is a Function object
    
    exMat["el:V_A"].dispersion.plot("el:V_A");
    
    
    // We can add properties to the bands. For instance we can store the spin polarisation
    
    exMat["el:V_A"].attribute.emplace_back("spin pol", +1);
    
    // which can be retrieved with
    
    cout << exMat["el:V_A"].attribute["spin pol"] << "\n";
    
    // Notice that user defined attributes are not used by TORTOISE automatically. They can be used by the user when, for instance, constructing scatterings, or when calculating expectation values.
    
    // We can also add attributes that are functions
    
    exMat["el:V_A"].functionAttribute.emplace_back("spin pol", Function<2>(exMat["el:V_A"].mesh, x*x+y));

    // Notice that the analytic expression will be discretised over the same mesh as the band's before being stored.
    // Similarly to the bands, these attributes can be built with, for instance, calls to ab initio codes.
    
    exMat["el:V_A"].functionAttribute["spin pol"].plot();
    
    
    
    // To complete the description of a material, we need to add scattering channels
    // Scattering channels can be added one by one, but that can quickly ebcome very cumbersome, as in a realistic material there
    // can be even hundreds of them. It is better to use provided functionalities.
    
    // We add all electron phonon scatterings with
    for (auto scatt : exMat.addScatteringChannelsWildcards<3>({"el:*", "pn:*", "el:*"}, {in, in, out}, 2)){
        scatt->numMCPoints = 50;            // Number of Monte Carlo Points used in the integration of each scattering element
        scatt->minimumDeterminant = 1.e-6;  // The inversion of the Diract delta required the solution of a linear system. Too
                                            // small determinants can lead to unstable results
        scatt->amplitude = 1;               // This is the scattering matrix element assumed constant (we will see later how to
                                            // do better)
    }
    
    // The method
    // exMat.addScatteringChannelsWildcards<3>({"el:*", "pn:*", "el:*"}, {in, in, out}, 2)
    // creates all the combinations of the abnds specified by the given wildcards and produces all the scattering channels.
    // The method then searches for scattering channels that are duplicates by exchange or time reversal symmetry, and does
    // not construct them.
    // The last integer in the method call represents the desired umklapp nearest neighbours that the method should check.
    // addScatteringChannelsWildcards explores the scattering channel for scattering elements up to the given umklapp nearest
    // neighbours and remembers which umklapp vectors can indeed lead to scattering events. If the method finds out that there
    // is no way for tehscattering channel to satisfy energy and momentum conservation, it will completely remove the scattering
    // channel.
    
    // After all the work is done, addScatteringChannelsWildcards returns a vector of pointers to the scattering channles that have
    // just been created. This is useful, since one can loop through that to assign some characteristics to the scattering channels.
    
    
    // We add now the absoprtion processes
    for (auto scatt : exMat.addScatteringChannelsWildcards<3>({"el:*", "pt:*", "el:*"}, {in, in, out}, 1)) {
        scatt->numMCPoints = 50;
        scatt->minimumDeterminant = 1.e-6;
        scatt->amplitude = 1;
    }
    
    // We now print the material, to see some data
    
    cout << exMat;
    
    return 0;
}
