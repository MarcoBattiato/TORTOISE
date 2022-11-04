//
//  main.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 27/10/22.
//
// Debug mode is activated by default
// #define NDEBUG            // Uncomment to deactivate debug mode

#include <TORTOISE>
#include <iostream>

using namespace Tortoise;
using std::cout;

template <int NDim> Function<1> functVsFunct(const Function<NDim>& funcx, const Function<NDim>& funcy, const Mesh<1> mesh){
    assert(funcx.mesh == funcy.mesh);
    
    Interpolator functInterp;
    for (auto elem = funcx.elementIterator(); elem.unfinished; ++elem) {
        functInterp.addData( funcx(elem, funcx.mesh->elemCentre(elem)), funcy(elem, funcx.mesh->elemCentre(elem)) );
        for (auto vert = funcx.vertexInElemIterator(); vert.unfinished; ++vert){
            Point<NDim> point = 0.98 * (*funcx.mesh)(elem, vert) + 0.02 * funcx.mesh->elemCentre(elem);
            functInterp.addData( funcx(elem, point), funcy(elem, point) );
            point = 0.5 * (*funcx.mesh)(elem, vert) + 0.5 * funcx.mesh->elemCentre(elem);
            functInterp.addData( funcx(elem, point), funcy(elem, point) );
        }
    }
    functInterp.prepareForInterpolation();
    return Function<1>(mesh, functInterp);
}

int main(int argc, const char * argv[]) {
    
    Real                pi = 3.1415;
    
    Point<2>            originBZ({0.0, 0.}), sideBZ0({1.0, 0.0}), sideBZ1({-.5, std::sqrt(3.0)/2});
    Point<2>            ortsideBZ0({.0, 1.0}), ortsideBZ1({std::sqrt(3.0)/2, .5});
    Region<2>           brillouinZone (originBZ,{sideBZ0,sideBZ1});
    Material<2>         exMat ("ExampleMaterial", brillouinZone);
    
    int                 res = 10;
    Mesh<2>             meshFull(brillouinZone, {0.,0.},{1.,1.},{res,res});
    Mesh<2>             meshA(brillouinZone, {0.1,0.1},{0.3,0.3},{res,res});
    Mesh<2>             meshB(brillouinZone, {0.6,0.6},{0.3,0.3},{res,res});
    Mesh<2>             meshPhot(brillouinZone, {0.,0.},{0.001,0.0001},{5*res,1});
    
    // *** Add Bands
    // Electrons
    exMat.addBand("el:V_A", -1, fermion, meshA, -1.1+sin(pi*dot(ortsideBZ0,k - 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2))*sin(pi*dot(ortsideBZ1,k- 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2)));
    exMat.addBand("el:C_A", -1, fermion, meshA, 1.1-sin(pi*dot(ortsideBZ0,k - 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2))*sin(pi*dot(ortsideBZ1,k- 0.101*sideBZ0 - 0.1*sideBZ1)/(0.3*std::sqrt(3.0)/2)));
    exMat.addBand("el:V_B", -1, fermion, meshB, -1.1+sin(pi*dot(ortsideBZ0,k - 0.601*sideBZ0 - 0.6*sideBZ1)/(0.3*std::sqrt(3.0)/2))*sin(pi*dot(ortsideBZ1,k- 0.601*sideBZ0 - 0.6*sideBZ1)/(0.3*std::sqrt(3.0)/2)));
    exMat.addBand("el:C_B", -1, fermion, meshB, 1.1-sin(pi*dot(ortsideBZ0,k - 0.601*sideBZ0 - 0.6*sideBZ1)/(0.3*std::sqrt(3.0)/2))*sin(pi*dot(ortsideBZ1,k- 0.601*sideBZ0 - 0.6*sideBZ1)/(0.3*std::sqrt(3.0)/2)));
    // Phonons
    exMat.addBand("pn:ac0", 0, boson, meshFull, 0.05);
    // Photons
    exMat.addBand("pt:Lin", 0, boson, meshPhot, 2.5*dot(ortsideBZ1,k)/(0.001*std::sqrt(3.0)/2));
    
    exMat.plot();
   
    // *** Add Scattering Channels
    for (auto scatt : exMat.addScatteringChannelsWildcards<3>({"el:*", "pn:*", "el:*"}, {in, in, out}, 2)){
        scatt->numMCPoints = 50;  scatt->minimumDeterminant = 1.e-6; scatt->amplitude = 1;
    }
    for (auto scatt : exMat.addScatteringChannelsWildcards<3>({"el:*", "pt:*", "el:*"}, {in, in, out}, 1)) {
        scatt->numMCPoints = 50;   scatt->minimumDeterminant = 1.e-6;
    }
    
    cout << exMat;

    // *** Create equilibrium population
    MaterialStatus<2>   equil(exMat);
    const Real chemPot = 0.0, temp = 0.1;
    for (auto band: exMat.matchListId("*")){
        if(exMat[band].statistics == fermion)
            equil[band] = 1/(exp((exMat[band].dispersion-chemPot)/temp)+1);
        else
            equil[band] = 1/(exp((exMat[band].dispersion-chemPot)/temp)-1);
    }
    
    equil["el:V_A"].plot();
    
    // *** Eletron-phonon scattering rate
    exMat.scatteringRates<3>({"el:*", "pn:*", "el:*"}, {in, in, out}, equil, "el:V_A").plot();
    
    // *** Absoprtion spectrum
    auto absSpectK = exMat.scatteringRates<3>({"el:*", "pt:*", "el:*"}, {in, in, out}, equil, "pt:Lin");
    Region<1>    energyRange (exMat["pt:Lin"].dispersion.min(),exMat["pt:Lin"].dispersion.max()-exMat["pt:Lin"].dispersion.min());
    Mesh<1>      energyMesh  (energyRange, 0., 1., 5*res);
    auto absSpect =  functVsFunct(exMat["pt:Lin"].dispersion, absSpectK, energyMesh);
    absSpect.plot();
    
    exMat["pt:Lin"].constrained = true;
    exMat["pt:Lin"].constrainedPopulation = [&exMat](Real t){return 1.e7*std::exp(-std::pow(t-0.3,2.)/0.1)*exp(-pow(exMat["pt:Lin"].dispersion-1.4,2)/0.04);};
    
    // *** Time evolution with laser excitation
    MaterialTimeStatus<2> thermDyn(equil);
    stopWatch.tic();
    thermDyn.propagateDeterministic(5, 1.);
    cout << "Time taken for time propagation: " << stopWatch.toc() << "s\n";
    
    for (auto elBand : exMat.matchListId("el:*")){
        thermDyn.back()[elBand].plot(elBand);
    }
    
    return 0;
}

