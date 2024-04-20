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
//  ScatteringChannel.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/11/20.
//

#ifndef ScatteringChannel_hpp
#define ScatteringChannel_hpp

#include <Physics/Material/Material.hpp>
#include <Physics/Status/MaterialStatus.hpp>
#include <Geometry/Structured/DirectSpace/MeshSubset.hpp>

#include <vector>
#include <functional>

namespace Tortoise {

// As ScatteringChannel uses data from the Material class, and the Material class stores and constructs ScatteringChannel,
// it is necessary to provide forward declaration of the other classes
template<int NDim> class Material;
template<int NDim> class MaterialStatus;

template <int NDim, int NLegs> class ScatteringChannel {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main properties
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    const Material<NDim>*                                                   material;
    
    const std::array<int,NLegs>                                             legBandNumber;
    const std::array<const Band<NDim> *, NLegs>                             legBand;        // Pointer to the band corresponding to each leg
    const std::array<scattLegDirection, NLegs>                              legDirection;

    VectorPoint<NDim>                                                       umklappUnitVectors;  // Vector of the umklapp vectors, expressed in terms of the sides of the BZ
                                                                                                 // Example: {{0.0,0.0}, {1.0,0.0}, {0.0,1.0}, {1.0,1.0}, {-1.0,0.0}, ...}
    
    Real                                                                    amplitude = 1.0;
//    std::function <Real(const std::array<Point<NDim>,NLegs>&)>              functAmplitude;   // The matrix element amplitude is given by amplitude*functAmplitude
//    std::function <Real(const std::vector<MeshSubsetElementIterator<NDim>>&)> functAmplitude;
    // amplitude is redundant, however it comes handy since often the user wants to rescale the whole scattering channel without having to modify the often complicated functAmplitude
    std::function<Real(const ArrayPoint<NDim, NLegs>&)>                     functAmplitude;
    // The k points per particle should be passed as columns of the matrix
    
    Real                                                                    minimumDeterminant = 1.e-8;
    Real                                                                    scattElemeThreshold = 1.e-8;

    int                                                                     numMCPoints = 100 ;

    int                                                                     numberNonZeroScatteringElements = -1;
    int                                                                     numberScatteringElementsOverThreshold = -1;

//private:
    
    unsigned long                                                           numberScatteringElements = 0;
    std::vector<std::array<int,NLegs>>                                      scatteringMatrixIndices;
    std::vector<int>                                                        scatteringMatrixIndicesCorrespondingUmklapp;
    std::vector<Eigen::Matrix<Real, NLegs*(NDim+1), Eigen::Dynamic>>        scatteringMatrix;
    const std::vector<std::vector<int>>                                     modalOrdersCombinations;
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main methods
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    //=======================================================
    // Constructors
    //===================
    ScatteringChannel(const Material<NDim>& t_material, const std::array<int, NLegs>& bandsNumbers, const std::array<scattLegDirection, NLegs>& t_legDirection);
    ScatteringChannel(const Material<NDim>& t_material, const std::array<std::string, NLegs>& bandsNames, const std::array<scattLegDirection, NLegs>& t_legDirection);
    ~ScatteringChannel(){};
    
    bool findUmklappVectors(const int nNearestBrillouinZones);
    
    //=======================================================
    // Scattering Analytic tools
    //===================
    // These methods compute a number of characteristcs of the scattering. They do not perform any integration: they simply analyse the list of scattering elements that will be calculated
    
    unsigned long                           getNumberScatteringElements() const {return numberScatteringElements;};
    
    // Deterministic
    unsigned long                           countScatteringElementIDsSingleUmklapp(const Point<NDim>& umklappUnitVect) const;     // Total number of scattering elements (can be a very large number)
    unsigned long                           countScatteringElementIDs();                  // Recounts total number of scattering elements (can be a very large number) and stores in numberScatteringElements
    std::vector<std::array<int,NLegs>>      scatteringElementIDs() const;                       // Vector of element indices of all the scattering elements (can be a very very large vector)
    std::array<std::vector<int>,NLegs>      scatteringPerElement() const;                       // Vector of number of scattering elements per mesh element (very rough estimation of joined density of state)
    Function<NDim>                          scatteringPerElementPlot(const int which) const;    // Same as above, but returns a Function<NDim> (for easy plotting and visualisation)
        
    //=======================================================
    // Scattering Integrators setup
    //===================
    
    //=======================================================
    // Scattering Precalculator
    //===================
    void precalculate();
    void precalculateIndices();
    
    void outputScatteringFile();
    
    //=======================================================
    // Scattering Propagators
    //===================
    MaterialStatus<NDim> propagate(const MaterialStatus<NDim>& status) const;
    
    //=======================================================
    // Scattering Rates
    //===================
    // This returns the scattering rate for the given particle/band, summing over all the legs where it appears
    Function<NDim> scatteringRates(const MaterialStatus<NDim>& status, const std::string &bandName) const;
    
    
//private:
    //=======================================================
    // Technicalities (not intended for final user)
    //===================
    void printmodalOrdersCombinations();
    std::string scatteringName() const;
    
    //=======================================================
    // Iterators
    //===================
    // The following methods require a rather technical understanding of how the scattering object works.
    // These are not intended for the average user
    template <typename functType> void deterministicIteration(const Point<NDim>& umklappUnitVect, const functType& functionToExecute) const;
    // functionToExecute has to take std::vector<MeshSubsetElementIterator<NDim>>& iter as argument and return void
    
    template <typename functType> void deterministicIterationParallel(const Point<NDim>& umklappUnitVect, const functType& functionToExecute) const;
    
    //=======================================================
    // Templated nested loops over legs
    //===================
    template<int depth, typename functType> void nestedForElements(std::vector<MeshSubsetElementIterator<NDim>>& iter, std::vector<MeshSubset<NDim>>& meshSubset, std::array<Real,NLegs>& minE, std::array<Real,NLegs>& maxE, const functType& functionToExecute) const;
    template<int depth, typename functType> void nestedForSections(std::vector<SectionIterator<NDim>>& iter, std::vector<MeshSubset<NDim>>& meshSubset, const Point<NDim>& umklappUnitVect, const functType& functionToExecute) const;
    
    //=======================================================
    // Calculation over a single scattering element
    //===================
    void singleScatteringRate(const Point<NDim>& umklappUnitVect, const MaterialStatus<NDim>& status, MaterialStatus<NDim>& rates, const int which, const std::vector<MeshSubsetElementIterator<NDim>>& elemID) const;
    void singleScatteringRate(const Point<NDim>& umklappUnitVect, const MaterialStatus<NDim>& status, MaterialStatus<NDim>& rates, const int which, const std::array<int,NLegs>& elemID) const;
    void singleScatteringPropagator(const Point<NDim>& umklappUnitVect, const MaterialStatus<NDim>& status, MaterialStatus<NDim>& propagation, const std::vector<MeshSubsetElementIterator<NDim>>& elemID) const;
    void singleScatteringPropagator(const Point<NDim>& umklappUnitVect, const MaterialStatus<NDim>& status, MaterialStatus<NDim>& propagation, const std::array<int,NLegs>& elemID) const;
    MaterialStatus<NDim> propagatePrecalculatedSingle(const MaterialStatus<NDim>& status, const unsigned long scatElem) const;
    
    //=======================================================
    // Non symmetrised scatterings (intermediate step in calculating the physically meaningful scattering)
    //===================
    // This returns the scattering rate only for the selected leg (this is used for technical reasons and does not have a physical meaning,
    // since it ignores particle exchange symmetry)
    Function<NDim> scatteringRates(const MaterialStatus<NDim>& status, const int which) const;
    
    
// Technicalities
    
    void constructMCPoints();
    VectorPoint<std::max(NDim*(NLegs-1)-1,2*NDim)>  mCPoints;
    // Note to myself: I made sure that the dimension is never smaller than 2*NDim, to avoid assert errors
    // Fix later
};




//=======================================================
// Implementation Details
//===================

//=======================================================
// Iterators
//===================

template<int NDim, int NLegs> template<typename functType> void ScatteringChannel<NDim,NLegs>::deterministicIteration(const Point<NDim>& umklappUnitVect, const functType& functionToExecute) const {
    std::vector<SectionIterator<NDim>> iter;
    for (int i=0; i < NLegs-1; ++i){ iter.push_back(legBand[i]->mesh->sectionIterator());}
    std::vector<MeshSubset<NDim>> meshSubset;
    for (int i=0; i < NLegs; ++i){ meshSubset.emplace_back(legBand[i]->mesh);}
    nestedForSections<NLegs-1>(iter, meshSubset, umklappUnitVect, functionToExecute);
};

//=======================================================
// Templated nested loops over legs
//===================
template<int NDim, int NLegs> template<int depth, typename functType> void ScatteringChannel<NDim,NLegs>::nestedForElements
(std::vector<MeshSubsetElementIterator<NDim>>& iter, std::vector<MeshSubset<NDim>>& meshSubset, std::array<Real,NLegs>& minE, std::array<Real,NLegs>& maxE, const functType& functionToExecute) const {
    if constexpr (depth>0){
        for (iter[NLegs-depth] = meshSubset[NLegs-depth].elementIterator(); iter[NLegs-depth].unfinished; ++iter[NLegs-depth]){
            minE[NLegs-depth] = (legDirection[NLegs-depth]>0)? legBand[NLegs-depth]->minElemEnergies(iter[NLegs-depth].currentElementID): -legBand[NLegs-depth]->maxElemEnergies(iter[NLegs-depth].currentElementID);
            maxE[NLegs-depth] = (legDirection[NLegs-depth]>0)? legBand[NLegs-depth]->maxElemEnergies(iter[NLegs-depth].currentElementID): -legBand[NLegs-depth]->minElemEnergies(iter[NLegs-depth].currentElementID);
            nestedForElements<depth-1>(iter, meshSubset, minE, maxE, functionToExecute);
        }
    } else {
        Real minEtot = 0., maxEtot = 0.;
        for (int i=0; i<NLegs; ++i) { minEtot += minE[i]; maxEtot += maxE[i];}
        if ( (minEtot<=0.0) && (maxEtot>=0.0)) { functionToExecute(iter);}
    }
}

template<int NDim, int NLegs> template<int depth, typename functType> void ScatteringChannel<NDim,NLegs>::nestedForSections
(std::vector<SectionIterator<NDim>>& iter, std::vector<MeshSubset<NDim>>& meshSubset, const Point<NDim>& umklappUnitVect, const functType& functionToExecute) const {
    
    if constexpr (depth>0){
        for (iter[NLegs-1-depth] = legBand[NLegs-1-depth]->mesh->sectionIterator(); iter[NLegs-1-depth].unfinished; ++iter[NLegs-1-depth]){
            meshSubset[NLegs-1-depth] = MeshSubset<NDim>(legBand[NLegs-1-depth]->mesh, iter[NLegs-1-depth]);
            nestedForSections<depth-1>(iter, meshSubset, umklappUnitVect, functionToExecute);
        }
    } else {
        if constexpr (NLegs == 2) {
            meshSubset[NLegs-1] = -static_cast<Real>(legDirection[NLegs-1])*(static_cast<Real>(legDirection[0]) * meshSubset[0]+ material->region->gVec * umklappUnitVect);
        }
        if constexpr (NLegs == 3) {
            meshSubset[NLegs-1] = -static_cast<Real>(legDirection[NLegs-1])*(static_cast<Real>(legDirection[0]) * meshSubset[0]+ static_cast<Real>(legDirection[1]) * meshSubset[1] + material->region->gVec * umklappUnitVect);
        }
        if constexpr (NLegs == 4) {
            meshSubset[NLegs-1] = -static_cast<Real>(legDirection[NLegs-1])*(static_cast<Real>(legDirection[0]) * meshSubset[0]+ static_cast<Real>(legDirection[1]) * meshSubset[1]+ static_cast<Real>(legDirection[2])*meshSubset[2] + material->region->gVec * umklappUnitVect);
        }
        std::vector<MeshSubsetElementIterator<NDim>> iterElem;
        for (int i=0; i < NLegs; ++i){ iterElem.emplace_back(meshSubset[i].elementIterator());}
        std::array<Real,NLegs> minE;
        std::array<Real,NLegs> maxE;
        nestedForElements<NLegs>(iterElem, meshSubset, minE, maxE, functionToExecute);
    }
}



} // namespace Tortoise



#endif /* ScatteringChannel_hpp */
