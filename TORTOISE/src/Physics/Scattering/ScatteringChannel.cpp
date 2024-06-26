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
//  ScatteringChannel.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/11/20.
//

#include <Physics/Scattering/ScatteringChannel.hpp>
#include <Physics/Scattering/IntegratorCore.hpp>
#include <Geometry/Structured/DirectSpace/MeshSubset.hpp>
#include <Generics/Utilities/RandomNGenerator.hpp>

#include <deque>
#include <fstream>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

using namespace Tortoise::GeometryCore;
using namespace Tortoise::PhysicsCore;

namespace Tortoise {

//=======================================================
// Constructors
//===================
template<int NDim, int Nlegs> std::vector<std::vector<int>> produceModalOrdersCombinations(){
    std::vector<std::vector<int>> orderSTDVEC;
    orderSTDVEC.emplace_back(Nlegs,0);
    for (int i=0; i<Nlegs; ++i){
        for (int l=1; l<NDim+1; ++l){
            std::vector<int> temp(Nlegs,0);
            temp[i]= l;
            orderSTDVEC.emplace_back(temp);
        }
    }
    for (int i=0; i<Nlegs; ++i){
        for (int j=i+1; j<Nlegs; ++j){
            for (int l=1; l<NDim+1; ++l){
                for (int m=1; m<NDim+1; ++m){
                    std::vector<int> temp(Nlegs,0);
                    temp[i]= l; temp[j]= m;
                    orderSTDVEC.emplace_back(temp);
                }
            }
        }
    }
//    for (int i=0; i<Nlegs; ++i){
//        for (int j=i+1; j<Nlegs; ++j){
//            for (int i1=j+1; i1<Nlegs; ++i1){
//                for (int l=1; l<NDim+1; ++l){
//                    for (int m=1; m<NDim+1; ++m){
//                        for (int l1=1; l1<NDim+1; ++l1){
//                            std::vector<int> temp(Nlegs,0);
//                            temp[i]= l; temp[j]= m; temp[i1]= l1;
//                            orderSTDVEC.emplace_back(temp);
//                        }
//                    }
//                }
//            }
//        }
//    }
//    for (int i=0; i<Nlegs; ++i){
//        for (int j=i+1; j<Nlegs; ++j){
//            for (int i1=j+1; i1<Nlegs; ++i1){
//                for (int j1=i1+1; j1<Nlegs; ++j1){
//                    for (int l=1; l<NDim+1; ++l){
//                        for (int m=1; m<NDim+1; ++m){
//                            for (int l1=1; l1<NDim+1; ++l1){
//                                for (int m1=1; m1<NDim+1; ++m1){
//                                    std::vector<int> temp(Nlegs,0);
//                                    temp[i]= l; temp[j]= m; temp[i1]= l1; temp[j1]= m1;
//                                    orderSTDVEC.emplace_back(temp);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
    return orderSTDVEC;
}
 
// Produces the pointers to the bands
template<int NDim, int Nlegs> std::array<const Band<NDim> *, Nlegs> assignlegBand (const Material<NDim>& t_material, std::array<int, Nlegs> bandsNumbers) {
    std::array<const Band<NDim> *, Nlegs>   legBand;
    for (int i=0; i<Nlegs; ++i){
        assert(bandsNumbers[i]>=0 && bandsNumbers[i]<t_material.size());
        legBand[i] = &(t_material[bandsNumbers[i]]);
    }
    return legBand;
};
template<int NDim, int Nlegs> std::array<const Band<NDim> *, Nlegs> assignlegBand (const Material<NDim>& t_material, const std::array<std::string, Nlegs>& bandsNames) {
    std::array<const Band<NDim> *, Nlegs>   legBand;
    for (int i=0; i<Nlegs; ++i){
        legBand[i] = &(t_material[bandsNames[i]]);
    }
    return legBand;
};
template<int NDim, int Nlegs> std::array<int, Nlegs> assignlegBandNumber (const Material<NDim>& t_material, const std::array<std::string, Nlegs>& bandsNames) {
    std::array<int, Nlegs>   legBand;
    for (int i=0; i<Nlegs; ++i){
        legBand[i] = t_material.id(bandsNames[i]);
    }
    return legBand;
};
template<int NDim, int Nlegs> VectorPoint<NDim> standardUmklappVectors (const std::array<scattLegDirection, Nlegs>& t_legDirection) {
    VectorPoint<NDim> toreturn;
    toreturn.resize(Eigen::NoChange,1);
    toreturn.col(0) = Point<NDim>::Zero();
    return toreturn;
};

template<int NDim, int Nlegs> ScatteringChannel<NDim,Nlegs>::ScatteringChannel(const Material<NDim>& t_material, const std::array<int, Nlegs>& bandsNumbers, const std::array<scattLegDirection, Nlegs>& t_legDirection): material(&t_material), legBand(assignlegBand<NDim,Nlegs>(t_material,bandsNumbers)), legDirection(t_legDirection), legBandNumber(bandsNumbers), modalOrdersCombinations(produceModalOrdersCombinations<NDim,Nlegs>()),
functAmplitude([](const ArrayPoint<NDim, Nlegs>&)->Real{return 1.;}),
//functAmplitude([](const std::array<Point<NDim>,Nlegs>&)->Real{return 1.;}),
//functAmplitude([](const std::vector<MeshSubsetElementIterator<NDim>>&)->Real{return 1.;}),
umklappUnitVectors(standardUmklappVectors<NDim,Nlegs>(t_legDirection)) {};

template<int NDim, int Nlegs> ScatteringChannel<NDim,Nlegs>::ScatteringChannel(const Material<NDim>& t_material, const std::array<std::string, Nlegs>& bandsNames, const std::array<scattLegDirection, Nlegs>& t_legDirection) : material(&t_material), legBand(assignlegBand<NDim,Nlegs>(t_material,bandsNames)), legDirection(t_legDirection), legBandNumber(assignlegBandNumber<NDim,Nlegs>(t_material,bandsNames)), modalOrdersCombinations(produceModalOrdersCombinations<NDim,Nlegs>()),
functAmplitude([](const ArrayPoint<NDim, Nlegs>&)->Real{return 1.;}),
//functAmplitude([](const std::array<Point<NDim>,Nlegs>&)->Real{return 1.;}),
//functAmplitude([](const std::vector<MeshSubsetElementIterator<NDim>>&)->Real{return 1.;}),
umklappUnitVectors(standardUmklappVectors<NDim,Nlegs>(t_legDirection)) {};
 
template<int NDim> VectorPoint<NDim> buildTestUmklappVects(const int nNearestBrillouinZones, const Material<NDim>& material); // Generates all umklapp vectors (including no umklapp) up to the chosen number of nearest brillouin zones

template<> VectorPoint<1> buildTestUmklappVects<1>(const int nNearestBrillouinZones, const Material<1>& material){
    VectorPoint<1> toreturn;
    const int numberTestUmklappVects = Utilities::intPower<1>(2*nNearestBrillouinZones+1);
    toreturn.resize(Eigen::NoChange, numberTestUmklappVects);
    int counter = 0;
    for (int i= - nNearestBrillouinZones; i <= nNearestBrillouinZones; ++i){
        toreturn.col(counter++) = Point<1>{static_cast<Real>(i)} ;
    }
    return toreturn;
}
template<> VectorPoint<2> buildTestUmklappVects<2>(const int nNearestBrillouinZones, const Material<2>& material){
    VectorPoint<2> toreturn;
    const int numberTestUmklappVects = Utilities::intPower<2>(2*nNearestBrillouinZones+1);
    toreturn.resize(Eigen::NoChange, numberTestUmklappVects);
    int counter = 0;
    for (int i= - nNearestBrillouinZones; i <= nNearestBrillouinZones; ++i){
        for (int j= - nNearestBrillouinZones; j <= nNearestBrillouinZones; ++j){
            toreturn.col(counter++) = Point<2>{static_cast<Real>(i), static_cast<Real>(j)};
        }
    }
    return toreturn;
}
template<> VectorPoint<3> buildTestUmklappVects<3>(const int nNearestBrillouinZones, const Material<3>& material){
    VectorPoint<3> toreturn;
    const int numberTestUmklappVects = Utilities::intPower<3>(2*nNearestBrillouinZones+1);
    toreturn.resize(Eigen::NoChange, numberTestUmklappVects);
    int counter = 0;
    for (int i= - nNearestBrillouinZones; i <= nNearestBrillouinZones; ++i){
        for (int j= - nNearestBrillouinZones; j <= nNearestBrillouinZones; ++j){
            for (int k= - nNearestBrillouinZones; k <= nNearestBrillouinZones; ++k){
                toreturn.col(counter++) = Point<3> {static_cast<Real>(i), static_cast<Real>(j), static_cast<Real>(k)};
            }
        }
    }
    return toreturn;
}

template<int NDim, int Nlegs> bool ScatteringChannel<NDim,Nlegs>::findUmklappVectors(const int nNearestBrillouinZones){
    // Below we construct all the umklapp vectors we want to try by hand
    numberScatteringElements = 0;
    VectorPoint<NDim> testUmklappVects;
    testUmklappVects = buildTestUmklappVects<NDim>(nNearestBrillouinZones, *material);
    // Let us do an analysis of the scattering, simply counting the number of combinations allowed by energy and momentum conservation
    // We need to do it for each possible umklapp scattering
    VectorPoint<NDim> umklappVects;       // We will accumulate the umklapp vectors that produce scatterings here
    for (int i=0; i<testUmklappVects.cols(); ++i) {
                    
        // Calculate the number of allowed combinations for the test umklapp vector
//        int nScattElem = countScatteringElementIDsSingleUmklapp(testUmklappVects.col(i));
                    
//        if (nScattElem > 0 ) { // If there are allowed combinations, we add the vector
            umklappVects.conservativeResize(Eigen::NoChange, umklappVects.cols()+1);
            umklappVects.col(umklappVects.cols()-1) = testUmklappVects.col(i);
//            numberScatteringElements += nScattElem;
//        }
    }
    umklappUnitVectors = umklappVects;
    return (umklappUnitVectors.cols() != 0 );
}

//=======================================================
// Scattering Analysers
//===================
// Deterministic
template<int NDim, int NLegs> unsigned long ScatteringChannel<NDim,NLegs>::countScatteringElementIDsSingleUmklapp(const Point<NDim>& umklappUnitVect) const{
    int toreturn = 0;
    deterministicIteration(umklappUnitVect, [&toreturn](const std::vector<MeshSubsetElementIterator<NDim>>& elemID){++toreturn;});
    return toreturn;
}

template<int NDim, int NLegs> unsigned long ScatteringChannel<NDim,NLegs>::countScatteringElementIDs() {
    int toreturn = 0;
    int nUmklVecs = umklappUnitVectors.cols();
    for (int i=0; i<nUmklVecs; ++i){
      deterministicIteration(umklappUnitVectors.col(i), [&toreturn](const std::vector<MeshSubsetElementIterator<NDim>>& elemID){++toreturn;});
    }
    numberScatteringElements = toreturn;
    return toreturn;
}
template<int NDim, int NLegs> std::vector<std::array<int,NLegs>> ScatteringChannel<NDim,NLegs>::scatteringElementIDs() const {
    std::vector<std::array<int,NLegs>> toreturn;
    int nUmklVecs = umklappUnitVectors.cols();
    for (int i=0; i<nUmklVecs; ++i){
        deterministicIteration(umklappUnitVectors.col(i),
            [&toreturn](const std::vector<MeshSubsetElementIterator<NDim>>& elemID){
                std::array<int,NLegs> indexList;
                for (int i=0; i<NLegs; ++i){
                    indexList[i] = elemID[i].currentElementID;
                }
                toreturn.emplace_back(indexList);
            });
    }
    return toreturn;
}
template<int NDim, int NLegs> std::array<std::vector<int>,NLegs> ScatteringChannel<NDim,NLegs>::scatteringPerElement() const {
    std::array<std::vector<int>,NLegs> toreturn;    // Will contain the number of scatterings involving that element
    for (int i=0; i<NLegs; ++i) {
        toreturn[i].reserve(this->legBand[i]->mesh->numberElements);
        std::fill(toreturn[i].begin(), toreturn[i].end(), 0);    // Initialises to 0
    }
    int nUmklVecs = umklappUnitVectors.cols();
    for (int i=0; i<nUmklVecs; ++i){
        deterministicIteration(umklappUnitVectors.col(i),
            [&toreturn](const std::vector<MeshSubsetElementIterator<NDim>>& elemID){
            for (int i=0; i<NLegs; ++i) {
                ++(toreturn[i][elemID[i].currentElementID]);
            }
        });
    }
    return toreturn;
}
template<int NDim, int NLegs> Function<NDim> ScatteringChannel<NDim,NLegs>::scatteringPerElementPlot(const int which) const{
    auto elemIds = this->scatteringElementIDs();
    Function<NDim> nScattNumb (this->legBand[which]->mesh);
    for (int i = 0; i<elemIds.size(); ++i) { nScattNumb += FunctionElement<NDim>(this->legBand[which]->mesh,elemIds[i][which],1.0);}
    return nScattNumb;
}


//=======================================================
// Scattering Integrators setup
//===================

//=======================================================
// Methods for the construction of the scattering integrals
//===================

// Constructs the D matrix
template<int NDim, int Nlegs> auto constructD(const ScatteringChannel<NDim,Nlegs>& scattChan, const std::vector<MeshSubsetElementIterator<NDim>>& elemID){
    // Construction of D
    //         leg 0     leg 1     ...   leg NLEGS-1
    //      [         |          |      |           ]    signed momentum conservation dim 0
    //      [         |          |      |           ]        ...
    // D =  [         |          |      |           ]    signe momentum conservation dim NDIMS - 1
    //      [ --------+----------+------+---------- ]
    //      [         |          |      |           ]    signed energy conservation
    Eigen::Matrix<Real, NDim+1,Nlegs*NDim> D;
    for (int leg=0; leg<Nlegs; ++leg) {
        for (int dim=0; dim<NDim; ++dim){
            D.template block<1,NDim>(dim,leg*NDim) = scattChan.legDirection[leg] * scattChan.legBand[leg]->crystalMomentum[dim].elementVec(elemID[leg]).template block<1,NDim>(0,1);
        }
        D.template block<1,NDim>(NDim,leg*NDim) = scattChan.legDirection[leg] * scattChan.legBand[leg]->dispersion.elementVec(elemID[leg]).template block<1,NDim>(0,1);
    }
    return D;
}
template<int NDim, int Nlegs> auto constructD(const ScatteringChannel<NDim,Nlegs>& scattChan, const std::array<int,Nlegs>& elemID){
    // Construction of D
    //         leg 0     leg 1     ...   leg NLEGS-1
    //      [         |          |      |           ]    signed momentum conservation dim 0
    //      [         |          |      |           ]        ...
    // D =  [         |          |      |           ]    signe momentum conservation dim NDIMS - 1
    //      [ --------+----------+------+---------- ]
    //      [         |          |      |           ]    signed energy conservation
    Eigen::Matrix<Real, NDim+1,Nlegs*NDim> D;
    for (int leg=0; leg<Nlegs; ++leg) {
        for (int dim=0; dim<NDim; ++dim){
            D.template block<1,NDim>(dim,leg*NDim) = scattChan.legDirection[leg] * scattChan.legBand[leg]->crystalMomentum[dim].elementVec(elemID[leg]).template block<1,NDim>(0,1);
        }
        D.template block<1,NDim>(NDim,leg*NDim) = scattChan.legDirection[leg] * scattChan.legBand[leg]->dispersion.elementVec(elemID[leg]).template block<1,NDim>(0,1);
    }
    return D;
}

// Constructs the N vector
template<int NDim, int Nlegs> auto constructN(const ScatteringChannel<NDim,Nlegs>& scattChan, const Point<NDim>& umklappVect, const std::vector<MeshSubsetElementIterator<NDim>>& elemID){
    Eigen::Matrix<Real, NDim+1, 1> N( Eigen::Matrix<Real, NDim+1, 1>::Zero());
    for (int leg=0; leg<Nlegs; ++leg) {
        for (int dim=0; dim<NDim; ++dim){
            N(dim,0) -= scattChan.legDirection[leg] * scattChan.legBand[leg]->crystalMomentum[dim].elementVec(elemID[leg])(0);
        }
        N(NDim,0) -= scattChan.legDirection[leg] * scattChan.legBand[leg]->dispersion.elementVec(elemID[leg])(0);
    }
    
    for (int dim=0; dim<NDim; ++dim){
        for (int dimumkl=0; dimumkl<NDim; ++dimumkl){
            N(dim,0)-= umklappVect(dimumkl) * scattChan.material->region->gVec(dim,dimumkl);
        }
    }
    return N;
}
template<int NDim, int Nlegs> auto constructN(const ScatteringChannel<NDim,Nlegs>& scattChan, const Point<NDim>& umklappVect, const std::array<int,Nlegs>& elemID){
    Eigen::Matrix<Real, NDim+1, 1> N( Eigen::Matrix<Real, NDim+1, 1>::Zero());
    for (int leg=0; leg<Nlegs; ++leg) {
        for (int dim=0; dim<NDim; ++dim){
            N(dim,0) -= scattChan.legDirection[leg] * scattChan.legBand[leg]->crystalMomentum[dim].elementVec(elemID[leg])(0);
        }
        N(NDim,0) -= scattChan.legDirection[leg] * scattChan.legBand[leg]->dispersion.elementVec(elemID[leg])(0);
    }
    
    for (int dim=0; dim<NDim; ++dim){
        for (int dimumkl=0; dimumkl<NDim; ++dimumkl){
            N(dim,0)-= umklappVect(dimumkl) * scattChan.material->region->gVec(dim,dimumkl);
        }
    }
    return N;
}

// Constructs the outputformsLinForm
template<int NDim, int Nlegs> constexpr auto constructOutputformsLinForm(){
    // Construction of outputformsLinForm
    //    Example of the structure for the outputformsLinForm when NDIMS = 2, NLEGS = 4
    //    outputformsLinForm <<   1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
    //                            0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
    //                            0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0,
    //                            1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
    //                            1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0,
    //                            1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0,
    //                            1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
    //                            1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0,
    //                            1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0,
    //                            1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
    //                            1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0,
    //                            1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1;
    ArrayMultiLegLinearForm<NDim,Nlegs,Nlegs*(NDim+1)> outputformsLinForm(ArrayMultiLegLinearForm<NDim,Nlegs,Nlegs*(NDim+1)>::Zero());
    for (int leg0=0; leg0<Nlegs; ++leg0){
        for (int leg1=0; leg1<leg0; ++leg1){
            outputformsLinForm.template block<NDim+1,1>(leg0*(NDim+1),leg1*(NDim+1)) = Eigen::Matrix<Real, NDim+1, 1>::Ones();
        }
        outputformsLinForm.template block<NDim+1,NDim+1>(leg0*(NDim+1),leg0*(NDim+1)) = Eigen::Matrix<Real, NDim+1, NDim+1>::Identity();
        for (int leg1=leg0+1; leg1<Nlegs; ++leg1){
            outputformsLinForm.template block<NDim+1,1>(leg0*(NDim+1),leg1*(NDim+1)) = Eigen::Matrix<Real, NDim+1, 1>::Ones();
        }
    }
    return outputformsLinForm;
}
template<int NDim, int Nlegs> Real constructAmplitude(const ScatteringChannel<NDim,Nlegs>& scattChan, const std::vector<MeshSubsetElementIterator<NDim>>& elemID){
    // Constructs the amplitude by including two effects
    Real fullAmplitude = scattChan.amplitude;
    //  -> Ratio between the actual integration volume and the integration over the reference elements
    for (int leg=0; leg<Nlegs; ++leg) {
        fullAmplitude *= scattChan.legBand[leg]->mesh->elemVolume * static_cast<double>(Utilities::factorial(NDim)) ;
    }
    //  -> Evaluation of a scattering element amplitude at the baricentres of the elements
    ArrayPoint<NDim, Nlegs> baricentres;
    for (int i=0; i<Nlegs; ++i ){ baricentres.col(i) =  scattChan.legBand[i]->mesh->elemCentre(elemID[i]);}
//    std::array<Point<NDim>,Nlegs> baricentres;
//    for (int i=0; i<Nlegs; ++i ){ baricentres[i] =  scattChan.legBand[i]->mesh->elemCentre(elemID[i]);}
    fullAmplitude *= scattChan.functAmplitude(baricentres);
    
//    fullAmplitude *= scattChan.functAmplitude(elemID);
    return fullAmplitude;
}
template<int NDim, int Nlegs> Real constructAmplitude(const ScatteringChannel<NDim,Nlegs>& scattChan, const std::array<int,Nlegs>& elemID){
    // Constructs the amplitude by including two effects
    Real fullAmplitude = scattChan.amplitude;
    //  -> Ratio between the actual integration volume and the integration over the reference elements
    for (int leg=0; leg<Nlegs; ++leg) {
        fullAmplitude *= scattChan.legBand[leg]->mesh->elemVolume * static_cast<double>(Utilities::factorial(NDim)) ;
    }
    //  -> Evaluation of a scattering element amplitude at the baricentres of the elements
    ArrayPoint<NDim, Nlegs> baricentres;
    for (int i=0; i<Nlegs; ++i ){ baricentres.col(i) =  scattChan.legBand[i]->mesh->elemCentre(elemID[i]);}
//    std::array<Point<NDim>,Nlegs> baricentres;
//    for (int i=0; i<Nlegs; ++i ){ baricentres[i] =  scattChan.legBand[i]->mesh->elemCentre(elemID[i]);}
    fullAmplitude *= scattChan.functAmplitude(baricentres);
    
//    fullAmplitude *= scattChan.functAmplitude(elemID);
    return fullAmplitude;
}

//=======================================================
// Scattering Precalculator
//===================

template<int NDim, int Nlegs> void ScatteringChannel<NDim,Nlegs>::precalculateIndices(){
    
    std::ifstream  scatteringMatrixIndicesFileIn(scatteringName()+".scattIndices", std::ios::in | std::ios::binary);
    
    if (scatteringMatrixIndicesFileIn.good()){
        bool same = true;
        for (int i=0; i<Nlegs; ++i){
            for (int d=0; d<NDim; ++d){
                int  numSecSplit;
                scatteringMatrixIndicesFileIn.read((char*) (&numSecSplit), sizeof(int));
                same = same && (numSecSplit == legBand[i]->mesh->nSecSplits[d]);
            }
        }
        if (same){
            unsigned long totalSize;
            scatteringMatrixIndicesFileIn.read((char*) (&totalSize), sizeof(unsigned long) );
            
            scatteringMatrixIndices.resize(totalSize);
            scatteringMatrixIndicesCorrespondingUmklapp.resize(totalSize);
            
            for(unsigned long s = 0; s < totalSize; ++s){
                for (int i=0; i<Nlegs; ++i){
                    int index ;
                    scatteringMatrixIndicesFileIn.read((char*) (&index), sizeof(int) );
                    scatteringMatrixIndices[s][i] = index;
                }
                int umklappid;
                scatteringMatrixIndicesFileIn.read((char*) (&umklappid), sizeof(int) );
                scatteringMatrixIndicesCorrespondingUmklapp[s] = umklappid;
            }
            scatteringMatrixIndicesFileIn.close();
            
            std::cout << "Reading " << scatteringName() <<" from file\n";
            
            return;
        }
    }
    
    scatteringMatrixIndicesFileIn.close();
    
    numberNonZeroScatteringElements = 0;
    numberScatteringElementsOverThreshold = 0;
            
    Eigen::Matrix<Real, 1, Nlegs*(NDim+1)> commonLinForm;
    for(int j=0; j<Nlegs; ++j){
        commonLinForm.template block<1,NDim+1>(0,(NDim+1)*j) = legBand[j]->dispersion.fromModalToLinFormTransf.row(modalOrdersCombinations[0][j]);
    }
    
    int nUmklVecs = umklappUnitVectors.cols();
    
    std::deque<std::array<int,Nlegs>>     dequeScatteringMatrixIndices;
    std::deque<int>                       dequeScatteringMatrixIndicesCorrespondingUmklapp;
    
    for (int i=0; i<nUmklVecs; ++i){
        deterministicIteration(umklappUnitVectors.col(i),
            [this, i, &commonLinForm, &dequeScatteringMatrixIndices, &dequeScatteringMatrixIndicesCorrespondingUmklapp]
                               (const std::vector<MeshSubsetElementIterator<NDim>>& elemID){
            
            Eigen::Matrix<Real, NDim+1,Nlegs*NDim> D (constructD<NDim, Nlegs>(*this, elemID));
            
            // If the determinant of the scattering is too small, then it does not calculate the scattering
            if (std::fabs(integrationInversionDeterminant<NDim,Nlegs>(D)) < minimumDeterminant) { return ; }
            
            // Construction of N
            const Eigen::Matrix<Real, NDim+1, 1> N (constructN<NDim, Nlegs>(*this, umklappUnitVectors.col(i), elemID));
           
            // Constructs the outputformsLinForm
            const ArrayMultiLegLinearForm<NDim,Nlegs,Nlegs*(NDim+1)> outputformsLinForm(constructOutputformsLinForm<NDim, Nlegs>());
            
            // Constructs the amplitude by including two effects
            const Real fullAmplitude (constructAmplitude<NDim, Nlegs>(*this, elemID));

            // Perform actual integration
            Eigen::Matrix<Real, Nlegs*(NDim+1), Eigen::Dynamic> result = scatteringIntegrationTypeC<NDim, Nlegs>(D, N, commonLinForm, outputformsLinForm, this->numMCPoints).eval();
            result *= fullAmplitude;
            
            if (std::fabs(result(0,0))>scattElemeThreshold){
                ++numberScatteringElementsOverThreshold;
                dequeScatteringMatrixIndices.emplace_back(std::array<int,Nlegs>());
                for (int i=0; i<Nlegs; ++i){ dequeScatteringMatrixIndices.back()[i] = elemID[i].currentElementID;}
                dequeScatteringMatrixIndicesCorrespondingUmklapp.emplace_back(i);
            }
            if (std::fabs(result(0,0))>0.) ++numberNonZeroScatteringElements;
        }
        );
    }
    scatteringMatrixIndices.resize(dequeScatteringMatrixIndices.size());
    scatteringMatrixIndicesCorrespondingUmklapp.resize(dequeScatteringMatrixIndices.size());
    for(unsigned long i = 0; i < dequeScatteringMatrixIndices.size(); ++i){
        scatteringMatrixIndices[i] = dequeScatteringMatrixIndices[i];
    }
    for(unsigned long i = 0; i < dequeScatteringMatrixIndices.size(); ++i){
        scatteringMatrixIndicesCorrespondingUmklapp[i] = dequeScatteringMatrixIndicesCorrespondingUmklapp[i];
    }
    
    // Save to file
    std::ofstream scatteringMatrixIndicesFileOut(scatteringName()+".scattIndices", std::ios::out | std::ios::binary | std::ios::trunc);
    for (int i=0; i<Nlegs; ++i){
        for (int d=0; d<NDim; ++d){
            int  numSecSplit = legBand[i]->mesh->nSecSplits[d];
            scatteringMatrixIndicesFileOut.write((char*) (&numSecSplit), sizeof(int) );
        }
    }
    unsigned long totalSize = scatteringMatrixIndices.size();
    scatteringMatrixIndicesFileOut.write((char*) (&totalSize), sizeof(unsigned long) );
    for(unsigned long s = 0; s < scatteringMatrixIndices.size(); ++s){
        for (int i=0; i<Nlegs; ++i){
            int index = scatteringMatrixIndices[s][i];
            scatteringMatrixIndicesFileOut.write((char*) (&index), sizeof(int) );
        }
        int umklappid = scatteringMatrixIndicesCorrespondingUmklapp[s];
        scatteringMatrixIndicesFileOut.write((char*) (&umklappid), sizeof(int) );
    }
    scatteringMatrixIndicesFileOut.close();
    
    std::cout << "Writing " << scatteringName() <<" to file\n";
    
}

//=======================================================
// Scattering Rates
//===================

template<int NDim, int Nlegs> void ScatteringChannel<NDim,Nlegs>::singleScatteringRate(const Point<NDim>& umklappVect, const MaterialStatus<NDim>& status, MaterialStatus<NDim>& rates, const int which, const std::array<int,Nlegs>& elemID) const {
    
    Eigen::Matrix<Real, NDim+1,Nlegs*NDim> D (constructD<NDim, Nlegs>(*this, elemID));
    
    // If the determinant of the scattering is too small, then it does not calculate the scattering
    if (std::fabs(integrationInversionDeterminant<NDim,Nlegs>(D)) < minimumDeterminant) { return ; }
    
    // Construction of N
    const Eigen::Matrix<Real, NDim+1, 1> N (constructN<NDim, Nlegs>(*this, umklappVect, elemID));
   
    // Constructs the outputformsLinForm
    const ArrayMultiLegLinearForm<NDim,Nlegs,Nlegs*(NDim+1)> outputformsLinForm(constructOutputformsLinForm<NDim, Nlegs>());
    
    // Construction of commonLinForm0 and commonLinForm1
    ArrayMultiLegLinearForm<NDim,Nlegs,1> commonLinForm0, commonLinForm1;
    LinearForm<NDim> unit(LinearForm<NDim>::Zero());
    unit(0) = 1.0;
    for (int leg=0; leg<Nlegs; ++leg) {
        if (leg == which) {
            commonLinForm0.template block<1,NDim+1>(0,leg*(NDim+1)) = (leg==0?-1.0:1.0) * ( (1. - legDirection[leg]) * legBand[leg]->statistics + (1. + legDirection[leg]) ) * unit; // in
            commonLinForm1.template block<1,NDim+1>(0,leg*(NDim+1)) = ( (1. + legDirection[leg]) * legBand[leg]->statistics + (1. - legDirection[leg]) ) * unit; // out
        } else {
            commonLinForm0.template block<1,NDim+1>(0,leg*(NDim+1)) = (leg==0?-1.0:1.0) * 0.5 * ((1. - legDirection[leg]) * unit
                + ( (1. - legDirection[leg]) * legBand[leg]->statistics + (1. + legDirection[leg]) ) * status[legBandNumber[leg]].elementVec(elemID[leg]));    // Direct transition
            commonLinForm1.template block<1,NDim+1>(0,leg*(NDim+1)) = 0.5 * ( (1. + legDirection[leg]) * unit
                + ( (1. + legDirection[leg]) * legBand[leg]->statistics + (1. - legDirection[leg]) ) * status[legBandNumber[leg]].elementVec(elemID[leg]));    // Time reversed transition
        }
    }

    // Constructs the amplitude by including two effects
    const Real fullAmplitude (constructAmplitude<NDim, Nlegs>(*this, elemID));
    
    // Perform actual integration
    Eigen::Matrix<Real, Nlegs*(NDim+1), 1> result = scatteringIntegrationTypeB<NDim, Nlegs>(D, N, commonLinForm0, commonLinForm1, outputformsLinForm, mCPoints).eval();
    result *= fullAmplitude;
    // Returns only the rate for the selected particle
    rates[legBandNumber[which]].elementVec(elemID[which]) -= legDirection[which] * result.template block<NDim+1,1>(which*(NDim+1),0) ;
};
template<int NDim, int Nlegs> Function<NDim> ScatteringChannel<NDim,Nlegs>::scatteringRates(const MaterialStatus<NDim>& status, const int which) const {
    assert( (which>=0 ) && (which<Nlegs) );
    
    auto     rangetbb = tbb::blocked_range<unsigned long>(0, scatteringMatrixIndices.size());
    MaterialStatus<NDim> identity_value(status.material);
    auto const    toExecute = [&](tbb::blocked_range<unsigned long> r, MaterialStatus<NDim> running_rates) {
        for (unsigned long i=r.begin(); i<r.end(); ++i) {
            singleScatteringRate(umklappUnitVectors.col(scatteringMatrixIndicesCorrespondingUmklapp[i]), status, running_rates, which, scatteringMatrixIndices[i]);
        }
        return running_rates;
    };
    auto const     reduction = std::plus<MaterialStatus<NDim>>();
    auto rates = tbb::parallel_reduce(rangetbb, identity_value, toExecute, reduction);
    
//    MaterialStatus<NDim> rates(status.material);
//
//    for(unsigned long i = 0; i < scatteringMatrixIndices.size(); ++i){
//        singleScatteringRate(umklappUnitVectors.col(scatteringMatrixIndicesCorrespondingUmklapp[i]), status, rates, which, scatteringMatrixIndices[i]);
//    }
    
    return rates[legBandNumber[which]].applyInverseMass();
}
template<int NDim, int Nlegs> Function<NDim> ScatteringChannel<NDim,Nlegs>::scatteringRates(const MaterialStatus<NDim>& status, const std::string &bandName) const{
    Function<NDim> toreturn(status[bandName].mesh);
    for (int i=0; i<Nlegs; ++i){
        if (legBandNumber[i] == status.id(bandName)) {
            toreturn += scatteringRates(status, i);
        }
    }
    return toreturn;
}

//=======================================================
// Scattering Propagators
//===================

template<int NDim, int Nlegs> void ScatteringChannel<NDim,Nlegs>::singleScatteringPropagator(const Point<NDim>& umklappVect, const MaterialStatus<NDim>& status, MaterialStatus<NDim>& propagation, const std::array<int,Nlegs>& elemID) const {
    
    Eigen::Matrix<Real, NDim+1,Nlegs*NDim> D (constructD<NDim, Nlegs>(*this, elemID));
    
    // If the determinant of the scattering is too small, then it does not calculate the scattering
    if (std::fabs(integrationInversionDeterminant<NDim,Nlegs>(D)) < minimumDeterminant) { return ; }
    
    // Construction of N
    const Eigen::Matrix<Real, NDim+1, 1> N (constructN<NDim, Nlegs>(*this, umklappVect, elemID));
   
    // Constructs the outputformsLinForm
    const ArrayMultiLegLinearForm<NDim,Nlegs,Nlegs*(NDim+1)> outputformsLinForm(constructOutputformsLinForm<NDim, Nlegs>());
    
    // Construction of commonLinForm0 and commonLinForm1
    ArrayMultiLegLinearForm<NDim,Nlegs,1> commonLinForm0, commonLinForm1;
    LinearForm<NDim> unit(LinearForm<NDim>::Zero());
    unit(0) = 1.0;
    for (int leg=0; leg<Nlegs; ++leg) {
        commonLinForm0.template block<1,NDim+1>(0,leg*(NDim+1)) = (leg==0?-1.0:1.0) * 0.5 * ( (1. - legDirection[leg]) * unit
            + ( (1. - legDirection[leg]) * legBand[leg]->statistics + (1. + legDirection[leg]) ) * status[legBandNumber[leg]].elementVec(elemID[leg]));    // Direct transition
        commonLinForm1.template block<1,NDim+1>(0,leg*(NDim+1)) = 0.5 * ( (1. + legDirection[leg]) * unit
            + ( (1. + legDirection[leg]) * legBand[leg]->statistics + (1. - legDirection[leg]) ) * status[legBandNumber[leg]].elementVec(elemID[leg]));    // Time reversed transition
    }

    // Constructs the amplitude by including two effects
    const Real fullAmplitude (constructAmplitude<NDim, Nlegs>(*this, elemID));
    
    // Perform actual integration
    Eigen::Matrix<Real, Nlegs*(NDim+1), 1> result = scatteringIntegrationTypeB<NDim, Nlegs>(D, N, commonLinForm0, commonLinForm1, outputformsLinForm, mCPoints).eval();
    result *= fullAmplitude;
    
    for (int i=0; i<Nlegs; ++i ){ propagation[legBandNumber[i]].elementVec(elemID[i]) += legDirection[i] * result.template block<NDim+1,1>(i*(NDim+1),0);}
};
template<int NDim, int Nlegs> MaterialStatus<NDim> ScatteringChannel<NDim,Nlegs>::propagate(const MaterialStatus<NDim>& status) const {

    auto                    rangetbb = tbb::blocked_range<unsigned long>(0, scatteringMatrixIndices.size());
    MaterialStatus<NDim>    identity_value(status.material);
    auto const              toExecute = [&](tbb::blocked_range<unsigned long> r, MaterialStatus<NDim> running_propagation) {
        for (unsigned long i=r.begin(); i<r.end(); ++i) {
            singleScatteringPropagator(umklappUnitVectors.col(scatteringMatrixIndicesCorrespondingUmklapp[i]), status, running_propagation, scatteringMatrixIndices[i]);
        }
        return running_propagation;
    };
    auto const              reduction = std::plus<MaterialStatus<NDim>>();
    auto                    propagation = tbb::parallel_reduce(rangetbb, identity_value, toExecute, reduction);
 
//    MaterialStatus<NDim> propagation(status.material);
//    for(unsigned long i = 0; i < scatteringMatrixIndices.size(); ++i){
//        singleScatteringPropagator(umklappUnitVectors.col(scatteringMatrixIndicesCorrespondingUmklapp[i]), status, propagation, scatteringMatrixIndices[i]);
//    }
    
    return propagation.applyInverseMass();
}

    
//=======================================================
// Monte Carlo Points constructor
//===================

template<int NDim, int Nlegs> void ScatteringChannel<NDim,Nlegs>::constructMCPoints(){
    mCPoints.resize(Eigen::NoChange, numMCPoints);
    for (int i = 0; i < Nlegs - 2; ++i){
        mCPoints.template middleRows<NDim>(i*NDim) = randomPointsReference<NDim>(numMCPoints);
    }
    // Notice that there is one row in mCPoints that contains uninitialised data
    // That was done as a horrible workaround since mCPoints would have 0 rows for NDim = 1 and NLegs = 2
    // Eigen does not accept empty matrices
    if constexpr(NDim>1){
        mCPoints.template bottomRows<NDim-1>() = randomPointsReference<NDim-1>(numMCPoints);
    }
}


//=======================================================
// Technicalities
//===================
template<int NDim, int NLegs>  std::string ScatteringChannel<NDim,NLegs>::scatteringName() const {
    std::string name("[");
    for (int i=0; i<NLegs; ++i){ if (legDirection[i]>0) { name += "("+material->id(legBandNumber[i])+")" ;} }
    name += "]-[";
    for (int i=0; i<NLegs; ++i){ if (legDirection[i]<0) { name += "("+material->id(legBandNumber[i])+")" ;} }
    name += "]";
    return name;
};

template<int NDim, int NLegs>  void ScatteringChannel<NDim,NLegs>::printmodalOrdersCombinations(){
    for (int i=0; i<modalOrdersCombinations.size(); ++i){
        std::cout << " [";
        for (int j=0; j<NLegs-1; ++j){ std::cout << modalOrdersCombinations[i][j]<<",";}
        std::cout << modalOrdersCombinations[i][NLegs-1]<<"]   ";
    }
    std::cout << "\n";
}


// Explicit template instantiation (For more info see: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl)
template class  ScatteringChannel<1,2>;
template class  ScatteringChannel<1,3>;
template class  ScatteringChannel<1,4>;
template class  ScatteringChannel<2,2>;
template class  ScatteringChannel<2,3>;
template class  ScatteringChannel<2,4>;
template class  ScatteringChannel<3,2>;
template class  ScatteringChannel<3,3>;
template class  ScatteringChannel<3,4>;

} // namespace Tortoise
