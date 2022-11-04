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
//  Material.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//

#include <Physics/Material/Material.hpp>
#include <Generics/Utilities/StringUtilities.hpp>

#include <cassert>
#include <algorithm>

namespace Tortoise {

//=======================================================
// Constructors
//===================
template<int NDim> Material<NDim>::Material(const std::string &t_id, const Region<NDim>& t_region): matID(t_id), region(&t_region){};

//=======================================================
// Material construction: Bands
//===================
template<int NDim> void Material<NDim>::addBand(const std::string &t_id, const Real t_charge, const Real t_statistics, const Function<NDim>& t_dispersion){
    assert(t_dispersion.mesh->region == region);    // Makes sure the band is defined onto the same region
    this->emplace_back(t_id, Band<NDim>(t_charge, t_statistics, t_dispersion));
}
template<int NDim> void Material<NDim>::addBand(const std::string &t_id, const Real t_charge, const Real t_statistics, const Mesh<NDim> &t_mesh, const std::function<Real(Point<NDim>)>& t_dispersion) {
    assert(t_mesh.region == region);    // Makes sure the band is defined onto the same region
    this->emplace_back(t_id, Band<NDim>(t_charge, t_statistics, t_mesh, t_dispersion));
}
template<int NDim> void Material<NDim>::addBand(const std::string &t_id, const Real t_charge, const Real t_statistics, const Mesh<NDim> &t_mesh, const Real t_dispersion) {
    assert(t_mesh.region == region);    // Makes sure the band is defined onto the same region
    this->emplace_back(t_id, Band<NDim>(t_charge, t_statistics, t_mesh, t_dispersion));
}

//=======================================================
// Material construction: Scattering Channels
//===================

template<int NDim> template<int NLegs> std::string Material<NDim>::scatteringName(const std::array<int, NLegs>& bandsNumbers, const std::array<scattLegDirection, NLegs>& t_legDirection) const {
    std::string name("[");
    for (int i=0; i<NLegs; ++i){ if (t_legDirection[i]>0) { name += "("+this->id(bandsNumbers[i])+")" ;} }
    name += "]-[";
    for (int i=0; i<NLegs; ++i){ if (t_legDirection[i]<0) { name += "("+this->id(bandsNumbers[i])+")" ;} }
    name += "]";
    return name;
}
template<int NDim> template<int NLegs> std::string Material<NDim>::scatteringName(const std::array<std::string, NLegs>& bandsNames, const std::array<scattLegDirection, NLegs>& t_legDirection) const {
    std::string name("[");
    for (int i=0; i<NLegs; ++i){ if (t_legDirection[i]>0) { name += "("+bandsNames[i]+")" ;} }
    name += "]-[";
    for (int i=0; i<NLegs; ++i){ if (t_legDirection[i]<0) { name += "("+bandsNames[i]+")" ;} }
    name += "]";
    return name;
}

#define CONSTRUCT_SCATTERINGCHANNELINT(_NDIM_,_NLEGS_)                                                          \
template<> template<> ScatteringChannel<_NDIM_,_NLEGS_>* Material<_NDIM_>::scatteringChannel<_NLEGS_>           \
 (const std::array<int, _NLEGS_>& bandsNumbers, const std::array<scattLegDirection, _NLEGS_>& t_legDirection) { \
    return &(scatteringChannelsList##_NLEGS_##Legs[scatteringName<_NLEGS_>(bandsNumbers,t_legDirection)]);      \
}
CONSTRUCT_SCATTERINGCHANNELINT(1,2)
CONSTRUCT_SCATTERINGCHANNELINT(1,3)
CONSTRUCT_SCATTERINGCHANNELINT(1,4)
CONSTRUCT_SCATTERINGCHANNELINT(2,2)
CONSTRUCT_SCATTERINGCHANNELINT(2,3)
CONSTRUCT_SCATTERINGCHANNELINT(2,4)
CONSTRUCT_SCATTERINGCHANNELINT(3,2)
CONSTRUCT_SCATTERINGCHANNELINT(3,3)
CONSTRUCT_SCATTERINGCHANNELINT(3,4)


#define CONSTRUCT_SCATTERINGCHANNELSTRING(_NDIM_,_NLEGS_)                                                               \
template<> template<> ScatteringChannel<_NDIM_,_NLEGS_>* Material<_NDIM_>::scatteringChannel<_NLEGS_>                   \
(const std::array<std::string, _NLEGS_>& bandsNames, const std::array<scattLegDirection, _NLEGS_>& t_legDirection) {    \
    return &(scatteringChannelsList##_NLEGS_##Legs[scatteringName<_NLEGS_>(bandsNames,t_legDirection)]);                \
}
CONSTRUCT_SCATTERINGCHANNELSTRING(1,2)
CONSTRUCT_SCATTERINGCHANNELSTRING(1,3)
CONSTRUCT_SCATTERINGCHANNELSTRING(1,4)
CONSTRUCT_SCATTERINGCHANNELSTRING(2,2)
CONSTRUCT_SCATTERINGCHANNELSTRING(2,3)
CONSTRUCT_SCATTERINGCHANNELSTRING(2,4)
CONSTRUCT_SCATTERINGCHANNELSTRING(3,2)
CONSTRUCT_SCATTERINGCHANNELSTRING(3,3)
CONSTRUCT_SCATTERINGCHANNELSTRING(3,4)


#define CONSTRUCT_SCATTERINGCHANNELMATCHLIST(_NDIM_)                                                                \
template<> template<> std::vector<ScatteringChannel<_NDIM_,2>*> Material<_NDIM_>::scatteringChannelMatchList<2>     \
 (const std::array<std::string, 2>& bandsNamesWildCards, const std::array<scattLegDirection, 2>& t_legDirection){   \
    std::vector<ScatteringChannel<_NDIM_,2>*> toreturn;                                                             \
    for (auto band0 : matchListIndex(bandsNamesWildCards[0])){                                                      \
        for (auto band1 : matchListIndex(bandsNamesWildCards[1])){                                                  \
            auto i = scatteringChannelsList2Legs.id(scatteringName<2>({band0,band1},t_legDirection));               \
            if (i>=0) toreturn.emplace_back(scatteringChannel<2>({band0,band1},t_legDirection));                    \
        }                                                                                                           \
    }                                                                                                               \
    return toreturn;                                                                                                \
}                                                                                                                   \
template<> template<> std::vector<ScatteringChannel<_NDIM_,3>*> Material<_NDIM_>::scatteringChannelMatchList<3>     \
 (const std::array<std::string, 3>& bandsNamesWildCards, const std::array<scattLegDirection, 3>& t_legDirection){   \
    std::vector<ScatteringChannel<_NDIM_,3>*> toreturn;                                                             \
    for (auto band0 : matchListIndex(bandsNamesWildCards[0])){                                                      \
        for (auto band1 : matchListIndex(bandsNamesWildCards[1])){                                                  \
            for (auto band2 : matchListIndex(bandsNamesWildCards[2])){                                              \
                auto i = scatteringChannelsList3Legs.id(scatteringName<3>({band0,band1,band2},t_legDirection));     \
                if (i>=0) toreturn.emplace_back(scatteringChannel<3>({band0,band1,band2},t_legDirection));          \
            }                                                                                                       \
        }                                                                                                           \
    }                                                                                                               \
    return toreturn;                                                                                                \
}                                                                                                                   \
template<> template<> std::vector<ScatteringChannel<_NDIM_,4>*> Material<_NDIM_>::scatteringChannelMatchList<4>     \
 (const std::array<std::string, 4>& bandsNamesWildCards, const std::array<scattLegDirection, 4>& t_legDirection){   \
    std::vector<ScatteringChannel<_NDIM_,4>*> toreturn;                                                             \
    for (auto band0 : matchListIndex(bandsNamesWildCards[0])){                                                      \
        for (auto band1 : matchListIndex(bandsNamesWildCards[1])){                                                  \
            for (auto band2 : matchListIndex(bandsNamesWildCards[2])){                                              \
                for (auto band3 : matchListIndex(bandsNamesWildCards[3])){                                          \
                    auto i = scatteringChannelsList4Legs.id(scatteringName<4>({band0,band1,band2,band3},t_legDirection));   \
                    if (i>=0) toreturn.emplace_back(scatteringChannel<4>({band0,band1,band2,band3},t_legDirection));        \
                }                                                                                                   \
            }                                                                                                       \
        }                                                                                                           \
    }                                                                                                               \
    return toreturn;                                                                                                \
}
CONSTRUCT_SCATTERINGCHANNELMATCHLIST(1)
CONSTRUCT_SCATTERINGCHANNELMATCHLIST(2)
CONSTRUCT_SCATTERINGCHANNELMATCHLIST(3)

#define CONSTRUCT_ADDSCATTERINGCHANNEL(_NDIM_,_NLEGS_)                                                                                                  \
template<> template<> ScatteringChannel<_NDIM_,_NLEGS_>* Material<_NDIM_>::addScatteringChannel<_NLEGS_>                                                \
    (const std::array<int, _NLEGS_>& bandsNumbers, const std::array<scattLegDirection, _NLEGS_>& t_legDirection){                                       \
    scatteringChannelsList##_NLEGS_##Legs.emplace_back(scatteringName<_NLEGS_>(bandsNumbers, t_legDirection),                                           \
                                             ScatteringChannel<_NDIM_,_NLEGS_>(*this, bandsNumbers, t_legDirection));                                   \
    scatteringChannelsList##_NLEGS_##Legs[scatteringChannelsList##_NLEGS_##Legs.size()-1].numMCPoints = defaultNMCPointsScatterings;                    \
    scatteringChannelsList##_NLEGS_##Legs[scatteringChannelsList##_NLEGS_##Legs.size()-1].minimumDeterminant = defaultMinimumDeterminantScatterings;    \
    scatteringChannelsList##_NLEGS_##Legs[scatteringChannelsList##_NLEGS_##Legs.size()-1].scattElemeThreshold = defaultScattElemeThresholdScatterings;  \
    return &scatteringChannelsList##_NLEGS_##Legs[scatteringChannelsList##_NLEGS_##Legs.size()-1];                                                      \
}
CONSTRUCT_ADDSCATTERINGCHANNEL(1,2)
CONSTRUCT_ADDSCATTERINGCHANNEL(1,3)
CONSTRUCT_ADDSCATTERINGCHANNEL(1,4)
CONSTRUCT_ADDSCATTERINGCHANNEL(2,2)
CONSTRUCT_ADDSCATTERINGCHANNEL(2,3)
CONSTRUCT_ADDSCATTERINGCHANNEL(2,4)
CONSTRUCT_ADDSCATTERINGCHANNEL(3,2)
CONSTRUCT_ADDSCATTERINGCHANNEL(3,3)
CONSTRUCT_ADDSCATTERINGCHANNEL(3,4)


#define CONSTRUCT_ADDSCATTERINGCHANNELSTRING(_NDIM_,_NLEGS_)                                                                                            \
template<> template<> ScatteringChannel<_NDIM_,_NLEGS_>* Material<_NDIM_>::addScatteringChannel<_NLEGS_>                                                \
 (const std::array<std::string, _NLEGS_>& bandsNames, const std::array<scattLegDirection, _NLEGS_>& t_legDirection){                                    \
    scatteringChannelsList##_NLEGS_##Legs.emplace_back(scatteringName<_NLEGS_>(bandsNames, t_legDirection),                                             \
                                             ScatteringChannel<_NDIM_,_NLEGS_>(*this, bandsNames, t_legDirection));                                     \
    scatteringChannelsList##_NLEGS_##Legs[scatteringChannelsList##_NLEGS_##Legs.size()-1].numMCPoints = defaultNMCPointsScatterings;                    \
    scatteringChannelsList##_NLEGS_##Legs[scatteringChannelsList##_NLEGS_##Legs.size()-1].minimumDeterminant = defaultMinimumDeterminantScatterings;    \
    scatteringChannelsList##_NLEGS_##Legs[scatteringChannelsList##_NLEGS_##Legs.size()-1].scattElemeThreshold = defaultScattElemeThresholdScatterings;  \
    return &scatteringChannelsList##_NLEGS_##Legs[scatteringChannelsList##_NLEGS_##Legs.size()-1];                                                      \
}
CONSTRUCT_ADDSCATTERINGCHANNELSTRING(1,2)
CONSTRUCT_ADDSCATTERINGCHANNELSTRING(1,3)
CONSTRUCT_ADDSCATTERINGCHANNELSTRING(1,4)
CONSTRUCT_ADDSCATTERINGCHANNELSTRING(2,2)
CONSTRUCT_ADDSCATTERINGCHANNELSTRING(2,3)
CONSTRUCT_ADDSCATTERINGCHANNELSTRING(2,4)
CONSTRUCT_ADDSCATTERINGCHANNELSTRING(3,2)
CONSTRUCT_ADDSCATTERINGCHANNELSTRING(3,3)
CONSTRUCT_ADDSCATTERINGCHANNELSTRING(3,4)


#define CONSTRUCT_ADDSCATTERINGCHANNELSWILDCARDS2LEGS(_NDIM_)                                                                                                                   \
template<> template<> std::vector<ScatteringChannel<_NDIM_,2>*> Material<_NDIM_>::addScatteringChannelsWildcards<2>                                                             \
 (const std::array<std::string, 2>& bandsNames, const std::array<scattLegDirection, 2>& t_legDirection, int numUmklappSteps){                                                   \
    for (auto i : matchListId(bandsNames[0])) {                                                                                                                                 \
        for (auto j : matchListId(bandsNames[1])) {                                                                                                                             \
        /* Check for exchange and time reversal symmetry */                                                                                                                     \
            bool foundSame = false;                                                                                                                                             \
            for (int l = 0; l < scatteringChannelsList2Legs.size(); ++l){                                                                                                       \
                std::array<std::string, 2> bandNAmesToAdd = {i,j};                                                                                                              \
                std::array<std::string, 2> bandNamesInScatt = {id(scatteringChannelsList2Legs[l].legBandNumber[0]),id(scatteringChannelsList2Legs[l].legBandNumber[1])};        \
                bool bandsPermuted = std::is_permutation(bandNamesInScatt.begin(), bandNamesInScatt.end(), bandNAmesToAdd.begin(), bandNAmesToAdd.end());                       \
                bool directionsPerm = (t_legDirection[0]*t_legDirection[1]) == (scatteringChannelsList2Legs[l].legDirection[0]*scatteringChannelsList2Legs[l].legDirection[1]); \
                if (bandsPermuted and directionsPerm) { foundSame = true; break;}                                                                                               \
            }                                                                                                                                                                   \
            if (foundSame) continue;                                                                                                                                            \
            scatteringChannelsList2Legs.emplace_back(scatteringName<2>({i,j}, t_legDirection), ScatteringChannel<_NDIM_,2>(*this, {i,j}, t_legDirection));                      \
            scatteringChannelsList2Legs[scatteringChannelsList2Legs.size()-1].numMCPoints = defaultNMCPointsScatterings;                                                        \
            scatteringChannelsList2Legs[scatteringChannelsList2Legs.size()-1].minimumDeterminant = defaultMinimumDeterminantScatterings;                                        \
            scatteringChannelsList2Legs[scatteringChannelsList2Legs.size()-1].scattElemeThreshold = defaultScattElemeThresholdScatterings;                                      \
            scatteringChannelsList2Legs[scatteringChannelsList2Legs.size()-1].findUmklappVectors(numUmklappSteps);                                                              \
            int nScattElem = scatteringChannelsList2Legs[scatteringChannelsList2Legs.size()-1].countScatteringElementIDs();                                                     \
            if (nScattElem == 0 ) { scatteringChannelsList2Legs.pop_back();}  /* If there are no allowed scatterings in the channel, we remove it */                            \
        }                                                                                                                                                                       \
    }                                                                                                                                                                           \
    return scatteringChannelMatchList<2>(bandsNames, t_legDirection);                                                                                                           \
}
CONSTRUCT_ADDSCATTERINGCHANNELSWILDCARDS2LEGS(1)
CONSTRUCT_ADDSCATTERINGCHANNELSWILDCARDS2LEGS(2)
CONSTRUCT_ADDSCATTERINGCHANNELSWILDCARDS2LEGS(3)

#define CONSTRUCT_ADDSCATTERINGCHANNELSWILDCARDS3LEGS(_NDIM_)                                                                                                                   \
template<> template<> std::vector<ScatteringChannel<_NDIM_,3>*> Material<_NDIM_>::addScatteringChannelsWildcards<3>                                                             \
(const std::array<std::string, 3>& bandsNames, const std::array<scattLegDirection, 3>& t_legDirection, int numUmklappSteps){                                                    \
    for (auto i : matchListId(bandsNames[0])) {                                                                                                                                 \
        for (auto j : matchListId(bandsNames[1])) {                                                                                                                             \
            for (auto k : matchListId(bandsNames[2])) {                                                                                                                         \
                /* Check for exchange and time reversal symmetry */                                                                                                             \
                bool foundSame = false;                                                                                                                                         \
                std::vector<std::string> bandNamesToAddGroupIn, bandNamesToAddGroupOut;                                                                                         \
                if (t_legDirection[0] == in) bandNamesToAddGroupIn.emplace_back(i); else bandNamesToAddGroupOut.emplace_back(i);                                                \
                if (t_legDirection[1] == in) bandNamesToAddGroupIn.emplace_back(j); else bandNamesToAddGroupOut.emplace_back(j);                                                \
                if (t_legDirection[2] == in) bandNamesToAddGroupIn.emplace_back(k); else bandNamesToAddGroupOut.emplace_back(k);                                                \
                for (int l = 0; l < scatteringChannelsList3Legs.size(); ++l){                                                                                                   \
                    std::vector<std::string> bandNamesGroupIn, bandNamesGroupOut;                                                                                               \
                    if (scatteringChannelsList3Legs[l].legDirection[0] == in)                                                                                                   \
                        bandNamesGroupIn.emplace_back(id(scatteringChannelsList3Legs[l].legBandNumber[0]));                                                                     \
                    else                                                                                                                                                        \
                        bandNamesGroupOut.emplace_back(id(scatteringChannelsList3Legs[l].legBandNumber[0]));                                                                    \
                    if (scatteringChannelsList3Legs[l].legDirection[1] == in)                                                                                                   \
                        bandNamesGroupIn.emplace_back(id(scatteringChannelsList3Legs[l].legBandNumber[1]));                                                                     \
                    else                                                                                                                                                        \
                        bandNamesGroupOut.emplace_back(id(scatteringChannelsList3Legs[l].legBandNumber[1]));                                                                    \
                    if (scatteringChannelsList3Legs[l].legDirection[2] == in)                                                                                                   \
                        bandNamesGroupIn.emplace_back(id(scatteringChannelsList3Legs[l].legBandNumber[2]));                                                                     \
                    else                                                                                                                                                        \
                        bandNamesGroupOut.emplace_back(id(scatteringChannelsList3Legs[l].legBandNumber[2]));                                                                    \
                    if ( bandNamesToAddGroupIn.size() == bandNamesGroupIn.size() &&                                                                                             \
                        std::is_permutation(bandNamesToAddGroupIn.begin(), bandNamesToAddGroupIn.end(), bandNamesGroupIn.begin(), bandNamesGroupIn.end()) &&                    \
                        std::is_permutation(bandNamesToAddGroupOut.begin(), bandNamesToAddGroupOut.end(), bandNamesGroupOut.begin(), bandNamesGroupOut.end())){                 \
                        foundSame = true; break;                                                                                                                                \
                    }                                                                                                                                                           \
                    if ( bandNamesToAddGroupIn.size() == bandNamesGroupOut.size() &&                                                                                            \
                        std::is_permutation(bandNamesToAddGroupIn.begin(), bandNamesToAddGroupIn.end(), bandNamesGroupOut.begin(), bandNamesGroupOut.end()) &&                  \
                        std::is_permutation(bandNamesToAddGroupOut.begin(), bandNamesToAddGroupOut.end(), bandNamesGroupIn.begin(), bandNamesGroupIn.end())){                   \
                        foundSame = true; break;                                                                                                                                \
                    }                                                                                                                                                           \
                }                                                                                                                                                               \
                if (foundSame) continue;                                                                                                                                        \
                scatteringChannelsList3Legs.emplace_back(scatteringName<3>({i,j,k}, t_legDirection), ScatteringChannel<_NDIM_,3>(*this, {i,j,k}, t_legDirection));              \
                scatteringChannelsList3Legs[scatteringChannelsList3Legs.size()-1].numMCPoints = defaultNMCPointsScatterings;                                                    \
                scatteringChannelsList3Legs[scatteringChannelsList3Legs.size()-1].minimumDeterminant = defaultMinimumDeterminantScatterings;                                    \
                scatteringChannelsList3Legs[scatteringChannelsList3Legs.size()-1].scattElemeThreshold = defaultScattElemeThresholdScatterings;                                  \
                scatteringChannelsList3Legs[scatteringChannelsList3Legs.size()-1].findUmklappVectors(numUmklappSteps);                                                          \
                int nScattElem = scatteringChannelsList3Legs[scatteringChannelsList3Legs.size()-1].countScatteringElementIDs();                                                 \
                if (nScattElem == 0 ) { scatteringChannelsList3Legs.pop_back();}  /* If there are no allowed scatterings in the channel, we remove it */                        \
            }                                                                                                                                                                   \
        }                                                                                                                                                                       \
    }                                                                                                                                                                           \
    return scatteringChannelMatchList<3>(bandsNames, t_legDirection);                                                                                                           \
}
CONSTRUCT_ADDSCATTERINGCHANNELSWILDCARDS3LEGS(1)
CONSTRUCT_ADDSCATTERINGCHANNELSWILDCARDS3LEGS(2)
CONSTRUCT_ADDSCATTERINGCHANNELSWILDCARDS3LEGS(3)

#define CONSTRUCT_ADDSCATTERINGCHANNELSWILDCARDS4LEGS(_NDIM_)                                                                                                                   \
template<> template<> std::vector<ScatteringChannel<_NDIM_,4>*> Material<_NDIM_>::addScatteringChannelsWildcards<4>                                                             \
(const std::array<std::string, 4>& bandsNames, const std::array<scattLegDirection, 4>& t_legDirection, int numUmklappSteps){                                                    \
    for (auto i : matchListId(bandsNames[0])) {                                                                                                                                 \
        for (auto j : matchListId(bandsNames[1])) {                                                                                                                             \
            for (auto k : matchListId(bandsNames[2])) {                                                                                                                         \
                for (auto l : matchListId(bandsNames[3])) {                                                                                                                     \
                    /* Check for exchange and time reversal symmetry */                                                                                                         \
                    bool foundSame = false;                                                                                                                                     \
                    std::vector<std::string> bandNamesToAddGroupIn, bandNamesToAddGroupOut;                                                                                     \
                    if (t_legDirection[0] == in) bandNamesToAddGroupIn.emplace_back(i); else bandNamesToAddGroupOut.emplace_back(i);                                            \
                    if (t_legDirection[1] == in) bandNamesToAddGroupIn.emplace_back(j); else bandNamesToAddGroupOut.emplace_back(j);                                            \
                    if (t_legDirection[2] == in) bandNamesToAddGroupIn.emplace_back(k); else bandNamesToAddGroupOut.emplace_back(k);                                            \
                    if (t_legDirection[3] == in) bandNamesToAddGroupIn.emplace_back(l); else bandNamesToAddGroupOut.emplace_back(l);                                            \
                    for (int m = 0; m < scatteringChannelsList4Legs.size(); ++m){                                                                                               \
                        std::vector<std::string> bandNamesGroupIn, bandNamesGroupOut;                                                                                           \
                        if (scatteringChannelsList4Legs[m].legDirection[0] == in)                                                                                               \
                            bandNamesGroupIn.emplace_back(id(scatteringChannelsList4Legs[m].legBandNumber[0]));                                                                 \
                        else                                                                                                                                                    \
                            bandNamesGroupOut.emplace_back(id(scatteringChannelsList4Legs[m].legBandNumber[0]));                                                                \
                        if (scatteringChannelsList4Legs[m].legDirection[1] == in)                                                                                               \
                            bandNamesGroupIn.emplace_back(id(scatteringChannelsList4Legs[m].legBandNumber[1]));                                                                 \
                        else                                                                                                                                                    \
                            bandNamesGroupOut.emplace_back(id(scatteringChannelsList4Legs[m].legBandNumber[1]));                                                                \
                        if (scatteringChannelsList4Legs[m].legDirection[2] == in)                                                                                               \
                            bandNamesGroupIn.emplace_back(id(scatteringChannelsList4Legs[m].legBandNumber[2]));                                                                 \
                        else                                                                                                                                                    \
                            bandNamesGroupOut.emplace_back(id(scatteringChannelsList4Legs[m].legBandNumber[2]));                                                                \
                        if (scatteringChannelsList4Legs[m].legDirection[3] == in)                                                                                               \
                            bandNamesGroupIn.emplace_back(id(scatteringChannelsList4Legs[m].legBandNumber[3]));                                                                 \
                        else                                                                                                                                                    \
                            bandNamesGroupOut.emplace_back(id(scatteringChannelsList4Legs[m].legBandNumber[3]));                                                                \
                        if ( bandNamesToAddGroupIn.size() == bandNamesGroupIn.size() &&                                                                                         \
                            std::is_permutation(bandNamesToAddGroupIn.begin(), bandNamesToAddGroupIn.end(), bandNamesGroupIn.begin(), bandNamesGroupIn.end()) &&                \
                            std::is_permutation(bandNamesToAddGroupOut.begin(), bandNamesToAddGroupOut.end(), bandNamesGroupOut.begin(), bandNamesGroupOut.end())){             \
                            foundSame = true; break;                                                                                                                            \
                        }                                                                                                                                                       \
                        if ( bandNamesToAddGroupIn.size() == bandNamesGroupOut.size() &&                                                                                        \
                            std::is_permutation(bandNamesToAddGroupIn.begin(), bandNamesToAddGroupIn.end(), bandNamesGroupOut.begin(), bandNamesGroupOut.end()) &&              \
                            std::is_permutation(bandNamesToAddGroupOut.begin(), bandNamesToAddGroupOut.end(), bandNamesGroupIn.begin(), bandNamesGroupIn.end())){               \
                            foundSame = true; break;                                                                                                                            \
                        }                                                                                                                                                       \
                    }                                                                                                                                                           \
                    if (foundSame) continue;                                                                                                                                    \
                    scatteringChannelsList4Legs.emplace_back(scatteringName<4>({i,j,k,l}, t_legDirection), ScatteringChannel<_NDIM_,4>(*this, {i,j,k,l}, t_legDirection));      \
                    scatteringChannelsList4Legs[scatteringChannelsList4Legs.size()-1].numMCPoints = defaultNMCPointsScatterings;                                                \
                    scatteringChannelsList4Legs[scatteringChannelsList4Legs.size()-1].minimumDeterminant = defaultMinimumDeterminantScatterings;                                \
                    scatteringChannelsList4Legs[scatteringChannelsList4Legs.size()-1].scattElemeThreshold = defaultScattElemeThresholdScatterings;                              \
                    scatteringChannelsList4Legs[scatteringChannelsList4Legs.size()-1].findUmklappVectors(numUmklappSteps);                                                      \
                    int nScattElem = scatteringChannelsList4Legs[scatteringChannelsList4Legs.size()-1].countScatteringElementIDs();                                             \
                    if (nScattElem == 0 ) { scatteringChannelsList4Legs.pop_back();}  /* If there are no allowed scatterings in the channel, we remove it */                    \
                }                                                                                                                                                               \
            }                                                                                                                                                                   \
        }                                                                                                                                                                       \
    }                                                                                                                                                                           \
    return scatteringChannelMatchList<4>(bandsNames, t_legDirection);                                                                                                           \
}
CONSTRUCT_ADDSCATTERINGCHANNELSWILDCARDS4LEGS(1)
CONSTRUCT_ADDSCATTERINGCHANNELSWILDCARDS4LEGS(2)
CONSTRUCT_ADDSCATTERINGCHANNELSWILDCARDS4LEGS(3)


template <int NDim> template<int NLegs> Function<NDim> Material<NDim>::scatteringRates(const std::array<std::string, NLegs>& bandsNamesWildCards, const std::array<scattLegDirection, NLegs>& t_legDirection, const MaterialStatus<NDim>& status, const std::string& bandName) {
    Function<NDim> toReturn((*this)[bandName].mesh);
    for (auto scatt : scatteringChannelMatchList<NLegs>(bandsNamesWildCards, t_legDirection)){
        toReturn += scatt->scatteringRatesDeterministic(status, bandName);
        }
    return toReturn;
};

//=======================================================
// I/O
//===================
Plotter2D& operator << (Plotter2D& plotter, const Material<1>& material) {
    // RED, BLUE, LIME, FUCHSIA, AQUA, YELLOW, MAROON, NAVY, GREEN
    static std::vector<std::string> colors  {"FF0000", "0000FF", "00FF00", "FF00FF", "00FFFF", "FFFF00", "800000", "000080", "008000"};
    plotter2d << *material.region  << LABEL("BZ") ; // << COLOR("black") ;
    for (int band = 0; band < material.size(); ++band ) {
        plotter2d << *material[band].mesh << NOLABEL << COLOR("#AA" + colors[band%colors.size()]) ;
        plotter2d << material[band].dispersion << LABEL(material.id(band)) << COLOR("#00" + colors[band%colors.size()]) ;
    }
    return plotter;
}
Plotter3D& operator << (Plotter3D& plotter, const Material<2>& material) {
    // RED, BLUE, LIME, FUCHSIA, AQUA, YELLOW, MAROON, NAVY, GREEN
    static std::vector<std::string> colors  {"FF0000", "0000FF", "00FF00", "FF00FF", "00FFFF", "FFFF00", "800000", "000080", "008000"};
    plotter3d << *material.region  << LABEL("BZ") << COLOR("black") ;
    for (int band = 0; band < material.size(); ++band ) {
        plotter3d << *material[band].mesh << NOLABEL << COLOR("#AA" + colors[band%colors.size()]) ;
        plotter3d << material[band].dispersion << LABEL(material.id(band)) << COLOR("#00" + colors[band%colors.size()]) ;
    }
    return plotter;
}

template<> void Material<1>::plot() const {plotter2d << *this << PLOT;};
template<> void Material<2>::plot() const {plotter3d << *this << PLOT;};
template<> void Material<1>::plot(const std::vector<std::string> &bandstoplot) const {
    // RED, BLUE, LIME, FUCHSIA, AQUA, YELLOW, MAROON, NAVY, GREEN
    static std::vector<std::string> colors  {"FF0000", "0000FF", "00FF00", "FF00FF", "00FFFF", "FFFF00", "800000", "000080", "008000"};
    plotter2d << *region  << LABEL("BZ") << COLOR("black") ;
    for (int band = 0; band < bandstoplot.size(); ++band ) {
        plotter2d << *((*this)[bandstoplot[band]].mesh) << NOLABEL << COLOR("#AA" + colors[band%colors.size()]) ;
        plotter2d << (*this)[bandstoplot[band]].dispersion << LABEL(bandstoplot[band]) << COLOR("#00" + colors[band%colors.size()]) ;
    }
    plotter2d << PLOT;
};
template<> void Material<2>::plot(const std::vector<std::string> &bandstoplot) const {
    // RED, BLUE, LIME, FUCHSIA, AQUA, YELLOW, MAROON, NAVY, GREEN
    static std::vector<std::string> colors  {"FF0000", "0000FF", "00FF00", "FF00FF", "00FFFF", "FFFF00", "800000", "000080", "008000"};
    plotter3d << *region  << LABEL("BZ") << COLOR("black") ;
    for (int band = 0; band < bandstoplot.size(); ++band ) {
        plotter3d << *((*this)[bandstoplot[band]].mesh) << NOLABEL << COLOR("#AA" + colors[band%colors.size()]) ;
        plotter3d << (*this)[bandstoplot[band]].dispersion << LABEL(bandstoplot[band]) << COLOR("#00" + colors[band%colors.size()]) ;
    }
    plotter3d << PLOT;
};

template<int NDim> std::ostream &operator<<(std::ostream &os, Material<NDim> const& material){
    os <<  "********************************************\n";
    os <<  " Material: " << material.matID << "\n";
    os <<  "********************************************\n";
    os <<  " Band (n elems) |  Total number of bands " << material.size() << "\n";
    for (auto i : material.matchListIndex("*")) {
        os << std::left << std::setw(13) <<  material.id(i) << "  (" << material[i].mesh->numberElements << ")\n";
    }
    os <<  "********************************************\n";
    os <<  " Scattering channel (n elems) \n";
    os <<  "********************************************\n";
    os <<  " 2 - Legs | Total number of scattering channels " << material.scatteringChannelsList2Legs.size() << "\n";
    for (int i = 0; i < material.scatteringChannelsList2Legs.size(); ++i){
        int nscatt = material.scatteringChannelsList2Legs[i].getNumberScatteringElements();
        os << material.scatteringChannelsList2Legs.id(i) <<  ")   umklapp vectors : " << material.scatteringChannelsList2Legs[i].umklappUnitVectors.cols() << "   scatt elements " << (nscatt>0?std::to_string(nscatt):" ? ") << "   scatt amplitude " << material.scatteringChannelsList2Legs[i].amplitude<<  "\n";
    }
    os <<  " 3 - Legs | Total number of scattering channels " << material.scatteringChannelsList3Legs.size() << "\n";
    for (int i = 0; i < material.scatteringChannelsList3Legs.size(); ++i){
        int nscatt = material.scatteringChannelsList3Legs[i].getNumberScatteringElements();
        os << material.scatteringChannelsList3Legs.id(i) <<  ")   umklapp vectors : " << material.scatteringChannelsList3Legs[i].umklappUnitVectors.cols() << "   scatt elements " << (nscatt>0?std::to_string(nscatt):" ? ") << "   scatt amplitude " << material.scatteringChannelsList3Legs[i].amplitude<<  "\n";
    }
    os <<  " 4 - Legs | Total number of scattering channels " << material.scatteringChannelsList4Legs.size() << "\n";
    for (int i = 0; i < material.scatteringChannelsList4Legs.size(); ++i){
        int nscatt = material.scatteringChannelsList4Legs[i].getNumberScatteringElements();
        os << material.scatteringChannelsList4Legs.id(i) <<  ")   umklapp vectors : " << material.scatteringChannelsList4Legs[i].umklappUnitVectors.cols() << "   scatt elements " << (nscatt>0?std::to_string(nscatt):" ? ") << "   scatt amplitude " << material.scatteringChannelsList4Legs[i].amplitude <<  "\n";
    }
    return os;
}


// Explicit template instantiation (For more info see: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl)

#define EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_TEMPLATES_AND_FUNCTIONS(_NDIM_)        \
template class Material<_NDIM_>;                                                        \
template  std::ostream &operator<<(std::ostream &os, const Material<_NDIM_> & material);
EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_TEMPLATES_AND_FUNCTIONS(1)
EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_TEMPLATES_AND_FUNCTIONS(2)
EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_TEMPLATES_AND_FUNCTIONS(3)


#define EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_METHODS_TEMPLATES(_NDIM_,_NLEGS_)                                                              \
template std::string Material<_NDIM_>::scatteringName<_NLEGS_>                                                                                  \
 (const std::array<int, _NLEGS_>& bandsNumbers, const std::array<scattLegDirection, _NLEGS_>& t_legDirection) const;                            \
template std::string Material<_NDIM_>::scatteringName<_NLEGS_>                                                                                  \
 (const std::array<std::string, _NLEGS_>& bandsNames, const std::array<scattLegDirection, _NLEGS_>& t_legDirection) const;                      \
template Function<_NDIM_> Material<_NDIM_>::scatteringRates<_NLEGS_>                                                                            \
 (const std::array<std::string, _NLEGS_>&, const std::array<scattLegDirection, _NLEGS_>&, const MaterialStatus<_NDIM_>&, const std::string&);
EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_METHODS_TEMPLATES(1,2)
EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_METHODS_TEMPLATES(1,3)
EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_METHODS_TEMPLATES(1,4)
EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_METHODS_TEMPLATES(2,2)
EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_METHODS_TEMPLATES(2,3)
EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_METHODS_TEMPLATES(2,4)
EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_METHODS_TEMPLATES(3,2)
EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_METHODS_TEMPLATES(3,3)
EXPLICIT_TEMPLATE_INSTANTIATION_MATERIAL_METHODS_TEMPLATES(3,4)



} // namespace Tortoise
