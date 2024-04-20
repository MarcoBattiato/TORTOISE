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
//  Material.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//
// All bands within a material must operate over the same region. This is to ensure that the meshes are aligned


#ifndef Material_hpp
#define Material_hpp

#include <Physics/Material/Band.hpp>
#include <Physics/Scattering/ScatteringChannel.hpp>
#include <Generics/Containers/StringMap.hpp>

#include <vector>
#include <iostream>

namespace Tortoise {

enum statisticType {fermion = -1, boson = 1};

// As ScatteringChannel uses data from the Material class, and the Material class stores and constructs ScatteringChannel,
// it is necessary to provide forward declaration of the other classes
template <int NDim, int NLegs> class ScatteringChannel;
template<int NDim> class MaterialStatus;

template<int NDim> class Material :                 // NDim = Number of Dimensions
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Inheritance
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
// Makes the Material a container of Bands indexed by strings (band names)
// Single bands can be accessed by materialname["bandname"]
// Bands can be iterated using wildcards
// for (auto i : materialname.matchListIndex("electr*")) { materialname[i];}
// or
// for (auto i : materialname.matchList("electr*")) { *i;}
// IMPORTANT: append bands only using the method provided
public Containers::StringMap<Band<NDim>>

{
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main properties
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    const std::string               matID;          // Stores the ID of the material
    const Region<NDim>*             region;         // Origin of the reciprocal lattice unit cell (which we will incorrectly call BZ)

    Containers::StringMap<ScatteringChannel<NDim,2>> scatteringChannelsList2Legs;
    Containers::StringMap<ScatteringChannel<NDim,3>> scatteringChannelsList3Legs;
    Containers::StringMap<ScatteringChannel<NDim,4>> scatteringChannelsList4Legs;
    
    // Default values used when creating a new scattering. The scattering will be constructed with these default values,
    // yet, the values can be modified within each scattering
    int                             defaultNMCPointsScatterings = 50;    // Default value of MC points for newly created scatterings
    Real                            defaultMinimumDeterminantScatterings = 1.e-9;
    Real                            defaultScattElemeThresholdScatterings = 1.e-9;

    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main methods
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
public:
    //=======================================================
    // Constructors
    //===================
    Material(const std::string &t_matid, const Region<NDim>& t_region);
    
    //=======================================================
    // Material construction: Bands
    //===================
    Band<NDim>& addBand(const std::string &t_id, const Real t_charge, const Real t_statistics, const Function<NDim>& t_dispersion);
    Band<NDim>& addBand(const std::string &t_id, const Real t_charge, const Real t_statistics, const Mesh<NDim> &t_mesh, const std::function<Real(Point<NDim>)>& t_dispersion);
    Band<NDim>& addBand(const std::string &t_id, const Real t_charge, const Real t_statistics, const Mesh<NDim> &t_mesh, const Real t_dispersion);
    
    //=======================================================
    // Material construction: Scattering Channels
    //===================
    template<int NLegs> ScatteringChannel<NDim,NLegs>* addScatteringChannel(const std::array<std::string, NLegs>& bandsNames, const std::array<scattLegDirection, NLegs>& t_legDirection);
    template<int NLegs> ScatteringChannel<NDim,NLegs>* addScatteringChannel(const std::array<int, NLegs>& bandsNumbers, const std::array<scattLegDirection, NLegs>& t_legDirection);
    template<int NLegs> std::vector<ScatteringChannel<NDim,NLegs>*> addScatteringChannelsWildcards(const std::array<std::string, NLegs>& bandsNames, const std::array<scattLegDirection, NLegs>& t_legDirection, int numUmklappSteps);
    
    template<int NLegs> std::string scatteringName(const std::array<int, NLegs>& bandsNumbers, const std::array<scattLegDirection, NLegs>& t_legDirection) const;
    template<int NLegs> std::string scatteringName(const std::array<std::string, NLegs>& bandsNames, const std::array<scattLegDirection, NLegs>& t_legDirection) const;
    
    template<int NLegs> ScatteringChannel<NDim,NLegs>* scatteringChannel(const std::array<int, NLegs>& bandsNumbers, const std::array<scattLegDirection, NLegs>& t_legDirection) ;
    template<int NLegs> ScatteringChannel<NDim,NLegs>* scatteringChannel(const std::array<std::string, NLegs>& bandsNames, const std::array<scattLegDirection, NLegs>& t_legDirection) ;
    
    template<int NLegs> std::vector<ScatteringChannel<NDim,NLegs>*> scatteringChannelMatchList(const std::array<std::string, NLegs>& bandsNamesWildCards, const std::array<scattLegDirection, NLegs>& t_legDirection);
    
    //=======================================================
    // Scattering Initialisation
    //===================
    void constructMCPoints();
    
    //=======================================================
    // Material Analysis
    //===================
    template<int NLegs> Function<NDim> scatteringRates(const std::array<std::string, NLegs>& bandsNamesWildCards, const std::array<scattLegDirection, NLegs>& t_legDirection, const MaterialStatus<NDim>& status, const std::string& bandName);
    template<int NLegs> Function<NDim> selectivePropagation(const std::array<std::string, NLegs>& bandsNamesWildCards, const std::array<scattLegDirection, NLegs>& t_legDirection, const MaterialStatus<NDim>& status, const std::string& bandName);
    
    //=======================================================
    // Utilities
    //===================
    // In case of typical applications in Materials and for light in the infrared or optical range, the light momentum is nearly 0 and the effect of the
    // light's finite momentum is negligible in typical scatterings. In these cases, the only important charactestic of light is its energy.
    // Therefore specifying the correct momentum of light becomes a waste of time. Furthermore it can create problems since the momentum can be so small
    // that it can hit the limits of numerical precision.
    // The methods below provide an easy way of setting up good meshes and dispersions for optical photons
    // WARNING: these methods should not be used when the momentum of the photon is relevant like in laser cooling or high energy photons.
    
    // The method below returns a mesh suitable for the construction of photon bands. The photon momentum is set to be very close to k = (0,0,0)
    // The user must save the returned mesh into a variable that must live as long as the material does
    Mesh<NDim> photonMesh(int resolution) const;
    // The method below adds a photon band to the material
    // Several photon bands can be added using the same mesh (for instance in the case different photon bands are required to describe different frequency ranges
    void addPhotonBand(const std::string &t_id, Mesh<NDim> const& photonMesh, Real lowestEnergy, Real highestEnergy);
    // This prevents that a temporary mesh is used 
    void addPhotonBand(const std::string &t_id, Mesh<NDim> const&& photonMesh, Real lowestEnergy, Real highestEnergy) = delete;
    
    //********************************
    //* I/O
    //********************************
    void plot() const;
    void plot(const std::vector<std::string> &bandstoplot) const;
    
    template<typename Derived, typename... Args> void addSeriesSpaghettiPoints(const Eigen::MatrixBase<Derived>& point, const std::string& namePoint, Args... args)  requires (NDim==2);
    
    void spaghettiPlot() const;
    void spaghettiPlot(std::string const & bandstoplotWildcards) const;
    void spaghettiPlot(const std::vector<std::string>& bandstoplot) const;
    
        
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Technical details
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
    // The user should not add bands except with the methods provided by Material so some of StringMap's methods are made private
    // https://stackoverflow.com/questions/4908539/a-way-in-c-to-hide-a-specific-function
    // Remove the possibility of passing temporaries to the constructor
    Material(const std::string &t_matid, const Region<NDim>&& t_region) = delete;
    
//private:
    template<typename Derived, typename... Args> void addSeriesSpaghettiPointsRecursion(const Eigen::MatrixBase<Derived>& point, const std::string& namePoint, Args... args) requires (NDim==2) ;
    template<typename Derived> void addSeriesSpaghettiPointsRecursion(const Eigen::MatrixBase<Derived>& point, const std::string& namePoint) requires (NDim==2);
    
    // Data for spaghetti plot
    std::vector<VectorPoint<NDim>>             seriesSpaghettiPoints;
    std::vector<std::vector<std::string>>      seriesSpaghettiPointsNames;

};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Friend methods
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~
Plotter3D& operator << (Plotter3D& plotter, const Material<2>& region);
template<int NDim> std::ostream &operator<<(std::ostream &os, Material<NDim> const& material);




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Technical details
//~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~

template<int NDim> template<typename Derived> void Material<NDim>::addSeriesSpaghettiPointsRecursion(const Eigen::MatrixBase<Derived>& point, const std::string& namePoint) requires (NDim==2) {
    seriesSpaghettiPoints.back().conservativeResize(Eigen::NoChange, seriesSpaghettiPoints.back().cols()+1);
    seriesSpaghettiPoints.back().col(seriesSpaghettiPoints.back().cols()-1) = point;
    seriesSpaghettiPointsNames.back().emplace_back(namePoint);
}

template<int NDim> template<typename Derived, typename... Args> void Material<NDim>::addSeriesSpaghettiPointsRecursion(const Eigen::MatrixBase<Derived>& point, const std::string& namePoint, Args... args) requires (NDim==2) {
    seriesSpaghettiPoints.back().conservativeResize(Eigen::NoChange, seriesSpaghettiPoints.back().cols()+1);
    seriesSpaghettiPoints.back().col(seriesSpaghettiPoints.back().cols()-1) = point;
    seriesSpaghettiPointsNames.back().emplace_back(namePoint);
    addSeriesSpaghettiPointsRecursion(args...);
}

template<int NDim> template<typename Derived, typename... Args> void Material<NDim>::addSeriesSpaghettiPoints(const Eigen::MatrixBase<Derived>& point, const std::string& namePoint, Args... args) requires (NDim==2) {
    VectorPoint<NDim> points;
    points.resize(Eigen::NoChange, 1);
    points.col(0) = point;
    seriesSpaghettiPoints.emplace_back(points);
    std::vector<std::string> names({namePoint});
    seriesSpaghettiPointsNames.emplace_back(names);
    addSeriesSpaghettiPointsRecursion(args...);
}



};   // namespace Tortoise

#endif /* Material_hpp */
