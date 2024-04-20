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
//  Mesh.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/7/20.
//

#include <Geometry/Structured/DirectSpace/Mesh.hpp>
#include <Generics/Utilities/RandomNGenerator.hpp>

#include <cassert>
#include <fstream>

using Tortoise::Utilities::randGen;
using namespace Tortoise::GeometryCore;

namespace Tortoise {

// Output for EigenMatrix types
inline std::string to_string(const Eigen::MatrixXd& mat){
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

void                                              addNextElementVertex(const std::vector<std::vector<Real>>& CurrentIncompleteVertices, std::vector<std::vector<std::vector<Real>>>& ListElementVertices, const int dim) {
    if (CurrentIncompleteVertices.size()==dim+1) {
        ListElementVertices.push_back(CurrentIncompleteVertices);
    }
    else {
        for (int i=0; i<dim; ++i){
            if((CurrentIncompleteVertices.back())[i]<0.001){
                std::vector<std::vector<Real>> newlistVertices(CurrentIncompleteVertices);
                newlistVertices.push_back(CurrentIncompleteVertices.back());
                (newlistVertices.back())[i]=1.0;
                addNextElementVertex(newlistVertices,ListElementVertices,dim);
            }
        }
    }
};
template<int NDim>  VectorElement<NDim>           calculatevertElemInRefSec(){
    VectorElement<NDim> toreturn;
    std::vector<std::vector<std::vector<Real>>> ListElementVertices;
    std::vector<Real> start(NDim,0.0);
    addNextElementVertex({start}, ListElementVertices, NDim);
    for (int i=0; i<ListElementVertices.size(); ++i){
        Element<NDim> temporary;
        for (int col=0; col<ListElementVertices[i].size(); ++col){
            for (int row=0; row<ListElementVertices[i][col].size(); ++row){
                temporary(row,col)= ListElementVertices[i][col][row];
            }
        }
        toreturn.push_back(temporary);
    }
    return toreturn;
}

template<int NDim> const VectorElement<NDim>                Mesh<NDim>::vertElemInRefSec        = calculatevertElemInRefSec<NDim>();


//=======================================================
// Constructors
//===================
// *** Functions to calculate the const attributes in the constructor

template<int NDim> LinTransform<NDim> calculate_massMatrix(const double elemVolume){
    LinTransform<NDim> toreturn;
    LinearForm<NDim> basisFunction0, basisFunction1;
    for (int i=0; i<NDim+1; ++i){
        basisFunction0 = LinearForm<NDim>::Zero();
        basisFunction0(i) = 1.0;
        for (int j=0; j<NDim+1; ++j){
            basisFunction1 = LinearForm<NDim>::Zero();
            basisFunction1(j) = 1.0;
            toreturn(i,j) = elemVolume*((apply(basisFunction0,gaussPoints<NDim,2>).array()*apply(basisFunction1,gaussPoints<NDim,2>).array()).matrix()*(gaussWeigths<NDim,2>.asDiagonal())).sum();
        }
    }
    return toreturn;
}

template<int NDim> LinTransform<NDim> calculate_inverseMassMatrix(const double elemVolume){
    LinTransform<NDim> toreturn;
    LinearForm<NDim> basisFunction0, basisFunction1;
    for (int i=0; i<NDim+1; ++i){
        basisFunction0 = LinearForm<NDim>::Zero();
        basisFunction0(i) = 1.0;
        for (int j=0; j<NDim+1; ++j){
            basisFunction1 = LinearForm<NDim>::Zero();
            basisFunction1(j) = 1.0;
            toreturn(i,j) = elemVolume*((apply(basisFunction0,gaussPoints<NDim,2>).array()*apply(basisFunction1,gaussPoints<NDim,2>).array()).matrix()*(gaussWeigths<NDim,2>.asDiagonal())).sum();
        }
    }
    return toreturn.inverse();
}

template<int NDim> VecOfDataVector calculate_verticesCoord(const Point<NDim>& t_OriginCoord, const ArrayPoint<NDim,NDim>& t_gVecfrac, const CartIndex<NDim>& t_nSecSplits);

template<int NDim> VectorLinearForm<NDim> calculate_basisFuncLF(const CartIndex<NDim>& t_nSecSplits){
    VectorLinearForm<NDim> toreturn;
    toreturn.resize(Mesh<NDim>::nVertElement * Mesh<NDim>::nElementSection * t_nSecSplits.prod(),Eigen::NoChange);
    for (int i=0; i<Mesh<NDim>::nElementSection * t_nSecSplits.prod(); ++i){
        toreturn.template block<Mesh<NDim>::nVertElement,Mesh<NDim>::nVertElement>(Mesh<NDim>::nVertElement*i,0).setIdentity();
    }
    return toreturn;
}

template<int NDim> std::vector<VectorLinearForm<NDim>> calculate_basisFuncDerivLF(const CartIndex<NDim>& t_nSecSplits, const VecOfDataVector& verticesCoord);

template<> std::vector<VectorLinearForm<1>> calculate_basisFuncDerivLF<1>(const CartIndex<1>& t_nSecSplits, const VecOfDataVector& verticesCoord){
    std::vector<VectorLinearForm<1>> toreturn;
    toreturn.emplace_back(VectorLinearForm<1>::Zero(Mesh<1>::nVertElement * Mesh<1>::nElementSection * t_nSecSplits.prod(),Mesh<1>::nVertElement));
    for (int elem=0; elem < Mesh<1>::nElementSection * t_nSecSplits.prod(); ++elem){  // loop over elements
        toreturn[0](elem * Mesh<1>::nVertElement + 1, 0) = 1. / (verticesCoord[0](1)-verticesCoord[0](0));
    }
    return toreturn;
}

template<> std::vector<VectorLinearForm<2>> calculate_basisFuncDerivLF<2>(const CartIndex<2>& t_nSecSplits, const VecOfDataVector& verticesCoord){
    std::vector<VectorLinearForm<2>> toreturn;
    for (int dim=0; dim<2; ++dim){  // Loop over derivation direction
        toreturn.emplace_back(VectorLinearForm<2>::Zero(Mesh<2>::nVertElement * Mesh<2>::nElementSection * t_nSecSplits.prod(),Mesh<2>::nVertElement));
    }
    for (int elem=0; elem < Mesh<2>::nElementSection * t_nSecSplits.prod(); ++elem){  // loop over elements
        const Real v0X = verticesCoord[0](elem * Mesh<2>::nVertElement),   v0Y = verticesCoord[1](elem * Mesh<2>::nVertElement);
        const Real v1X = verticesCoord[0](elem * Mesh<2>::nVertElement+1), v1Y = verticesCoord[1](elem * Mesh<2>::nVertElement+1);
        const Real v2X = verticesCoord[0](elem * Mesh<2>::nVertElement+2), v2Y = verticesCoord[1](elem * Mesh<2>::nVertElement+2);
        
        toreturn[0](elem * Mesh<2>::nVertElement + 1, 0) = (v0Y-v2Y)/( v0Y*(v1X-v2X) + v1Y*(v2X-v0X) + v2Y*(v0X-v1X));
        toreturn[0](elem * Mesh<2>::nVertElement + 2, 0) = - (v0Y-v1Y)/( v0Y*(v1X-v2X) + v1Y*(v2X-v0X) + v2Y*(v0X-v1X));
        toreturn[1](elem * Mesh<2>::nVertElement + 1, 0) = (v0X-v2X)/( v0X*(v1Y-v2Y) + v1X*(v2Y-v0Y) + v2X*(v0Y-v1Y));
        toreturn[1](elem * Mesh<2>::nVertElement + 2, 0) = - (v0X-v1X)/( v0X*(v1Y-v2Y) + v1X*(v2Y-v0Y) + v2X*(v0Y-v1Y));
    }
    return toreturn;
}

template<> std::vector<VectorLinearForm<3>> calculate_basisFuncDerivLF<3>(const CartIndex<3>& t_nSecSplits, const VecOfDataVector& verticesCoord){
    std::vector<VectorLinearForm<3>> toreturn;
    for (int dim=0; dim<3; ++dim){  // Loop over derivation direction
        toreturn.emplace_back(VectorLinearForm<3>::Zero(Mesh<3>::nVertElement * Mesh<3>::nElementSection * t_nSecSplits.prod(),Mesh<3>::nVertElement));
    }
    for (int elem=0; elem < Mesh<3>::nElementSection * t_nSecSplits.prod(); ++elem){  // loop over elements
        const Real v0X = verticesCoord[0](elem * Mesh<3>::nVertElement),   v0Y = verticesCoord[1](elem * Mesh<3>::nVertElement),   v0Z = verticesCoord[2](elem * Mesh<3>::nVertElement);
        const Real v1X = verticesCoord[0](elem * Mesh<3>::nVertElement+1), v1Y = verticesCoord[1](elem * Mesh<3>::nVertElement+1), v1Z = verticesCoord[2](elem * Mesh<3>::nVertElement+1);
        const Real v2X = verticesCoord[0](elem * Mesh<3>::nVertElement+2), v2Y = verticesCoord[1](elem * Mesh<3>::nVertElement+2), v2Z = verticesCoord[2](elem * Mesh<3>::nVertElement+2);
        const Real v3X = verticesCoord[0](elem * Mesh<3>::nVertElement+3), v3Y = verticesCoord[1](elem * Mesh<3>::nVertElement+3), v3Z = verticesCoord[2](elem * Mesh<3>::nVertElement+3);
        
        // The expressions below have been constructed using the following python program (and modifications of it)
        // import sympy as sym
        // v0x = sym.Symbol('v0X'); v0y = sym.Symbol('v0Y'); v0z = sym.Symbol('v0Z')
        // v1x = sym.Symbol('v1X'); v1y = sym.Symbol('v1Y'); v1z = sym.Symbol('v1Z')
        // v2x = sym.Symbol('v2X'); v2y = sym.Symbol('v2Y'); v2z = sym.Symbol('v2Z')
        // v3x = sym.Symbol('v3X'); v3y = sym.Symbol('v3Y'); v3z = sym.Symbol('v3Z')
        // a = sym.Symbol('a'); b = sym.Symbol('b'); c = sym.Symbol('c'); d = sym.Symbol('d');
        // solution = sym.solve((a + b * v0x + c * v0y + d * v0z - 0,
        //                       a + b * v1x + c * v1y + d * v1z - 0,
        //                       a + b * v2x + c * v2y + d * v2z - 0,
        //                       a + b * v3x + c * v3y + d * v3z - 1), (a, b, c, d))
        // print(solution[d])
        
        toreturn[0](elem * Mesh<3>::nVertElement + 1, 0) = (-v0Y*v2Z + v0Y*v3Z + v0Z*v2Y - v0Z*v3Y - v2Y*v3Z + v2Z*v3Y)/(v0X*v1Y*v2Z - v0X*v1Y*v3Z - v0X*v1Z*v2Y + v0X*v1Z*v3Y + v0X*v2Y*v3Z - v0X*v2Z*v3Y - v0Y*v1X*v2Z + v0Y*v1X*v3Z + v0Y*v1Z*v2X - v0Y*v1Z*v3X - v0Y*v2X*v3Z + v0Y*v2Z*v3X + v0Z*v1X*v2Y - v0Z*v1X*v3Y - v0Z*v1Y*v2X + v0Z*v1Y*v3X + v0Z*v2X*v3Y - v0Z*v2Y*v3X - v1X*v2Y*v3Z + v1X*v2Z*v3Y + v1Y*v2X*v3Z - v1Y*v2Z*v3X - v1Z*v2X*v3Y + v1Z*v2Y*v3X);
        toreturn[0](elem * Mesh<3>::nVertElement + 2, 0) = (v0Y*v1Z - v0Y*v3Z - v0Z*v1Y + v0Z*v3Y + v1Y*v3Z - v1Z*v3Y)/(v0X*v1Y*v2Z - v0X*v1Y*v3Z - v0X*v1Z*v2Y + v0X*v1Z*v3Y + v0X*v2Y*v3Z - v0X*v2Z*v3Y - v0Y*v1X*v2Z + v0Y*v1X*v3Z + v0Y*v1Z*v2X - v0Y*v1Z*v3X - v0Y*v2X*v3Z + v0Y*v2Z*v3X + v0Z*v1X*v2Y - v0Z*v1X*v3Y - v0Z*v1Y*v2X + v0Z*v1Y*v3X + v0Z*v2X*v3Y - v0Z*v2Y*v3X - v1X*v2Y*v3Z + v1X*v2Z*v3Y + v1Y*v2X*v3Z - v1Y*v2Z*v3X - v1Z*v2X*v3Y + v1Z*v2Y*v3X);
        toreturn[0](elem * Mesh<3>::nVertElement + 3, 0) = (-v0Y*v1Z + v0Y*v2Z + v0Z*v1Y - v0Z*v2Y - v1Y*v2Z + v1Z*v2Y)/(v0X*v1Y*v2Z - v0X*v1Y*v3Z - v0X*v1Z*v2Y + v0X*v1Z*v3Y + v0X*v2Y*v3Z - v0X*v2Z*v3Y - v0Y*v1X*v2Z + v0Y*v1X*v3Z + v0Y*v1Z*v2X - v0Y*v1Z*v3X - v0Y*v2X*v3Z + v0Y*v2Z*v3X + v0Z*v1X*v2Y - v0Z*v1X*v3Y - v0Z*v1Y*v2X + v0Z*v1Y*v3X + v0Z*v2X*v3Y - v0Z*v2Y*v3X - v1X*v2Y*v3Z + v1X*v2Z*v3Y + v1Y*v2X*v3Z - v1Y*v2Z*v3X - v1Z*v2X*v3Y + v1Z*v2Y*v3X);
        
        toreturn[1](elem * Mesh<3>::nVertElement + 1, 0) = (v0X*v2Z - v0X*v3Z - v0Z*v2X + v0Z*v3X + v2X*v3Z - v2Z*v3X)/(v0X*v1Y*v2Z - v0X*v1Y*v3Z - v0X*v1Z*v2Y + v0X*v1Z*v3Y + v0X*v2Y*v3Z - v0X*v2Z*v3Y - v0Y*v1X*v2Z + v0Y*v1X*v3Z + v0Y*v1Z*v2X - v0Y*v1Z*v3X - v0Y*v2X*v3Z + v0Y*v2Z*v3X + v0Z*v1X*v2Y - v0Z*v1X*v3Y - v0Z*v1Y*v2X + v0Z*v1Y*v3X + v0Z*v2X*v3Y - v0Z*v2Y*v3X - v1X*v2Y*v3Z + v1X*v2Z*v3Y + v1Y*v2X*v3Z - v1Y*v2Z*v3X - v1Z*v2X*v3Y + v1Z*v2Y*v3X);
        toreturn[1](elem * Mesh<3>::nVertElement + 2, 0) = (-v0X*v1Z + v0X*v3Z + v0Z*v1X - v0Z*v3X - v1X*v3Z + v1Z*v3X)/(v0X*v1Y*v2Z - v0X*v1Y*v3Z - v0X*v1Z*v2Y + v0X*v1Z*v3Y + v0X*v2Y*v3Z - v0X*v2Z*v3Y - v0Y*v1X*v2Z + v0Y*v1X*v3Z + v0Y*v1Z*v2X - v0Y*v1Z*v3X - v0Y*v2X*v3Z + v0Y*v2Z*v3X + v0Z*v1X*v2Y - v0Z*v1X*v3Y - v0Z*v1Y*v2X + v0Z*v1Y*v3X + v0Z*v2X*v3Y - v0Z*v2Y*v3X - v1X*v2Y*v3Z + v1X*v2Z*v3Y + v1Y*v2X*v3Z - v1Y*v2Z*v3X - v1Z*v2X*v3Y + v1Z*v2Y*v3X);
        toreturn[1](elem * Mesh<3>::nVertElement + 3, 0) = (v0X*v1Z - v0X*v2Z - v0Z*v1X + v0Z*v2X + v1X*v2Z - v1Z*v2X)/(v0X*v1Y*v2Z - v0X*v1Y*v3Z - v0X*v1Z*v2Y + v0X*v1Z*v3Y + v0X*v2Y*v3Z - v0X*v2Z*v3Y - v0Y*v1X*v2Z + v0Y*v1X*v3Z + v0Y*v1Z*v2X - v0Y*v1Z*v3X - v0Y*v2X*v3Z + v0Y*v2Z*v3X + v0Z*v1X*v2Y - v0Z*v1X*v3Y - v0Z*v1Y*v2X + v0Z*v1Y*v3X + v0Z*v2X*v3Y - v0Z*v2Y*v3X - v1X*v2Y*v3Z + v1X*v2Z*v3Y + v1Y*v2X*v3Z - v1Y*v2Z*v3X - v1Z*v2X*v3Y + v1Z*v2Y*v3X);
        
        toreturn[2](elem * Mesh<3>::nVertElement + 1, 0) = (-v0X*v2Y + v0X*v3Y + v0Y*v2X - v0Y*v3X - v2X*v3Y + v2Y*v3X)/(v0X*v1Y*v2Z - v0X*v1Y*v3Z - v0X*v1Z*v2Y + v0X*v1Z*v3Y + v0X*v2Y*v3Z - v0X*v2Z*v3Y - v0Y*v1X*v2Z + v0Y*v1X*v3Z + v0Y*v1Z*v2X - v0Y*v1Z*v3X - v0Y*v2X*v3Z + v0Y*v2Z*v3X + v0Z*v1X*v2Y - v0Z*v1X*v3Y - v0Z*v1Y*v2X + v0Z*v1Y*v3X + v0Z*v2X*v3Y - v0Z*v2Y*v3X - v1X*v2Y*v3Z + v1X*v2Z*v3Y + v1Y*v2X*v3Z - v1Y*v2Z*v3X - v1Z*v2X*v3Y + v1Z*v2Y*v3X);
        toreturn[2](elem * Mesh<3>::nVertElement + 2, 0) = (v0X*v1Y - v0X*v3Y - v0Y*v1X + v0Y*v3X + v1X*v3Y - v1Y*v3X)/(v0X*v1Y*v2Z - v0X*v1Y*v3Z - v0X*v1Z*v2Y + v0X*v1Z*v3Y + v0X*v2Y*v3Z - v0X*v2Z*v3Y - v0Y*v1X*v2Z + v0Y*v1X*v3Z + v0Y*v1Z*v2X - v0Y*v1Z*v3X - v0Y*v2X*v3Z + v0Y*v2Z*v3X + v0Z*v1X*v2Y - v0Z*v1X*v3Y - v0Z*v1Y*v2X + v0Z*v1Y*v3X + v0Z*v2X*v3Y - v0Z*v2Y*v3X - v1X*v2Y*v3Z + v1X*v2Z*v3Y + v1Y*v2X*v3Z - v1Y*v2Z*v3X - v1Z*v2X*v3Y + v1Z*v2Y*v3X);
        toreturn[2](elem * Mesh<3>::nVertElement + 3, 0) = (-v0X*v1Y + v0X*v2Y + v0Y*v1X - v0Y*v2X - v1X*v2Y + v1Y*v2X)/(v0X*v1Y*v2Z - v0X*v1Y*v3Z - v0X*v1Z*v2Y + v0X*v1Z*v3Y + v0X*v2Y*v3Z - v0X*v2Z*v3Y - v0Y*v1X*v2Z + v0Y*v1X*v3Z + v0Y*v1Z*v2X - v0Y*v1Z*v3X - v0Y*v2X*v3Z + v0Y*v2Z*v3X + v0Z*v1X*v2Y - v0Z*v1X*v3Y - v0Z*v1Y*v2X + v0Z*v1Y*v3X + v0Z*v2X*v3Y - v0Z*v2Y*v3X - v1X*v2Y*v3Z + v1X*v2Z*v3Y + v1Y*v2X*v3Z - v1Y*v2Z*v3X - v1Z*v2X*v3Y + v1Z*v2Y*v3X);
    }
    return toreturn;
}

template<int NDim> VecOfDataVector calculate_basisFuncDerivLF0thOrd(const std::vector<VectorLinearForm<NDim>> basisFuncDerivLF){
    VecOfDataVector toreturn;
    for (int dim=0; dim<NDim; ++dim){
        toreturn.emplace_back(basisFuncDerivLF[dim].col(0));
    }
    return toreturn;
}


template<> VecOfDataVector calculate_verticesCoord<1>(const Point<1>& t_OriginCoord, const ArrayPoint<1,1>& t_dGVec, const CartIndex<1>& t_nSecSplits){
    VecOfDataVector toreturn;
    toreturn.emplace_back(DataVector::Zero(2*t_nSecSplits(0)));
    for (int i = 0; i< t_nSecSplits[0]; ++i){
        toreturn[0](2*i) = t_OriginCoord(0) + static_cast<Real>(i) * t_dGVec(0);
        toreturn[0](2*i+1) = t_OriginCoord(0) + static_cast<Real>(i+1) * t_dGVec(0);
    }
    return toreturn;
}
template<> VecOfDataVector calculate_verticesCoord<2>(const Point<2>& t_OriginCoord, const ArrayPoint<2,2>& t_dGVec, const CartIndex<2>& t_nSecSplits){
    VecOfDataVector toreturn;
    const int  dataVectorDim = t_nSecSplits.prod() * Mesh<2>::nElementSection * Mesh<2>::nVertElement;
    toreturn.emplace_back(DataVector::Zero(1,dataVectorDim));   // X coordinate
    toreturn.emplace_back(DataVector::Zero(1,dataVectorDim));   // Y coordinate
    for (int secY = 0; secY< t_nSecSplits(1); ++secY){
        for (int secX = 0; secX< t_nSecSplits(0); ++secX){
            for (int elemInSec = 0; elemInSec < Mesh<2>::nElementSection; ++elemInSec){
                for (int vertexinElem = 0; vertexinElem < Mesh<2>::nVertElement; ++vertexinElem){
                    const Point<2>  sectcartindouble({static_cast<Real>(secX),static_cast<Real>(secY)});
                    const Point<2>  currvert = t_OriginCoord + t_dGVec * (sectcartindouble + Mesh<2>::vertElemInRefSec[elemInSec].col(vertexinElem));
                    const int       curreIndex = vertexinElem + (elemInSec + (secX + secY * t_nSecSplits(0)) * Mesh<2>::nElementSection) * Mesh<2>::nVertElement;
                    toreturn[0](curreIndex) = currvert(0);
                    toreturn[1](curreIndex) = currvert(1);
                }
            }
        }
    }
    return toreturn;
}
template<> VecOfDataVector calculate_verticesCoord<3>(const Point<3>& t_OriginCoord, const ArrayPoint<3,3>& t_dGVec, const CartIndex<3>& t_nSecSplits){
    VecOfDataVector toreturn;
    const int  dataVectorDim = t_nSecSplits.prod() * Mesh<3>::nElementSection * Mesh<3>::nVertElement;
    toreturn.emplace_back(DataVector::Zero(1,dataVectorDim));   // X coordinate
    toreturn.emplace_back(DataVector::Zero(1,dataVectorDim));   // Y coordinate
    toreturn.emplace_back(DataVector::Zero(1,dataVectorDim));   // Z coordinate
    for (int secZ = 0; secZ< t_nSecSplits(2); ++secZ){
        for (int secY = 0; secY< t_nSecSplits(1); ++secY){
            for (int secX = 0; secX< t_nSecSplits(0); ++secX){
                for (int elemInSec = 0; elemInSec < Mesh<3>::nElementSection; ++elemInSec){
                    for (int vertexinElem = 0; vertexinElem < Mesh<3>::nVertElement; ++vertexinElem){
                        const Point<3>  sectcartindouble({static_cast<Real>(secX),static_cast<Real>(secY),static_cast<Real>(secZ)});
                        const Point<3>  currvert = t_OriginCoord + t_dGVec * (sectcartindouble + Mesh<3>::vertElemInRefSec[elemInSec].col(vertexinElem));
                        const int       curreIndex = vertexinElem + (elemInSec + (secX + (secY + secZ * t_nSecSplits(1)) * t_nSecSplits(0)) * Mesh<3>::nElementSection) * Mesh<3>::nVertElement;
                        toreturn[0](curreIndex) = currvert(0);
                        toreturn[1](curreIndex) = currvert(1);
                        toreturn[2](curreIndex) = currvert(2);
                    }
                }
            }
        }
    }
    return toreturn;
}

// Constructor
template<int NDim> Mesh<NDim>::Mesh(Region<NDim> const & t_region, const Point<NDim>& t_OriginCoord, const Point<NDim>& t_gVecfrac, const CartIndex<NDim>& t_nSecSplits):
    region(&t_region),
    origin(t_region.origin+t_region.gVec * t_OriginCoord),
    gVec(t_region.gVec * t_gVecfrac.asDiagonal()),
    invgVec(gVec.inverse()),
    relativeOrigin(t_OriginCoord),
    relativegVec(t_gVecfrac),
    relativedGVec(t_gVecfrac.array() * t_nSecSplits.template cast<Real>().cwiseInverse().array()),
    nSecSplits(t_nSecSplits),
    numberSections(t_nSecSplits.prod()),
    numberElements(t_nSecSplits.prod()*nElementSection),
    dataVectorDim(t_nSecSplits.prod()*nElementSection*nVertElement),
    dGVec(gVec.array().rowwise() * nSecSplits.template cast<Real>().cwiseInverse().transpose().array()),
    invdGVec(dGVec.inverse()),
    elemVolume(std::fabs(dGVec.determinant()/nElementSection)),
    massMatrix(calculate_massMatrix<NDim>(elemVolume)),
    invMassMatrix(calculate_inverseMassMatrix<NDim>(elemVolume)),
    verticesCoord(calculate_verticesCoord<NDim>(t_OriginCoord, dGVec, t_nSecSplits)) ,
    basisFuncLF(calculate_basisFuncLF<NDim>(t_nSecSplits)) ,
    basisFuncDerivLF(calculate_basisFuncDerivLF<NDim>(t_nSecSplits, verticesCoord)) ,
    basisFuncDerivLF0thOrd(calculate_basisFuncDerivLF0thOrd<NDim>(basisFuncDerivLF))
{
// ******* Assertions
        assert( (t_nSecSplits.array() > 0).all() );
        assert( (t_gVecfrac.array() > 0.).all() );
};


template<int NDim> Mesh<NDim>::Mesh(Region<NDim> const & t_region, const Real t_relativeOrigin, const Real t_relativegVec, const int t_nSecSplits) requires (NDim == 1) :
    region(&t_region) ,
    origin(t_region.origin + t_region.gVec * Point<NDim>::Constant(t_relativeOrigin) ) ,
    gVec(t_region.gVec * Point<NDim>::Constant(t_relativegVec).asDiagonal()) ,
    invgVec(gVec.inverse()),
    relativeOrigin(Point<NDim>::Constant(t_relativeOrigin)) ,
    relativegVec(Point<NDim>::Constant(t_relativegVec)) ,
    relativedGVec(Point<NDim>::Constant(t_relativegVec).array() * CartIndex<NDim>::Constant(t_nSecSplits).template cast<Real>().cwiseInverse().array()),
    nSecSplits(CartIndex<NDim>::Constant(t_nSecSplits)),
    numberSections(t_nSecSplits) ,
    numberElements(t_nSecSplits*nElementSection) ,
    dataVectorDim(t_nSecSplits*nElementSection*nVertElement) ,
    dGVec(gVec.array().rowwise() * nSecSplits.template cast<Real>().cwiseInverse().transpose().array()) ,
    invdGVec(dGVec.inverse()),
    elemVolume(std::fabs(dGVec.determinant()/nElementSection)) ,
    massMatrix(calculate_massMatrix<NDim>(elemVolume)) ,
    invMassMatrix(calculate_inverseMassMatrix<NDim>(elemVolume)),
    verticesCoord(calculate_verticesCoord<NDim>(Point<NDim>::Constant(t_relativeOrigin), dGVec, CartIndex<NDim>::Constant(t_nSecSplits))) ,
    basisFuncLF(calculate_basisFuncLF<NDim>(CartIndex<NDim>::Constant(t_nSecSplits))),
    basisFuncDerivLF(calculate_basisFuncDerivLF<NDim>(CartIndex<NDim>::Constant(t_nSecSplits), verticesCoord)) ,
    basisFuncDerivLF0thOrd(calculate_basisFuncDerivLF0thOrd<NDim>(basisFuncDerivLF))
{
// ******* Assertions
    assert( t_nSecSplits > 0 );
    assert( t_relativegVec > 0. );
};

template<int NDim> Mesh<NDim>::Mesh(const Region<NDim>* t_region, const Point<NDim>& t_OriginCoord, const Point<NDim>& t_gVecfrac, const CartIndex<NDim>& t_nSecSplits):
    region(t_region),
    origin(t_region->origin+t_region->gVec * t_OriginCoord),
    gVec(t_region->gVec * t_gVecfrac.asDiagonal()),
    invgVec(gVec.inverse()),
    relativeOrigin(t_OriginCoord),
    relativegVec(t_gVecfrac),
    relativedGVec(t_gVecfrac.array() * t_nSecSplits.template cast<Real>().cwiseInverse().array()),
    nSecSplits(t_nSecSplits),
    numberSections(t_nSecSplits.prod()),
    numberElements(t_nSecSplits.prod()*nElementSection),
    dataVectorDim(t_nSecSplits.prod()*nElementSection*nVertElement),
    dGVec(gVec.array().rowwise() * nSecSplits.template cast<Real>().cwiseInverse().transpose().array()),
    invdGVec(dGVec.inverse()),
    elemVolume(std::fabs(dGVec.determinant()/nElementSection)),
    massMatrix(calculate_massMatrix<NDim>(elemVolume)),
    invMassMatrix(calculate_inverseMassMatrix<NDim>(elemVolume)),
    verticesCoord(calculate_verticesCoord<NDim>(t_OriginCoord, dGVec, t_nSecSplits)) ,
    basisFuncLF(calculate_basisFuncLF<NDim>(t_nSecSplits)) ,
    basisFuncDerivLF(calculate_basisFuncDerivLF<NDim>(t_nSecSplits, verticesCoord)) ,
    basisFuncDerivLF0thOrd(calculate_basisFuncDerivLF0thOrd<NDim>(basisFuncDerivLF))
{
// ******* Assertions
    assert( (t_nSecSplits.array() > 0).all() );
    assert( (t_gVecfrac.array() > 0.).all() );
};

template<int NDim>  Mesh<NDim>::Mesh(const Region<NDim>* t_region, const Real t_relativeOrigin, const Real t_relativegVec, const int t_nSecSplits) requires (NDim == 1) :
region(t_region),
origin(t_region->origin + t_region->gVec * Point<NDim>::Constant(t_relativeOrigin) ) ,
gVec(t_region->gVec * Point<NDim>::Constant(t_relativegVec).asDiagonal()) ,
invgVec(gVec.inverse()),
relativeOrigin(Point<NDim>::Constant(t_relativeOrigin)) ,
relativegVec(Point<NDim>::Constant(t_relativegVec)) ,
relativedGVec(Point<NDim>::Constant(t_relativegVec).array() * CartIndex<NDim>::Constant(t_nSecSplits).template cast<Real>().cwiseInverse().array()),
nSecSplits(CartIndex<NDim>::Constant(t_nSecSplits)),
numberSections(t_nSecSplits) ,
numberElements(t_nSecSplits*nElementSection) ,
dataVectorDim(t_nSecSplits*nElementSection*nVertElement) ,
dGVec(gVec.array().rowwise() * nSecSplits.template cast<Real>().cwiseInverse().transpose().array()) ,
invdGVec(dGVec.inverse()),
elemVolume(std::fabs(dGVec.determinant()/nElementSection)) ,
massMatrix(calculate_massMatrix<NDim>(elemVolume)) ,
invMassMatrix(calculate_inverseMassMatrix<NDim>(elemVolume)),
verticesCoord(calculate_verticesCoord<NDim>(Point<NDim>::Constant(t_relativeOrigin), dGVec, CartIndex<NDim>::Constant(t_nSecSplits))) ,
basisFuncLF(calculate_basisFuncLF<NDim>(CartIndex<NDim>::Constant(t_nSecSplits))),
basisFuncDerivLF(calculate_basisFuncDerivLF<NDim>(CartIndex<NDim>::Constant(t_nSecSplits), verticesCoord)) ,
basisFuncDerivLF0thOrd(calculate_basisFuncDerivLF0thOrd<NDim>(basisFuncDerivLF))
{
// ******* Assertions
    assert( t_nSecSplits > 0 );
    assert( t_relativegVec > 0. );
};

//********************************
//* Mesh structure
//********************************
//=======================================================
// Mesh Iterators and accessor
//===================
// Indices conversions
template<> CartIndex<1> Mesh<1>::sectionID(const int t_elementID) const{
    assert(t_elementID>=0 && t_elementID<numberElements);
    return CartIndex<1>(t_elementID);
}
template<> CartIndex<2> Mesh<2>::sectionID(const int t_elementID) const{
    assert(t_elementID>=0 && t_elementID<numberElements);
    return CartIndex<2>({(t_elementID/2)%nSecSplits[0],(t_elementID/2)/nSecSplits[0]});
}
template<> CartIndex<3> Mesh<3>::sectionID(const int t_elementID) const{
    assert(t_elementID>=0 && t_elementID<numberElements);
    return CartIndex<3>({(t_elementID/6)%nSecSplits[0], (t_elementID/(6*nSecSplits[0]))%nSecSplits[1], (t_elementID/(6*nSecSplits[0]*nSecSplits[1]))});
}
template<> int Mesh<1>::sectionID(const CartIndex<1>& SectionCartID) const{
    assert( (SectionCartID.array()>=0).all()); assert( (SectionCartID.array()<nSecSplits.array()).all());
    return SectionCartID(0);
}
template<> int Mesh<2>::sectionID(const CartIndex<2>& SectionCartID) const{
    assert( (SectionCartID.array()>=0).all()); assert( (SectionCartID.array()<nSecSplits.array()).all());
    return SectionCartID(0) + SectionCartID(1)*nSecSplits(0);
}
template<> int Mesh<3>::sectionID(const CartIndex<3>& SectionCartID) const{
    assert( (SectionCartID.array()>=0).all()); assert( (SectionCartID.array()<nSecSplits.array()).all());
    return SectionCartID(0) + (SectionCartID(1) + SectionCartID(2)*nSecSplits(1))*nSecSplits(0);
}
template<> int Mesh<1>::elemID(const CartIndex<1>& SectionCartID, const int elemInSect) const{
    assert(elemInSect == 0); assert( (SectionCartID.array()>=0).all()); assert( (SectionCartID.array()<nSecSplits.array()).all());
    return SectionCartID(0);
}
template<> int Mesh<2>::elemID(const CartIndex<2>& SectionCartID, const int elemInSect) const{
    assert(elemInSect >= 0 && elemInSect < 2); assert( (SectionCartID.array()>=0).all());
    assert( (SectionCartID.array()<nSecSplits.array()).all());
    return 2*(SectionCartID(0) + SectionCartID(1)*nSecSplits(0)) + elemInSect;
}
template<> int Mesh<3>::elemID(const CartIndex<3>& SectionCartID, const int elemInSect) const{
    assert(elemInSect >= 0 && elemInSect < 6); assert( (SectionCartID.array()>=0).all());
    assert( (SectionCartID.array()<nSecSplits.array()).all());
    return 6*(SectionCartID(0) + (SectionCartID(1) + SectionCartID(2)*nSecSplits(1))*nSecSplits(0)) + elemInSect;
}
template<> int Mesh<1>::elemInSectionID(const int t_elementID) const {
    assert(t_elementID>=0 && t_elementID<numberElements);
    return 0;
};
template<> int Mesh<2>::elemInSectionID(const int t_elementID) const {
    assert(t_elementID>=0 && t_elementID<numberElements);
    return t_elementID%2;
};
template<> int Mesh<3>::elemInSectionID(const int t_elementID) const {
    assert(t_elementID>=0 && t_elementID<numberElements);
    return t_elementID%6;
};
// Mesh iterators
template<int NDim> SectionIterator<NDim> Mesh<NDim>::sectionIterator() const {
    return SectionIterator<NDim>(nSecSplits, CartIndex<NDim>::Zero(), nSecSplits);
}
template<int NDim> MeshElementIterator<NDim> Mesh<NDim>::elementIterator() const {
    return MeshElementIterator<NDim>(numberElements, 0);
}
template<int NDim> VertexInElemIterator<NDim> Mesh<NDim>::vertexInElemIterator() const {
    return VertexInElemIterator<NDim>();
}
// Random accessors
template<> CartIndex<1> Mesh<1>::randomSection() const {
    return CartIndex<1>(randGen.randDiscr<int>(0,nSecSplits(0)-1));
};
template<> CartIndex<2> Mesh<2>::randomSection() const {
    return CartIndex<2>({randGen.randDiscr<int>(0,nSecSplits(0)-1), randGen.randDiscr<int>(0,nSecSplits(1)-1)});
};
template<> CartIndex<3> Mesh<3>::randomSection() const {
    return CartIndex<3>({randGen.randDiscr<int>(0,nSecSplits(0)-1), randGen.randDiscr<int>(0,nSecSplits(1)-1), randGen.randDiscr<int>(0,nSecSplits(2)-1)});
};
template<int NDim> SectionIterator<NDim> Mesh<NDim>::randomSectionIterator() const {
    return SectionIterator<NDim>(nSecSplits, randomSection(), nSecSplits);
}
template<int NDim> int Mesh<NDim>::randomElement() const {
    return randGen.randDiscr<int>(0,numberElements-1);
};


//=======================================================
// Mesh geometry
//===================
// NOTE TO SELF: To be optimised
// Mesh Points from section ID
template<int NDim> Point<NDim> Mesh<NDim>::operator()(const std::array<int, NDim>& SectionCartID) const {
    Point<NDim> sectcartindouble;
    for (int i=0; i<NDim; ++i){sectcartindouble(i)=static_cast<Real>(SectionCartID[i]);}
    return origin + dGVec * sectcartindouble;
}
template<int NDim> Point<NDim> Mesh<NDim>::operator()(const SectionIterator<NDim>& t_sectionID) const {
    Point<NDim> toreturn(origin);
    for (int i=0; i<NDim; ++i){
        toreturn += dGVec.col(i)*static_cast<Real>(t_sectionID.currentcoord[i]);
    }
    return toreturn;
}
template<int NDim> Point<NDim> Mesh<NDim>::operator()(const CartIndex<NDim>& t_sectionID) const{
    Point<NDim> toreturn(origin);
    for (int i=0; i<NDim; ++i){
        toreturn += dGVec.col(i)*static_cast<Real>(t_sectionID[i]);
    }
    return toreturn;
}

// Vertex from element ID and vertex in element ID
template<int NDim> Point<NDim> Mesh<NDim>::operator()(const int t_elementID, const int t_nVertex) const {
    Point<NDim> sectcartindouble;
    CartIndex<NDim> currentSection = sectionID(t_elementID);
    int elemInSec = elemInSectionID(t_elementID);
    for (int i=0; i<NDim; ++i){sectcartindouble(i)=static_cast<Real>(currentSection(i));}
    return origin + dGVec * (sectcartindouble + vertElemInRefSec[elemInSec].col(t_nVertex));
}
template<int NDim> Point<NDim> Mesh<NDim>::operator()(const MeshElementIterator<NDim>& t_elementID, const VertexInElemIterator<NDim>& t_nVertex) const {
    return (*this)(t_elementID.currentElementID, t_nVertex.currentVertID);
}
template<int NDim> Point<NDim> Mesh<NDim>::operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const VertexInElemIterator<NDim>& t_nVertex) const {
    return (*this)( t_elementID.currentElementID , t_nVertex.currentVertID);
}

// Point from element ID and coordinates in reference element
template<int NDim> Point<NDim> Mesh<NDim>::operator()(const int t_elementID, const Point<NDim>& refCoord) const {
    Point<NDim> sectcartindouble;
    CartIndex<NDim> currentSection = sectionID(t_elementID);
    int elemInSec = elemInSectionID(t_elementID);
    for (int i=0; i<NDim; ++i){sectcartindouble(i)=static_cast<Real>(currentSection(i));}
    return origin + dGVec * (sectcartindouble + vertElemInRefSec[elemInSec].template block<NDim,NDim>(0,1) * refCoord);
}
template<int NDim> Point<NDim> Mesh<NDim>::operator()(const int t_elementID, const Real refCoord) const requires (NDim == 1){
    Point<NDim> sectcartindouble;
    CartIndex<NDim> currentSection = sectionID(t_elementID);
    int elemInSec = elemInSectionID(t_elementID);
    for (int i=0; i<NDim; ++i){sectcartindouble(i)=static_cast<Real>(currentSection(i));}
    return origin + dGVec * (sectcartindouble + vertElemInRefSec[elemInSec].template block<NDim,NDim>(0,1) * Point<NDim>::Constant(refCoord));
}
template<int NDim> Point<NDim> Mesh<NDim>::operator()(const MeshElementIterator<NDim>& t_elementID, const Point<NDim>& refCoord) const {
    return (*this)(t_elementID.currentElementID, refCoord);
}
template<int NDim> Point<NDim> Mesh<NDim>::operator()(const MeshElementIterator<NDim>& t_elementID, const Real refCoord) const requires (NDim == 1){
    return (*this)(t_elementID.currentElementID, refCoord);
}
template<int NDim> Point<NDim> Mesh<NDim>::operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const Point<NDim>& refCoord) const {
    return (*this)(t_elementID.currentElementID, refCoord);
}
template<int NDim> Point<NDim> Mesh<NDim>::operator()(const MeshSubsetElementIterator<NDim>& t_elementID, const Real refCoord) const requires (NDim == 1){
    return (*this)(t_elementID.currentElementID, refCoord);
}
// Baricentre of a given element
template<> Point<1> Mesh<1>::elemCentre(const int t_elementID) const {
    return (*this)(t_elementID,Point<1>(0.5));
}
template<> Point<2> Mesh<2>::elemCentre(const int t_elementID) const {
    return (*this)(t_elementID,Point<2>({0.333333333,0.333333333}));
}
template<> Point<3> Mesh<3>::elemCentre(const int t_elementID) const {
    return (*this)(t_elementID,Point<3>({0.25,0.25,0.25}));
}
template<int NDim> Point<NDim> Mesh<NDim>::elemCentre(const MeshElementIterator<NDim>& t_elementID) const {
    return elemCentre(t_elementID.currentElementID);
}
template<int NDim> Point<NDim> Mesh<NDim>::elemCentre(const MeshSubsetElementIterator<NDim>& t_elementID) const {
    return elemCentre(t_elementID.currentElementID);
}

//=======================================================
// Associated Meshes
//===================
template<> Mesh<1> Mesh<1>::convolutionMesh() const {
    Real  halfSizeConvMesh;
    if (((relativegVec(0, 0) > 0) && (region->gVec(0) < 0)) || ((relativegVec(0, 0) < 0) && (region->gVec(0) > 0))) { halfSizeConvMesh = -(gVec(0, 0)); }
    else { halfSizeConvMesh = (gVec(0, 0)); }
    
    int numberOfElemConvMesh        = 2 * numberElements;
    return Mesh<1>( region, Point<1>((- halfSizeConvMesh - region->origin(0))/region->gVec(0)), Point<1>(2.0*halfSizeConvMesh/region->gVec(0)), CartIndex<1>(numberOfElemConvMesh));
}
template<> Mesh<2> Mesh<2>::convolutionMesh() const {
    assert(false && "Function<2>::convolutionMesh not implemented");
    return *(this);
}
template<> Mesh<3> Mesh<3>::convolutionMesh() const {
    assert(false && "Function<3>::convolutionMesh not implemented");
    return *(this);
}

//=======================================================
// Adjacency Map
//===================
template<> int Mesh<1>::adjacentElement(const int elementID, const int face) const{
    assert(false);
    return 0;
}
template<> int Mesh<2>::adjacentElement(const int elementID, const int face) const{
    assert(face>=0 && face<3);
    if (elementID % 2) { // elementID is odd
        switch (face) {
            case 0:
                return (elementID + 2 * nSecSplits(0) - 1 + 2 * nSecSplits(0)*nSecSplits(1)) % (2 * nSecSplits(0)*nSecSplits(1));
            case 1:
                return elementID - 1;
            case 2:
                if ( (elementID-1)%(2*nSecSplits(0)) == 0) return elementID + 2 * nSecSplits(0) - 3;
                else return elementID - 3;
        }
    } else {             // elementID is even
        switch (face) {
            case 0:
                if ( (elementID+2)%(2*nSecSplits(0)) == 0) return elementID - 2 * nSecSplits(0) + 3;
                else return elementID + 3;
            case 1:
                return elementID + 1;
            case 2:
                return (elementID - 2 * nSecSplits(0) + 1 + 2 * nSecSplits(0)*nSecSplits(1)) % (2 * nSecSplits(0)*nSecSplits(1));
        }
    }
    return -1;
}
template<> int Mesh<3>::adjacentElement(const int elementID, const int face) const{
    assert(false);
    return 0;
}
template<> int Mesh<1>::adjacentFace(const int elementID, const int face) const{
    assert(false);
    return 0;
}
template<> int Mesh<2>::adjacentFace(const int elementID, const int face) const{
    assert(face>=0 && face<3);
    return 2 - face;
}
template<> int Mesh<3>::adjacentFace(const int elementID, const int face) const{
    assert(false);
    return 0;
}

template<> Point<1> Mesh<1>::faceNormal(const int elementID, const int face) const {
    assert(false);
    return Point<1>();
};
template<> Point<2> Mesh<2>::faceNormal(const int elementID, const int face) const {
    assert(face>=0 && face<3);
    Point<2> toreturn;
    if (elementID % 2) { // elementID is odd
        switch (face) {
            case 0: return Point<2>({dGVec(1,0),-dGVec(0,0)})/Point<2>({dGVec(1,0),-dGVec(0,0)}).norm() * ((Point<2>({dGVec(1,0),-dGVec(0,0)}).dot(dGVec.col(1))>0.)?1.:-1.);
            case 1: return Point<2>({dGVec(1,1)+dGVec(1,0),-dGVec(0,1)-dGVec(0,0)}) / Point<2>({dGVec(1,1)+dGVec(1,0),-dGVec(0,1)-dGVec(0,0)}).norm() * ((Point<2>({dGVec(1,1)+dGVec(1,0),-dGVec(0,1)-dGVec(0,0)}).dot(dGVec.col(0))>0.)?1.:-1.);
            case 2: return Point<2>({dGVec(1,1),-dGVec(0,1)})/Point<2>({dGVec(1,1),-dGVec(0,1)}).norm() * ((Point<2>({dGVec(1,1),-dGVec(0,1)}).dot(dGVec.col(0))>0.)?-1.:1.);
        }
    } else {             // elementID is even
        switch (face) {
            case 0: return Point<2>({dGVec(1,1),-dGVec(0,1)})/Point<2>({dGVec(1,1),-dGVec(0,1)}).norm() * ((Point<2>({dGVec(1,1),-dGVec(0,1)}).dot(dGVec.col(0))>0.)?1.:-1.);
            case 1: return Point<2>({dGVec(1,1)+dGVec(1,0),-dGVec(0,1)-dGVec(0,0)}) / Point<2>({dGVec(1,1)+dGVec(1,0),-dGVec(0,1)-dGVec(0,0)}).norm() * ((Point<2>({dGVec(1,1)+dGVec(1,0),-dGVec(0,1)-dGVec(0,0)}).dot(dGVec.col(0))>0.)?-1.:1.);
            case 2: return Point<2>({dGVec(1,0),-dGVec(0,0)})/Point<2>({dGVec(1,0),-dGVec(0,0)}).norm() * ((Point<2>({dGVec(1,0),-dGVec(0,0)}).dot(dGVec.col(1))>0.)?-1.:1.);
        }
    }
    return Point<2>();
};
template<> Point<3> Mesh<3>::faceNormal(const int elementID, const int face) const {
    assert(false);
    return Point<3>();
};

template<> Real Mesh<1>::faceArea(const int elementID, const int face) const{
    assert(false);
    return 1.;
}
template<> Real Mesh<2>::faceArea(const int elementID, const int face) const{
    assert(face>=0 && face<3);
    if (elementID % 2) { // elementID is odd
        switch (face) {
            case 0: return dGVec.col(0).norm();
            case 1: return (dGVec.col(0)+dGVec.col(1)).norm();
            case 2: return dGVec.col(1).norm();
        }
    } else {             // elementID is even
        switch (face) {
            case 0: return dGVec.col(1).norm();
            case 1: return (dGVec.col(0)+dGVec.col(1)).norm();
            case 2: return dGVec.col(0).norm();
        }
    }
    return 0.;
}
template<> Real Mesh<3>::faceArea(const int elementID, const int face) const{
    assert(false);
    return 0.;
}
//********************************
//* I/O
//********************************
template<int NDim> void Mesh<NDim>::writeToTxtFile(const std::string &FileName) const {
    std::ofstream myfile (FileName);
    for (int i=0; i<NDim; ++i) { myfile << origin(i,0) << " ";}
    myfile << "\n";
    for (int j=0; j<NDim; ++j) {
        for (int i=0; i<NDim; ++i) { myfile << gVec(i,j) << " ";}
            myfile << "\n";
    }
    myfile.close();
}

template<int NDim> std::ostream &operator<<(std::ostream &os, Mesh<NDim> const& mesh) {
    os <<  "********************************************\n";
    os <<  " Region: \n";
    os <<  "origin  : " + to_string(mesh.region->origin.transpose()) + "\n";
    for (int i=0; i<NDim; ++i) {
        os <<  "gVec[" + std::to_string(i) + "] : " + to_string(mesh.region->gVec.col(i).transpose()) + "\n";
    }
    os <<  "-----------------------------------\n";
    os <<  " Mesh : \n";
    os <<  "relative origin : " + to_string(mesh.relativeOrigin.transpose()) + "\n";
    os <<  "relative gVecs  : " + to_string(mesh.relativegVec.transpose()) + "\n";
    os <<  "origin  : " + to_string(mesh.origin.transpose()) + "\n";
    for (int i=0; i<NDim; ++i) {
        os <<  "gVec[" + std::to_string(i) + "] : " + to_string(mesh.gVec.col(i).transpose()) + "\n";
    }
    os <<  "Number of Sectionings per dimension : " ;
    for (int i=0; i<NDim; ++i) {
        os <<  std::to_string(mesh.nSecSplits[i]) + " ";
    }
    os <<  "\n";
    os <<  "Number of Sections                  : " + std::to_string(mesh.numberSections) + "\n";
    for (int i=0; i<NDim; ++i) {
        os <<  "dgVec[" + std::to_string(i) + "] : " + to_string(mesh.dGVec.col(i).transpose()) + "\n";
    }
    os <<  "Number of Elements                  : " + std::to_string(mesh.numberElements) + "\n";
    os <<  "Number of Basis Functions           : " + std::to_string(mesh.dataVectorDim) + "\n";
//    output+= "Volume of each element              : " + std::to_string(mesh.sqrtElemVolume*mesh.sqrtElemVolume) + "\n";
    os <<  "********************************************\n";
    return os ;
}

Plotter3D& operator << (Plotter3D& plotter, const Mesh<3>& mesh){
    plotter << OPENSTREAM << POLYGONMODE;
    for (auto elem = mesh.elementIterator(); elem.unfinished; ++elem){
        ArrayPoint<3,8> line; int i=0;
        for (auto vertex = mesh.vertexInElemIterator(); vertex.unfinished ; ++vertex){
            line.col(i++) = mesh(elem,vertex);
        }
        line.col(i++) = mesh(elem.currentElementID,0);
        line.col(i++) = mesh(elem.currentElementID,2);
        line.col(i++) = mesh(elem.currentElementID,1);
        line.col(i++) = mesh(elem.currentElementID,3);
        plotter << line;
    }
    return plotter << CLOSESTREAM ;
}
Plotter3D& operator << (Plotter3D& plotter, const Mesh<2>& mesh){
    plotter << OPENSTREAM << POLYGONMODE;
    for (auto elem = mesh.elementIterator(); elem.unfinished; ++elem){
        ArrayPoint<2,2+1> triang; int i=0;
        for (auto vertex = mesh.vertexInElemIterator(); vertex.unfinished ; ++vertex){
            triang.col(i++) = mesh(elem,vertex);
        }
        plotter << triang;
    }
    return plotter << CLOSESTREAM ;
}
Plotter2D& operator << (Plotter2D& plotter, const Mesh<2>& mesh){
    plotter << OPENSTREAM << POLYGONMODE;
    for (auto elem = mesh.elementIterator(); elem.unfinished; ++elem){
        ArrayPoint<2,2+1> triang; int i=0;
        for (auto vertex = mesh.vertexInElemIterator(); vertex.unfinished ; ++vertex){
            triang.col(i++) = mesh(elem,vertex);
        }
        plotter << triang;
    }
    return plotter << CLOSESTREAM ;
}
Plotter2D& operator << (Plotter2D& plotter, const Mesh<1>& mesh){
    VectorPoint<1> vertices;
    vertices.resize(1, mesh.numberSections + 1);
    int i = 0;
    for (auto elem = mesh.elementIterator(); elem.unfinished; ++elem){
        vertices.col(i++) = mesh(elem,mesh.vertexInElemIterator());
    }
    vertices.col(mesh.numberSections) = mesh.origin + mesh.gVec.col(0);
    return plotter << ASPOINTS(vertices) ;
}

template<> void Mesh<3>::plot(const std::string& t_title) const {
    plotter3d << *region << NOLABEL << COLOR("blue")  << *this ;
    plotter3d << NOLABEL << COLOR("red");
    if (t_title!= "") plotter3d << PLOTTITLE(t_title);
    plotter3d << PLOT;
}
template<> void Mesh<2>::plot(const std::string& t_title) const {
    plotter3d << *region << NOLABEL << COLOR("blue")  << *this ;
    plotter3d << NOLABEL << COLOR("red");
    if (t_title!= "") plotter3d << PLOTTITLE(t_title);
    plotter3d << PLOT;
}
template<> void Mesh<1>::plot(const std::string& t_title) const {
    plotter2d << *region ;
    plotter2d << *this << NOLABEL << COLOR("red");
    if (t_title!= "") plotter2d << PLOTTITLE(t_title);
    plotter2d << PLOT;
}

template<> void Mesh<3>::plotDetailed(const std::string& t_title) const {
    plotter3d << *this << NOLABEL << COLOR("red");
    if (t_title!= "") plotter3d << PLOTTITLE(t_title);
    for (int i = 0; i<this->numberElements; ++i){
        plotter3d << TEXT(std::to_string(i), this->elemCentre(i));
    }
    plotter3d << PLOT;
}
template<> void Mesh<2>::plotDetailed(const std::string& t_title) const {
    plotter3d << *this << NOLABEL << COLOR("red");
    if (t_title!= "") plotter3d << PLOTTITLE(t_title);
    for (int i = 0; i<this->numberElements; ++i){
        plotter3d << TEXT(std::to_string(i), this->elemCentre(i));
    }
    plotter3d << PLOT;
}
template<> void Mesh<1>::plotDetailed(const std::string& t_title) const {
    plotter2d <<  *this << NOLABEL << COLOR("red");
    if (t_title!= "") plotter2d << PLOTTITLE(t_title);
    for (int i = 0; i<this->numberElements; ++i){
        plotter2d << TEXT(std::to_string(i), this->elemCentre(i));
    }
    plotter2d << PLOT;
}

template<int NDim> bool operator==(const Mesh<NDim>& lhs, const Mesh<NDim>& rhs) {
    bool toreturn = true;
    toreturn = toreturn && (lhs.region == rhs.region);
    toreturn = toreturn && (lhs.relativeOrigin == rhs.relativeOrigin);
    toreturn = toreturn && (lhs.relativegVec == rhs.relativegVec);
    toreturn = toreturn && (lhs.nSecSplits == rhs.nSecSplits);
    return toreturn;
}

// Explicit template instantiations (For more info see: https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl)

template class Mesh<1>;
template class Mesh<2>;
template class Mesh<3>;

template  std::ostream &operator<<(std::ostream &os, const Mesh<1> & mesh);
template  std::ostream &operator<<(std::ostream &os, const Mesh<2> & mesh);
template  std::ostream &operator<<(std::ostream &os, const Mesh<3> & mesh);

template bool operator==(const Mesh<1>& lhs, const Mesh<1>& rhs);
template bool operator==(const Mesh<2>& lhs, const Mesh<2>& rhs);
template bool operator==(const Mesh<3>& lhs, const Mesh<3>& rhs);

} // namespace Tortoise
