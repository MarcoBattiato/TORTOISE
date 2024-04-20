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
//  Vlasov.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 31/7/21.
//

#ifndef Vlasov_hpp
#define Vlasov_hpp


#include <Geometry/Structured/FunctionSpace/Function.hpp>

using namespace Tortoise::GeometryCore;

namespace Tortoise {

//template <int NDim, int order> const ArrayPoint<NDim,ngausspoint(NDim,order)> gaussPoints;

template <int NDim, int order> const std::vector< Eigen::Matrix<Real, NDim, GeometryCore::ngausspoint(NDim-1,order)>, Eigen::aligned_allocator<Eigen::Matrix<Real, NDim, GeometryCore::ngausspoint(NDim-1,order)>> > gaussPointsSide;

// 2D
// ==============
template<int order> const inline std::vector< Eigen::Matrix<Real, 2, GeometryCore::ngausspoint(1,order)>, Eigen::aligned_allocator<Eigen::Matrix<Real, 2, GeometryCore::ngausspoint(1,order)>> > gaussPointsSide<2,order> = []{
    std::vector< Eigen::Matrix<Real, 2, GeometryCore::ngausspoint(1,order)>, Eigen::aligned_allocator<Eigen::Matrix<Real, 2, GeometryCore::ngausspoint(1,order)>> > tmp;
    tmp.emplace_back(Eigen::Matrix<Real, 2, GeometryCore::ngausspoint(1,order)>::Zero());
    tmp.emplace_back(Eigen::Matrix<Real, 2, GeometryCore::ngausspoint(1,order)>::Zero());
    tmp.emplace_back(Eigen::Matrix<Real, 2, GeometryCore::ngausspoint(1,order)>::Zero());
    tmp[0].row(0) = (1.-GeometryCore::gaussPoints<1,order>.array()).matrix();
    tmp[0].row(1) = GeometryCore::gaussPoints<1,order>;
    tmp[1].row(1) = GeometryCore::gaussPoints<1,order>;
    tmp[2].row(0) = GeometryCore::gaussPoints<1,order>;
    return tmp;}();




template<int NDim> class Vlasov {
public:
    static Function<NDim> vlasovFlux(const Function<NDim>& density, const Point<NDim>& velocity){
        return (volumeFlux(density,velocity)-surfaceFlux(density,velocity)).applyInverseMass();
    }
    
    static Function<NDim> vlasovCentralFlux(const Function<NDim>& density, const Point<NDim>& velocity){
        return (volumeFlux(density,velocity)-surfaceCentralFlux(density,velocity)).applyInverseMass();
    }
    
//private:
    static Function<NDim> volumeFlux(const Function<NDim>& density, const Point<NDim>& velocity){
        // The method below relies on the specific choise of the basis functions
        // The derivative of the basis functions have always a representation of the type (a 0 0) (the derivative is always a constant)
        // We can therefore only include the first coefficient
        Function<NDim> toreturn(density.mesh);
        for (int dim=0; dim < NDim; ++dim){
            Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> 
                elementviewderiv(density.mesh->basisFuncDerivLF0thOrd[dim].data(),  density.mesh->numberElements, density.mesh->nVertElement);
            Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> 
                intermediate(elementviewderiv.array().colwise() * (apply(density.elementview(),GeometryCore::gaussPoints<NDim,2>)* GeometryCore::gaussWeigths<NDim,2>.transpose()).array());
            DataVector partialVolIntegral(Eigen::Map<DataVector>(intermediate.data(), intermediate.cols()*intermediate.rows()));
            toreturn.vec += density.mesh->elemVolume * velocity(dim) * partialVolIntegral;
        }
        return toreturn;
    }
    
    static Function<NDim> surfaceFlux(const Function<NDim>& density, const Point<NDim>& velocity){
        Function<NDim> toreturn(density.mesh);
        for (int elem=0; elem < density.mesh->numberElements; ++elem){
            for (int basis=0; basis < density.mesh->nVertElement; ++basis){
                for (int side=0; side < density.mesh->nVertElement; ++side){
                    GeometryCore::LinearForm<NDim> basisFunc (GeometryCore::LinearForm<NDim>::Zero() );
                    basisFunc(basis)=1.;
                    Real vel_SurfNorm_scalarProd = velocity.dot(density.mesh->faceNormal(elem,side));
                    int elemToInteg, sideToIntegr;
                    if ( vel_SurfNorm_scalarProd >0. ) {  // Outgoing flux
                        elemToInteg = elem;
                        sideToIntegr = side;
                    } else {                              // Ingoing flux
                        elemToInteg = density.mesh->adjacentElement(elem, side);
                        sideToIntegr = density.mesh->adjacentFace(elem, side);
                    }
                    toreturn.vec(0,elem*(NDim+1)+basis) += vel_SurfNorm_scalarProd *
                    ( apply(basisFunc, gaussPointsSide<NDim,2>[side]).array() *
                     apply(density.vec.template block<1,Mesh<NDim>::nVertElement>(0,elemToInteg*(NDim+1)), gaussPointsSide<NDim,2>[sideToIntegr]).array() *
                     GeometryCore::gaussWeigths<NDim-1,2>.array() ).sum()  * density.mesh->faceArea(elem, side);
                }
            }
        }
        return toreturn;
    }
    
    static Function<NDim> surfaceCentralFlux(const Function<NDim>& density, const Point<NDim>& velocity){
        Function<NDim> toreturn(density.mesh);
        for (int elem=0; elem < density.mesh->numberElements; ++elem){
            for (int basis=0; basis < density.mesh->nVertElement; ++basis){
                for (int side=0; side < density.mesh->nVertElement; ++side){
                    GeometryCore::LinearForm<NDim> basisFunc (GeometryCore::LinearForm<NDim>::Zero() );
                    basisFunc(basis)=1.;
                    Real vel_SurfNorm_scalarProd = velocity.dot(density.mesh->faceNormal(elem,side));
                    // Outgoing flux
                    int elemToInteg = elem, sideToIntegr = side;
                    toreturn.vec(0,elem*(NDim+1)+basis) += 0.5 * vel_SurfNorm_scalarProd *
                    ( apply(basisFunc, gaussPointsSide<NDim,2>[side]).array() *
                    apply(density.vec.template block<1,Mesh<NDim>::nVertElement>(0,elemToInteg*(NDim+1)), gaussPointsSide<NDim,2>[sideToIntegr]).array() *
                     GeometryCore::gaussWeigths<NDim-1,2>.array() ).sum()  * density.mesh->faceArea(elem, side);
                    // Ingoing flux
                    elemToInteg = density.mesh->adjacentElement(elem, side);
                    sideToIntegr = density.mesh->adjacentFace(elem, side);
                    toreturn.vec(0,elem*(NDim+1)+basis) += vel_SurfNorm_scalarProd *
                    ( apply(basisFunc, gaussPointsSide<NDim,2>[side]).array() *
                     apply(density.vec.template block<1,Mesh<NDim>::nVertElement>(0,elemToInteg*(NDim+1)), gaussPointsSide<NDim,2>[sideToIntegr]).array() *
                     GeometryCore::gaussWeigths<NDim-1,2>.array() ).sum()  * density.mesh->faceArea(elem, side);
                }
            }
        }
        return toreturn;
    }
    
    
};

} // namespace Tortoise
#endif /* Vlasov_hpp */



// NOTE: the flux due to a velocity that is a function has not been included, since for the flux to make sense, the velocity must be continuous.
// That condition cannot be enforced without changing the library
// This prevents from modeling some exotic cases, where the electric field is k dependent, but it includes all common cases

//static Function<NDim> volumeFlux(const Function<NDim>& density, const std::vector<Function<NDim>>& velocity){
//    // The method below relies on the specific choise of the basis functions
//    // The derivative of the basis functions have always a representation of the type (a 0 0) (the derivative is always a constant)
//    // We can therefore only include the first coefficient
//    Function<NDim> toreturn(density.mesh);
//    for (int dim=0; dim < NDim; ++dim){
//        Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> elementviewderiv(density.mesh->basisFuncDerivLF0thOrd[dim].data(),  density.mesh->numberElements, density.mesh->nVertElement);
//        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> intermediate
//        (elementviewderiv.array().colwise() *
//        ( (apply(density.elementview(),gaussPoints<NDim,3>).array() * apply(velocity[dim].elementview(),gaussPoints<NDim,3>).array() ).matrix() *
//        gaussWeigths<NDim,3>.transpose()).array() );
//        DataVector partialVolIntegral(Map<DataVector>(intermediate.data(), intermediate.cols()*intermediate.rows()));
//        toreturn.vec += density.mesh->elemVolume * partialVolIntegral;
//    }
//    return toreturn;
//}
