//
//  NearestNeighbourInterpolator.hpp
//  Tortoise
//
//  Created by Marco Battiato on 3/1/24.
//

#ifndef NearestNeighbourInterpolator_hpp
#define NearestNeighbourInterpolator_hpp

#include <Eigen/Dense>
#include <cassert>

template <int Dim> class NearestNeighbourInterpolator {
public:
    Eigen::Matrix<double, Eigen::Dynamic, 1>        data;
    Eigen::Matrix<double, 2, Eigen::Dynamic>        limits;
    Eigen::Matrix<Eigen::Index, 1, Eigen::Dynamic>  resolution;
    
    Eigen::Index                                    totPoints;
    Eigen::Matrix<Eigen::Index, 1, Eigen::Dynamic>  cumulativeResolution;
    
    template <typename FunctType, typename Derived0, typename Derived1>
        NearestNeighbourInterpolator(FunctType f, const Eigen::DenseBase<Derived0>& limits, const Eigen::DenseBase<Derived1>& resolution);
    
    template <typename FunctType, typename Derived0, typename Derived1, typename Derived2>
        NearestNeighbourInterpolator(FunctType f, const Eigen::DenseBase<Derived0>& lows, const Eigen::DenseBase<Derived1>& highs, const Eigen::DenseBase<Derived2>& resolution);

    template <typename Derived> auto operator()(const Eigen::DenseBase<Derived>& x) const;
    
    template <typename Derived> auto fromCartIndexToId(const Eigen::DenseBase<Derived>& cartIndex) const;
    auto fromIdToCartIndex(Eigen::Index id) const;
    
};

template <int Dim> template <typename FunctType, typename Derived0, typename Derived1>
NearestNeighbourInterpolator<Dim>::NearestNeighbourInterpolator(FunctType f, const Eigen::DenseBase<Derived0>& _limits, const Eigen::DenseBase<Derived1>& _resolution):
NearestNeighbourInterpolator(f, _limits.row(0), _limits.row(1), _resolution){}
template <int Dim> template <typename FunctType, typename Derived0, typename Derived1, typename Derived2>
NearestNeighbourInterpolator<Dim>::NearestNeighbourInterpolator(FunctType f, const Eigen::DenseBase<Derived0>& lows, const Eigen::DenseBase<Derived1>& highs, const Eigen::DenseBase<Derived2>& _resolution){
    assert( (lows.rows() == 1 && lows.cols() == Dim) || (lows.cols() == 1 && lows.rows() == Dim));
    assert( (highs.rows() == 1 && highs.cols() == Dim) || (highs.cols() == 1 && highs.rows() == Dim));
    assert( (_resolution.rows() == 1 && _resolution.cols() == Dim) || (_resolution.cols() == 1 && _resolution.rows() == Dim));
    
    limits.resize(2,Dim);
    
    if ( lows.rows() == 1) {
        limits.row(0) = lows;
    } else {
        limits.row(0) = lows.transpose();
    }
    if ( highs.rows() == 1) {
        limits.row(1) = highs;
    } else {
        limits.row(1) = highs.transpose();
    }
    
    if ( _resolution.rows() == 1) {
        resolution = _resolution.template cast<Eigen::Index>();
    } else {
        resolution = _resolution.transpose().template cast<Eigen::Index>();
    }
    
    cumulativeResolution.resize(Eigen::NoChange, Dim);
    cumulativeResolution(0) = 1;
    for (Eigen::Index i = 1; i < Dim; ++i){
        cumulativeResolution(i) = cumulativeResolution(i-1)*resolution(i-1);
    }
    totPoints = cumulativeResolution(Dim-1) * resolution(Dim-1);
    
    data.resize(totPoints);
    
    for (Eigen::Index i = 0; i < totPoints; ++i){
        data(i) = f( limits.row(0).array() + fromIdToCartIndex(i).template cast<double>().derived().array() * (limits.row(1) - limits.row(0)).array() / (resolution.array().template cast<double>()-1) );
    }
    
}

template <int Dim> template <typename Derived> auto NearestNeighbourInterpolator<Dim>::operator()(const Eigen::DenseBase<Derived>& x) const{
    auto   cart = (resolution.array().template cast<double>()-1) * (x.derived().array()-limits.row(0).array()) / (limits.row(1) - limits.row(0)).array();
    auto   cartRounded = cart.round().template cast<Eigen::Index>();
    auto   cartLimBott = cartRounded.cwiseMax(Eigen::Array<Eigen::Index, 1, Eigen::Dynamic>::Zero(1,Dim));
    auto   cartLimCompl = cartLimBott.cwiseMin(resolution.array());
    return data(fromCartIndexToId(cartLimCompl));
}

template <int Dim> template <typename Derived> auto NearestNeighbourInterpolator<Dim>::fromCartIndexToId(const Eigen::DenseBase<Derived>& cartIndex) const{
    assert( (cartIndex.rows() == 1 && cartIndex.cols() == Dim));
    return cartIndex.derived().matrix().template cast<Eigen::Index>().dot(cumulativeResolution.derived().matrix());
}
template <int Dim> auto NearestNeighbourInterpolator<Dim>::fromIdToCartIndex(Eigen::Index id) const{
    Eigen::Matrix<Eigen::Index, 1, Eigen::Dynamic> toReturn(1,Dim);
    for (int i =0 ; i < Dim; ++i){
        toReturn(Dim-i-1) = id/cumulativeResolution(Dim-1-i);
        id %= cumulativeResolution(Dim-1-i);
    }
    return toReturn;
}


#endif /* NearestNeighbourInterpolator_hpp */

