//
//  Interpolator2D.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 16/11/22.
//

#ifndef Interpolator2D_hpp
#define Interpolator2D_hpp

#include <Eigen/Dense>
#include <vector>
#include <cassert>
#include <iostream>

namespace Tortoise {

namespace Utilities {

namespace Iterpolation {

class RegularGridMultipleData {
    
public:
    const Eigen::Matrix<double, Eigen::Dynamic,1>  xs;
    const Eigen::Matrix<double, Eigen::Dynamic,1>  ys;
    std::vector< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > data;
    
    template<typename Derived0, typename Derived1>
    RegularGridMultipleData(const Eigen::MatrixBase<Derived0>& xs, const Eigen::MatrixBase<Derived1>& ys, int nDataSets): xs(xs), ys(ys) {
        for (int i=0; i< nDataSets; ++i){
            data.emplace_back(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(xs.size(), ys.size()));
        }
    }
    
    long size() const {
        return xs.size() * ys.size();
    }
    
    std::pair<long, long> index(long nPoint) const {
        assert (nPoint>=0 && nPoint < size());
        return std::pair<long,long>{nPoint%xs.size(), nPoint/xs.size() };
    }
    
    Eigen::Matrix<double, 2,1> operator[](long nPoint) const {
        assert (nPoint>=0 && nPoint < size());
        auto [xn, yn] = index(nPoint);
        return Eigen::Matrix<double, 2,1>{xs(xn), ys(yn) };
    }
    
    auto& operator()(long nPoint, long nDataSet){
        assert (nPoint>=0 && nPoint < size());
        assert (nDataSet >= 0 && nDataSet < data.size());
        auto [xn, yn] = index(nPoint);
        return data[nDataSet](xn, yn);
    }
    
    auto operator()(long nPoint, long nDataSet) const {
        assert (nPoint>=0 && nPoint < size());
        assert (nDataSet >= 0 && nDataSet < data.size());
        auto [xn, yn] = index(nPoint);
        return data[nDataSet](xn, yn);
    }
    
    void addData (long nPoint, long nDataSet, double datap){
        (*this)(nPoint,nDataSet) = datap;
    }
    
    template<typename Derived0> void addData (long nPoint, const Eigen::MatrixBase<Derived0>& datas){
        assert(datas.size() == data.size());
        for (int i=0; i < data.size(); ++i){
            (*this)(nPoint,i) = datas(i);
        }
    }
    
    friend std::ostream& operator<<(std::ostream& os, const RegularGridMultipleData& data);
    
};

std::ostream& operator<<(std::ostream& os, const RegularGridMultipleData& data);

template<typename Derived0, typename Derived1>
double interpolate1D(const Eigen::MatrixBase<Derived0>& xs, const Eigen::MatrixBase<Derived1>& ys, double x, bool extrapolate = false){
    // Assumes that the values are sorded with respect to x
    assert(xs.size() == ys.size());
    assert(xs.size() >= 2);
    long size = xs.size();
    long i = 0;
    if ( x >= xs(size - 2) ){
        i = size - 2;
    } else {
        while ( x > xs(i+1))
            i++;
    }
    double xL = xs(i), yL = ys(i), xR = xs(i+1), yR = ys(i+1);
    if ( !extrapolate )
    {
        if ( x < xL ) yR = yL;
        if ( x > xR ) yL = yR;
    }
    double dydx = ( yR - yL ) / ( xR - xL );
    return yL + dydx * ( x - xL );
}

class Interpolator2D : public RegularGridMultipleData {
public:
    
    template<typename Derived0, typename Derived1>
    Interpolator2D(const Eigen::MatrixBase<Derived0>& xs, const Eigen::MatrixBase<Derived1>& ys, int nDataSets): RegularGridMultipleData(xs, ys, nDataSets){}
    
    auto dataInterpolator(long nDataSet, bool extrapolate = false) const {
        return [this, nDataSet, extrapolate](Eigen::Matrix<double,2,1> point){
            Eigen::Matrix<double, Eigen::Dynamic,1> partialInterp;
            partialInterp.resize(this->ys.size());
            for (int i=0; i< this->ys.size(); ++i){
                partialInterp(i) = interpolate1D(this->xs, this->data[nDataSet].col(i), point(0), extrapolate);
            }
            return interpolate1D(this->ys, partialInterp, point(1), extrapolate);
        };
    }
    
};

} // namespace Iterpolation
} // namespace Utilities
} // Tortoise

#endif /* Interpolator2D_hpp */
