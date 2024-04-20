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
//  DataFromFile.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 15/10/20.
//
// Usage:
// > DataFromFile funct("datafile.txt");
// > std::cout << funct(3.4) << "\n";    // Notice that it can be used as a function
//
// "datafile.txt"
//  0.0  3.21
//  0.1  4.34
//  0.5  4.12
//   ...


#ifndef DataFromFile_hpp
#define DataFromFile_hpp

#include <Geometry/Structured/FunctionSpace/Function.hpp>

#include <Eigen/Dense>
#include <stdio.h>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>      // std::setw

namespace Tortoise {

namespace Utilities {

double interpolate1D(const std::vector<std::array<double,2>>& data, double x, bool extrapolate = false){
    // Assumes that the values are sorded with respect to x
    unsigned long size = data.size();
    unsigned long i = 0;
    if ( x >= data[size - 2][0] ){ i = size - 2;}
    else {while ( x > data[i+1][0] ) i++;}
    double xL = data[i][0], yL = data[i][1], xR = data[i+1][0], yR = data[i+1][1];
    if ( !extrapolate )
    {
        if ( x < xL ) yR = yL;
        if ( x > xR ) yL = yR;
    }
    double dydx = ( yR - yL ) / ( xR - xL );
    return yL + dydx * ( x - xL );
}
std::vector<double> interpolate1D(const std::vector<std::array<double,2>>& data, const std::vector<double>& x, bool extrapolate = false){
    // Assumes that the values are sorded with respect to x
    std::vector<double> toreturn(x.size());
    for (int i=0; i<x.size(); ++i){
        toreturn[i] = interpolate1D(data, x[i], extrapolate);
    }
    return toreturn;
}

class Interpolator1D {
public:
    bool extrapolate = true;
private:
    std::vector<std::array<double,2>> data;
public:
    Interpolator1D(){};
    Interpolator1D(const std::vector<std::array<double,2>>& _data): data(_data){
        prepareForInterpolation();
    };
    Interpolator1D(const std::vector<double>& _x, const std::vector<double>& _y){
        assert(_x.size() == _y.size());
        data.resize(_x.size());
        for (int i=0; i<_x.size(); ++i){
            data[i][0] = _x[i]; data[i][1] = _y[i];
        }
        prepareForInterpolation();
    };
    template <int NDim> Interpolator1D(const Function<NDim>& _x, const Function<NDim>& _y){
        assert(_x.mesh == _y.mesh);
        auto xs = _x.toNodalFromLinForm();
        std::vector<double> xdata(std::vector<double>(xs.data(), xs.data() + xs.rows() * xs.cols()));
        auto ys = _y.toNodalFromLinForm();
        std::vector<double> ydata(std::vector<double>(ys.data(), ys.data() + ys.rows() * ys.cols()));
        data.resize(xdata.size());
        for (int i=0; i<xdata.size(); ++i){
            data[i][0] = xdata[i]; data[i][1] = ydata[i];
        }
        prepareForInterpolation();
    };
    void addData(const double _x, const double _y) {
        data.emplace_back(std::array<double,2>({_x,_y}));
    }
    void prepareForInterpolation(){
            std::sort(data.begin(), data.end(), [](const std::array<double,2>& a, const std::array<double,2>& b) {
                return a[0] < b[0];
            });
        }
    double operator()(double x) const {
        return interpolate1D(data,x,extrapolate);
    }
    double operator()(Eigen::Matrix<double, 1, 1> p) const {
        return interpolate1D(data,p(0),extrapolate);
    }
    std::vector<double> operator()(const std::vector<double>& p) const {
        return interpolate1D(data,p,extrapolate);
    }
    friend std::ostream &operator<<(std::ostream &os, Interpolator1D const& interp){
        for (int i=0; i<interp.data.size(); ++i){
            std::cout << std::setw(10) << interp.data[i][0]<< " ";
        }
        std::cout << "\n";
        for (int i=0; i<interp.data.size(); ++i){
            std::cout << std::setw(10) << interp.data[i][1]<< " ";
        }
        std::cout << "\n";
        return os;
    };
};





double interpolate(const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate ){
    // Assumes that the values are sorded with respect to x
    unsigned long size = xData.size();
    
    unsigned long i = 0;
    if ( x >= xData[size - 2] ){ i = size - 2;}
    else {while ( x > xData[i+1] ) i++;}
    double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];
    if ( !extrapolate )
    {
        if ( x < xL ) yR = yL;
        if ( x > xR ) yL = yR;
    }
    double dydx = ( yR - yL ) / ( xR - xL );
    return yL + dydx * ( x - xL );
}

double interpolate(const std::vector<std::pair<double,double>> &data, double x, bool extrapolate ){
    // Assumes that the values are sorded with respect to x
    unsigned long size = data.size();
    
    unsigned long i = 0;
    if ( x >= data[size - 2].first ){ i = size - 2;}
    else {while ( x > data[i+1].first ) i++;}
    double xL = data[i].first, yL = data[i].second, xR = data[i+1].first, yR = data[i+1].second;
    if ( !extrapolate )
    {
        if ( x < xL ) yR = yL;
        if ( x > xR ) yL = yR;
    }
    double dydx = ( yR - yL ) / ( xR - xL );
    return yL + dydx * ( x - xL );
}


struct xydata { std::vector<double> x; std::vector<double> y; };

class Interpolator {
public:
    bool extrapolate = true;
public:
    std::vector<std::pair<double,double>> data;
public:
    Interpolator(){};
    Interpolator(const std::vector<double>& _x, const std::vector<double>& _y){
        assert(_x.size() == _y.size());
        for (int i = 0; i < _x.size(); ++i){
            data.push_back(std::make_pair(_x[i], _y[i]));
        }
    };
    void addData(const double _x, const double _y) { data.push_back(std::make_pair(_x, _y));}
    void prepareForInterpolation(){
        auto comparison = [](const std::pair<double,double> a, const std::pair<double,double> b){return a.first < b.first; };
        std::sort(data.begin(), data.end(), comparison);
    }
    double operator()(double p) const {return interpolate(data,p,extrapolate);}
    double operator()(Eigen::Matrix<double, 1, 1> p) const {return interpolate(data,p(0),extrapolate);}
};

xydata loadTxt(const std::string filename){
    std::string linedata;
    std::ifstream inputfile;
    inputfile.open(filename);
    xydata data;
    while (std::getline(inputfile, linedata)){
        std::stringstream linestream(linedata);
        float v1, v2;
        if (linestream >> v1 >> v2){
            data.x.push_back(v1);
            data.y.push_back(v2);
        }
    }
    inputfile.close();
    return data;
}

class DataFromFile {
public:
    bool extrapolate = true;
private:
    xydata data;
public:
    DataFromFile(const std::string &filename): data(loadTxt(filename)){};
    double operator()(double x){return interpolate(data.x,data.y,x,extrapolate);}
    double operator()(Eigen::Matrix<double, 1, 1> p ){return interpolate(data.x,data.y,p(0),extrapolate);}
};

} // namespace Utilities

} // namespace Tortoise
#endif /* DataFromFile_hpp */
