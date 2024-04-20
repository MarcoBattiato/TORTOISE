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
#include <vector>
#include <array>
#include <tuple>
#include <cmath>
#include <NearestNeighbourInterpolator.hpp>

using std::cout, std::vector, std::string, std::array;
using namespace Tortoise;
using Tortoise::AnalyticExpression::kx, Tortoise::AnalyticExpression::ky, Tortoise::AnalyticExpression::kz, Tortoise::AnalyticExpression::k;
using Tortoise::AnalyticExpression::x, Tortoise::AnalyticExpression::y, Tortoise::AnalyticExpression::z;
using Tortoise::Utilities::stopWatch, Tortoise::Utilities::Interpolator, Tortoise::Utilities::toString;

Real const pi = 3.1415;

int main(int argc, const char * argv[]) {
    
    const int dim = 3;
    Eigen::Matrix<double, 2, Eigen::Dynamic> limits;
    Eigen::Matrix<double, 1, Eigen::Dynamic> resolution;
    limits.resize(2, dim);
    limits << 1, 1, 1, 5, 5, 5;
    resolution.resize(1, dim);
    resolution << 5,5,5;
        
    auto f = [](Eigen::Array<double, 1, dim> p){return std::sin(p(0)+p(1)*p(2));};
        
    NearestNeighbourInterpolator<dim> bla(f, limits, resolution);
        
    Eigen::Matrix<double, 1, Eigen::Dynamic> point(dim);
    point(0)=2; point(1)=2; point(2)=3;
            
    Eigen::Matrix<double, 1, Eigen::Dynamic> point2(dim);
    point2(0)=2.2; point2(1)=2.1; point2(2)=3.1;
    
    cout << f(point) << " " << bla(point) << "\n\n";
    cout << f(point) << " " << bla(point2) << "\n\n";
  
    
    return 0;
}

