//
//  main.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 15/2/23.
//

// #define NDEBUG            // Uncomment to deactivate debug mode

#include <TORTOISE>
#include <iostream>
using namespace Tortoise;
using std::cout;
using Tortoise::AnalyticExpression::x, Tortoise::AnalyticExpression::y, Tortoise::AnalyticExpression::z;

int main(int argc, const char * argv[]) {
    
    Real                pi = 3.1415;
    Point<2>            originBZ({0.0, 0.}),
                        sideBZ0({1.0, 0.0}),
                        sideBZ1({-.5, std::sqrt(3.0)/2});
    Region<2>           brillouinZone (originBZ,{sideBZ0,sideBZ1});
    
    int                 resolutionFull = 15, resolutionPart = 7;
    Mesh<2>             meshFull(brillouinZone, {0.,0.},{1.,1.},{resolutionFull,resolutionFull});
    Mesh<2>             meshPart(brillouinZone, {0.3,0.3},{.4,.4},{resolutionPart,resolutionPart});
    
    Function<2>         fun1(meshFull, x*sin(x+3.*y));
    fun1.plot();
    
    plotter3d << brillouinZone << COLOR("black") << NOLABEL;
    plotter3d << meshFull << COLOR("green") << NOLABEL;
    plotter3d << fun1 << COLOR("green") << LABEL("My function");
    plotter3d << PLOTTITLE("New way of plotting") << PLOT;
    
    
    Region<1>           region1D(-1, 2.);
    Mesh<1>             meshFull1D(region1D, 0, 1, 50);
    Function<1>         fun1D(meshFull1D, exp(sin(5*x))-1);
    Function<1>         fun1Dother(meshFull1D, exp(sin(3*x))-1);
    
    plotter2d << region1D << COLOR("black") << NOLABEL;
    plotter2d << fun1D << COLOR("red") << LABEL("First function");
    plotter2d << fun1Dother << COLOR("blue") << LABEL("Second function");
    plotter2d << PLOTTITLE("2D plotting") << PLOT;
    
    return 0;
}
