//
//  main.cpp
//  Workshop
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
    brillouinZone.plot("1: Region");
    
    int                 resolution = 4;
    Mesh<2>             meshFull(brillouinZone, {0.,0.},{1.,1.},{resolution,resolution});
    meshFull.plot("2: MeshFull");

    Mesh<2>             meshPart(brillouinZone, {0.3,0.3},{.6,.6},{resolution,resolution});
    meshPart.plot("3: MeshPart");
    
    Function<2>         fun1(meshFull, x*x+sin(y+1.2*x/(y+0.1)));
    fun1.plot("4: fun1");
    
    Function<2>         fun2(meshFull, x*x+y+3);
    Function<2>         fun3(meshFull, x*x+cos(y)+3);
    (fun1 + 5 * cos(fun2)/ fun3).plot("5: Operations on functions");
    fun3 = fun1 + 5 * cos(fun2)/ fun3;
    
    Function<2>         fun5(meshPart, x*x+sin(y+1.2*x/(y+0.1)));
    
    auto lambdaExpress = [](Point<2> p){ return p(0) + 3.*p(1); };
    Function<2>         fun4(meshFull, lambdaExpress);
    
    fun4.plot("6: Function from lambda");
    
    fun3.derivative(0).plot("7: Derivative with respect to x");
    
    cout << "Max: " << fun4.max() << "   Min: " << fun4.min() << "\n";
    cout << "Integral f4: " << fun4.integrate() << "   Integral f4*f3: " << fun4.integrate(fun3)  << "   Integral f4*f3*f2: " << fun4.integrate(fun3,fun2) << "\n";
    
    
    plotter3d << meshFull;   // We send the mesh to the plotter
    for (int i=0; i <meshFull.numberElements; ++i){
        plotter3d << TEXT( i, meshFull.elemCentre(i));  // We send a text to plot
    }
    plotter3d << PLOTTITLE("Example 3.1: Element indexing") << PLOT;

    return 0;
}
