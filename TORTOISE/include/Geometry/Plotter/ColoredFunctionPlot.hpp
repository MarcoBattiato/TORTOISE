//
//  ColoredFunctionPlot.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 25/1/24.
//

#ifndef ColoredFunctionPlot_hpp
#define ColoredFunctionPlot_hpp

#include <Geometry/Plotter/Plotter.hpp>
#include <Geometry/Structured/FunctionSpace/Function.hpp>

namespace Tortoise {

class COLOREDPLOT{
    const Function<2> & comparison_function_one;
    const Function<2> & comparison_function_two;
public:
    COLOREDPLOT(const Function<2>& function_one, const Function<2>& function_two);
};

Plotter3D& operator<<(Plotter3D& plotter3D, const COLOREDPLOT& colourPlot);

} // namespace Tortoise

#endif /* ColoredFunctionPlot_hpp */
