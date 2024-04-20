//
//  Point.hpp
//  TORTOISE_Development
//
//  Created by Marco Battiato on 2/11/22.
//

#ifndef Point_hpp
#define Point_hpp

#include <Eigen/Dense>
#include <Generics/Utilities/UsefulMathFunctions.hpp>

namespace Tortoise {

using Real = double;

//*******************************
// Point   P
//*******************************
// ===  Types
// ================

template <int NDim> using Point                             = Eigen::Matrix<Real, NDim, 1>;
// A point can be initialised easily     Point<2> exmplPoint = {1.2,-2.5};


template <int NDim> using VectorPoint                       = Eigen::Matrix<Real, NDim, Eigen::Dynamic>;
// VectorPoint is designed to hold a large number of points usually generated randomly by the functions defined below)
// In VectorPoint the dimension related to the number of points is left to be Dynamic since there is a hard limit to the size of statically stored matrices


template <int NDim, int NPoints> using ArrayPoint           = Eigen::Matrix<Real, NDim, NPoints>;
// A vector of points with statically assigned dimension (to be used for small number of points)

namespace GeometryCore {

// ===  Creators
// ================

template <int NDim> VectorPoint<NDim> randomPoints(const int numberPoints){
    return 0.5*Eigen::Array<Real, NDim, Eigen::Dynamic>::Random(NDim,numberPoints)+0.5;
}                // Constructs random coordinates in the range [0,1]
template <int NDim> VectorPoint<NDim> randomPointsReference(const int numberPoints){
    VectorPoint<NDim> toreturn(NDim,numberPoints);
    toreturn = 0.5*Eigen::Array<Real, NDim, Eigen::Dynamic>::Random(NDim,numberPoints)+0.5;
    switch (NDim) {
        case 2: {
    // cut'n fold the square into the reference triangle
            for (int i= 0; i<numberPoints; ++i){
                if (toreturn(0,i)+toreturn(1,i)>1.0) {toreturn(0,i) = 1.0 -toreturn(0,i); toreturn(1,i) = 1.0 -toreturn(1,i);}
            }
            break;
        }
        case 3: {
    // cut'n fold the square into the reference tetrahedron   see http://vcg.isti.cnr.it/jgt/tetra.htm
            for (int i= 0; i<numberPoints; ++i){
                if (toreturn(0,i)+toreturn(1,i)>1.0) {toreturn(0,i) = 1.0 -toreturn(0,i);toreturn(1,i) = 1.0 -toreturn(1,i);}                                       // cut'n fold the cube into a prism
                if(toreturn(1,i)+toreturn(2,i)>1.0) {Real tmp = toreturn(2,i); toreturn(2,i) = 1.0 - toreturn(0,i) - toreturn(1,i); toreturn(1,i) = 1.0 - tmp;    // cut'n fold the prism into a tetrahedron
                } else if(toreturn(0,i)+toreturn(1,i)+toreturn(2,i)>1.0) {Real tmp = toreturn(2,i); toreturn(2,i) = toreturn(0,i) + toreturn(1,i) + toreturn(2,i) - 1.0; toreturn(0,i) = 1 - toreturn(1,i) - tmp;}
            }
            break;
        }
        default:
            break;
    }
    return  toreturn;
}       // Constructs random points in the reference tetrahedron

} // namespace GeometryCore 

} // namespace Tortoise

#endif /* Point_hpp */
