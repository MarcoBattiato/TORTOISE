//
//  ReferenceElement.hpp
//  TORTOISE_Development
//
//  Created by Marco Battiato on 2/11/22.
//

#ifndef ReferenceElement_hpp
#define ReferenceElement_hpp

#include <Geometry/GeometryCore/Point.hpp>

namespace Tortoise {

template <int NDim> using CartIndex = Eigen::Matrix<int, NDim, 1> ;

namespace GeometryCore {

// Specialised types with special meaning
template <int NDim> using ElemNodalData             = Eigen::Matrix<Real, 1, NDim+1>;
// Values at nodes of a tetrahedron
template <int NDim> using Element                   = ArrayPoint<NDim, NDim+1>;
// Stores the vertices of a given element


// ===  Constants
// ================
template <int NDim> Element<NDim> refNodes = [] {
    ArrayPoint<NDim, NDim+1> tmp = ArrayPoint<NDim, NDim+1>::Zero();
    for (int i=0; i<NDim; ++i){ tmp(i,i+1) = 1.0; }
    return tmp; }();

} // namespace GeometryCore 

} // namespace Tortoise

#endif /* ReferenceElement_hpp */
