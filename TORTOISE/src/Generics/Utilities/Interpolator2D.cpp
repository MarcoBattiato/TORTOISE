//
//  Interpolator2D.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 16/11/22.
//

#include "Interpolator2D.hpp"

namespace Tortoise {

namespace Utilities {

namespace Iterpolation {

std::ostream& operator<<(std::ostream& os, const RegularGridMultipleData& data){
    os << " xs: " << data.xs.transpose() << "\n";
    os << " ys: " << data.ys.transpose() << "\n";
    for (int i=0; i< data.data.size(); ++i){
        os << "--- data n: " << i << " ---\n";
        os << data.data[i] << "\n";
    }
    os << "-----------------\n";
    return  os;
}

} // namespace Iterpolation
} // namespace Utilities
} // Tortoise
