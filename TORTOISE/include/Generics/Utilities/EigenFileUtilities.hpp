//
//  EigenFileUtilities.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 4/8/23.
//

#ifndef EigenFileUtilities_hpp
#define EigenFileUtilities_hpp

#include <Eigen/Dense>
#include <string>
#include <iostream>
#include <fstream>

namespace Tortoise {

namespace Utilities {

template<typename Derived> void writeEigenToTxtFile(const Eigen::DenseBase<Derived>& eigenObject, const std::string& fileName){
    std::ofstream outFile (fileName);
    outFile << std::setprecision( std::numeric_limits<typename Eigen::DenseBase<Derived>::Scalar>::digits10+2) <<  eigenObject.reshaped().transpose();
    outFile.close();
}

template<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
void readEigenFromTxtFile(Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>& eigenObject, const std::string& fileName){
    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic,1>> map(eigenObject.data(), eigenObject.size());
    const long int numberAllocatedElements = eigenObject.rows() * eigenObject.cols();
    
    std::ifstream inFile (fileName);
    std::string line;
    long int i = 0;
    while (std::getline(inFile >> std::ws, line, ' ') ) {
        assert(i < numberAllocatedElements);
        map(i++) = std::stod(line);
    }
    inFile.close();
}

template<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
void readEigenFromTxtFile(Eigen::Array<Scalar, RowsAtCompileTime, ColsAtCompileTime>& eigenObject, const std::string& fileName){
    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic,1>> map(eigenObject.data(), eigenObject.size());
    const long int numberAllocatedElements = eigenObject.rows() * eigenObject.cols();
    
    std::ifstream inFile (fileName);
    std::string line;
    long int i = 0;
    while (std::getline(inFile >> std::ws, line, ' ') ) {
        assert(i < numberAllocatedElements);
        map(i++) = std::stod(line);
    }
    inFile.close();
}

} // namespace Utilities
} // namespace Tortoise

#endif /* EigenFileUtilities_hpp */
