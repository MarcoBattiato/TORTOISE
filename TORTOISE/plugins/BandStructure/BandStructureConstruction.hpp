//
//  BandStructureConstruction.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 26/5/23.
//

#ifndef BandStructureConstruction_hpp
#define BandStructureConstruction_hpp

#include <stdio.h>

template <typename IDType, typename EigenDerived0, typename EigenDerived1>
std::vector<IDType> assignIndicesBySimilarity(const Eigen::MatrixBase<EigenDerived0>& vectorsToClassify,
                                              const Eigen::MatrixBase<EigenDerived1>& vectorsToCompareTo,
                                              const std::vector<IDType>& IDsToCompareTo){
    
    assert(vectorsToClassify.rows() == vectorsToCompareTo.rows());
    assert(vectorsToClassify.cols() == vectorsToCompareTo.cols());
    assert(vectorsToClassify.cols() == IDsToCompareTo.size());
    
    std::vector<IDType> IDClassified;
    
    auto overlapMatrix = (vectorsToClassify.transpose() * vectorsToCompareTo).array().abs();
    
    for (int i = 0; i < vectorsToClassify.cols(); ++i){
        Eigen::Index maxPos;
        overlapMatrix.row(i).maxCoeff(&maxPos);
        IDClassified.emplace_back(IDsToCompareTo[maxPos]);
    }
    
    return IDClassified;
}

#endif /* BandStructureConstruction_hpp */


//Eigen::Matrix<Real, 5, 5> mat = Eigen::Matrix<Real, 5, 5>::Random();
//mat(0,0) += 1.2;
//cout << mat << "\n\n =============== \n\n";
//
//Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real, 5, 5>> eigensolver(mat);
//
//cout << eigensolver.eigenvalues().transpose() << "\n\n";
//cout << eigensolver.eigenvectors() << "\n\n =============== \n\n";
//
//Eigen::Matrix<Real, 5, 5> eigenV = eigensolver.eigenvectors();
//
//Eigen::Matrix<Real, 5, 5> mat2(mat);
//mat2(0,0) += .2;
//
//Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real, 5, 5>> eigensolver2(mat2);
//
//cout << eigensolver2.eigenvalues().transpose() << "\n\n";
//cout << eigensolver2.eigenvectors() << "\n\n =============== \n\n";
//
//Eigen::Matrix<Real, 5, 5> eigenV2 = eigensolver2.eigenvectors();
//
//Eigen::Matrix<Real, 5, 1> temp2 = eigenV2.col(2);
//Eigen::Matrix<Real, 5, 1> temp3 = eigenV2.col(3);
//eigenV2.col(2) = eigenV2.col(0);
//eigenV2.col(3) = temp2;
//eigenV2.col(0) = temp3;
//
//cout << eigenV2.transpose() * eigenV << "\n\n";
//
//auto newIds = assignIndicesBySimilarity(eigenV2,eigenV,std::vector<int>{0,1,2,3,4});
//
//for (auto id : newIds)
//     cout << id << " ";
//
//cout << "\n\n";
