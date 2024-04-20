//
//  GaussQuadratureCartConst.hpp
//  TORTOISE_Development
//
//  Created by Marco Battiato on 2/11/22.
//

#ifndef GaussQuadratureCartConst_hpp
#define GaussQuadratureCartConst_hpp

namespace Tortoise {

namespace GeometryCore {

// Gauss Points for cartesian product of simplexes
// Authors: Chua Jia Yi, Marco Battiato
// It creates the Gauss points and weights for the cartesian product of nFullSimplexes simplexes of NDim dimensions and, if NDim>1, a simplex of NDim-1 dimensions

// The order of the Gauss integration for the cartesian product of simplexes is hardcoded below.
// This is done since it is not really useful to modify it, after the initial testing.
// Moreover TORTOISE does not need the flexibility of changing it in other parts of the library.
// This is because this will only be used in the integration of the scattering integrals. In that case
// the integrands will never be simple polynomials, so the Gauss integration will never be exact.

const int integrOrderForGaussCartSimplexes = 2;   // The long name has been chosen to avoid naming conflicts

template <int NDim, int integrOrder> void recursiveForGenerateGaussPointsCartFull(int numberOfSimplexes, std::vector<int>& indexVec, std::vector<VectorPoint<NDim>>& gaussPointsCart) {
// Authors: Chua Jia Yi, Marco Battiato
    if (numberOfSimplexes != 0) {
        for (int i = 0; i < ngausspoint(NDim, integrOrder); i++) {
            indexVec[numberOfSimplexes - 1] = i;             // indices are stored i, j, n from the back of indexVec
            recursiveForGenerateGaussPointsCartFull<NDim, integrOrder>(numberOfSimplexes - 1, indexVec, gaussPointsCart);
        };
    } else {
        // Iteration over the Gauss points of the partial leg. Notice that in 1D the partial leg does not exist, and the for-loop simply runs a single iteration
        for (int j = 0; j < (NDim >= 2 ? ngausspoint(NDim - 1, integrOrder) : 1); j++) {
            int colIndex = 0;
            for (int i = 0; i < indexVec.size(); i++) {
                colIndex += indexVec[i] * std::pow(ngausspoint(NDim, integrOrder), i);
            }
            colIndex = (colIndex * (NDim >= 2 ? ngausspoint(NDim - 1, integrOrder) : 1)) + j;

            for (int i = 0; i < indexVec.size(); i++) {
                gaussPointsCart[i].col(colIndex) = gaussPoints<NDim, integrOrder>.col(indexVec[indexVec.size() - 1 - i]);
            }
        }
    }
};
template <int NDim, int integrOrder> void recursiveForGenerateGaussPointsCartPartial(int numberOfSimplexes, std::vector<int>& indexVec, VectorPoint<NDim-1>& gaussPointsCartPart) {
// Authors: Chua Jia Yi, Marco Battiato
    if (numberOfSimplexes != 0) {
        for (int i = 0; i < ngausspoint(NDim, integrOrder); i++) {
            indexVec[numberOfSimplexes - 1] = i;             // indices are stored i, j, n from the back of indexVec
            recursiveForGenerateGaussPointsCartPartial<NDim, integrOrder>(numberOfSimplexes - 1, indexVec, gaussPointsCartPart);
        };
    } else {
        // Iteration over the Gauss points of the partial leg. Notice that in 1D the partial leg does not exist, and the for loop simply runs a single iteration
        for (int j = 0; j < (NDim >= 2 ? ngausspoint(NDim - 1, integrOrder) : 1); j++) {
            int colIndex = 0;
            for (int i = 0; i < indexVec.size(); i++) {
                colIndex += indexVec[i] * std::pow(ngausspoint(NDim, integrOrder), i);
            }
            colIndex = (colIndex * (NDim >= 2 ? ngausspoint(NDim - 1, integrOrder) : 1)) + j;

            gaussPointsCartPart.col(colIndex) = gaussPoints<NDim - 1, integrOrder>.col(j);
        }
    }
};
template <int NDim, int integrOrder> void recursiveForGenerateGaussPointsCartWeights(int numberOfSimplexes, std::vector<int>& indexVec, VectorPoint<1>& gaussWeightsCart){
// Authors: Chua Jia Yi, Marco Battiato
    if (numberOfSimplexes != 0) {
        for (int i = 0; i < ngausspoint(NDim, integrOrder); i++) {
            indexVec[numberOfSimplexes - 1] = i;             // indices are stored i, j, n from the back of indexVec
            recursiveForGenerateGaussPointsCartWeights<NDim, integrOrder>(numberOfSimplexes - 1, indexVec, gaussWeightsCart);
        };
    } else {
        // Iteration over the Gauss points of the partial leg. Notice that in 1D the partial leg does not exist, and the for-loop simply runs a single iteration
        for (int j = 0; j < (NDim >= 2 ? ngausspoint(NDim - 1, integrOrder) : 1); j++) {
            int colIndex = 0;
            for (int i = 0; i < indexVec.size(); i++) {
                colIndex += indexVec[i] * std::pow(ngausspoint(NDim, integrOrder), i);
            }
            colIndex = (colIndex * (NDim >= 2 ? ngausspoint(NDim - 1, integrOrder) : 1)) + j;

            for (int i = 0; i < indexVec.size(); i++) {
                gaussWeightsCart[colIndex] *= gaussWeigths<NDim, integrOrder>(indexVec[i]);
            }
            if (NDim >= 2){ gaussWeightsCart[colIndex] *= gaussWeigths<NDim - 1, integrOrder>(j);}
        }
    }
};

template <int NDim, int fullLegsNumber> std::vector<VectorPoint<NDim>> constructGaussPointsCartFullLegs(){
    // Authors: Chua Jia Yi, Marco Battiato
    std::vector<VectorPoint<NDim>>  gaussPointsCart;           // This will be returned
    int const integrOrder = integrOrderForGaussCartSimplexes;  // Defined to be able to use a shorter name
    // Allocate a number of Eigen::Matrix<Real, NDim, Eigen::Dynamic> equal to the number of simplexes at full dimensions
    gaussPointsCart.resize(fullLegsNumber);
    for (int i = 0; i < fullLegsNumber; i++) {
        gaussPointsCart[i].resize(Eigen::NoChange, std::pow(ngausspoint(NDim, integrOrder), fullLegsNumber) * (NDim >= 2 ? ngausspoint(NDim - 1, integrOrder) : 1));
    };
    std::vector<int>  indexVec;
    indexVec.resize(fullLegsNumber);
    recursiveForGenerateGaussPointsCartFull<NDim, integrOrder>(fullLegsNumber, indexVec, gaussPointsCart);
    return gaussPointsCart;
}
template <int NDim, int fullLegsNumber> VectorPoint<NDim-1> constructGaussPointsCartPartialLeg(){
    // Authors: Chua Jia Yi, Marco Battiato
    VectorPoint<NDim-1> gaussPointsCartPart;                  // This will be returned
    int const integrOrder = integrOrderForGaussCartSimplexes; // Defined to be able to use a shorter name
    // Notice that the following will not work for NDim == 1, however it will never be executed for that number of dimensions
    gaussPointsCartPart.resize(Eigen::NoChange, std::pow(ngausspoint(NDim, integrOrder), fullLegsNumber) * ngausspoint(NDim - 1, integrOrder));
    std::vector<int>  indexVec;
    indexVec.resize(fullLegsNumber);
    recursiveForGenerateGaussPointsCartPartial<NDim, integrOrder>(fullLegsNumber, indexVec, gaussPointsCartPart);
    return gaussPointsCartPart;
}
template <int NDim, int fullLegsNumber> VectorPoint<1> constructGaussWeigthsCart(){
    // Authors: Chua Jia Yi, Marco Battiato
    VectorPoint<1> gaussWeigthsCart; // This will be returned
    int const integrOrder = integrOrderForGaussCartSimplexes; // Defined to be able to use a shorter name
    gaussWeigthsCart.resize(Eigen::NoChange, std::pow(ngausspoint(NDim, integrOrder), fullLegsNumber) * (NDim >= 2 ? ngausspoint(NDim - 1, integrOrder) : 1));
    // Start with gaussWeights = 1
    for (int i = 0; i < gaussWeigthsCart.size(); i++) { gaussWeigthsCart[i] = 1.; };
    std::vector<int>  indexVec;
    indexVec.resize(fullLegsNumber);
    recursiveForGenerateGaussPointsCartWeights<NDim, integrOrder>(fullLegsNumber, indexVec, gaussWeigthsCart);
    return gaussWeigthsCart;
}

template <int NDim, int fullLegsNumber> const std::vector<VectorPoint<NDim>>   gaussPointsCartFullLegs;
template <int NDim, int fullLegsNumber> const VectorPoint<NDim-1> gaussPointsCartPartialLeg;  // Notice that this variable will not exist for NDim == 1
template <int NDim, int fullLegsNumber> const VectorPoint<1>      gaussWeigthsCart ;

template <> const inline std::vector<VectorPoint<1>> gaussPointsCartFullLegs<1,1>      = constructGaussPointsCartFullLegs<1,1>();
template <> const inline VectorPoint<1>              gaussWeigthsCart<1,1>             = constructGaussWeigthsCart<1,1>();

template <> const inline std::vector<VectorPoint<1>> gaussPointsCartFullLegs<1,2>      = constructGaussPointsCartFullLegs<1,2>();
template <> const inline VectorPoint<1>              gaussWeigthsCart<1,2>             = constructGaussWeigthsCart<1,2>();

template <> const inline std::vector<VectorPoint<1>> gaussPointsCartFullLegs<1,3>      = constructGaussPointsCartFullLegs<1,3>();
template <> const inline VectorPoint<1>              gaussWeigthsCart<1,3>             = constructGaussWeigthsCart<1,3>();

template <> const inline VectorPoint<1>              gaussPointsCartPartialLeg<2,0>    = constructGaussPointsCartPartialLeg<2,0>();
template <> const inline VectorPoint<1>              gaussWeigthsCart<2,0>             = constructGaussWeigthsCart<2,0>();

template <> const inline std::vector<VectorPoint<2>> gaussPointsCartFullLegs<2,1>      = constructGaussPointsCartFullLegs<2,1>();
template <> const inline VectorPoint<1>              gaussPointsCartPartialLeg<2,1>    = constructGaussPointsCartPartialLeg<2,1>();
template <> const inline VectorPoint<1>              gaussWeigthsCart<2,1>             = constructGaussWeigthsCart<2,1>();

template <> const inline std::vector<VectorPoint<2>> gaussPointsCartFullLegs<2,2>      = constructGaussPointsCartFullLegs<2,2>();
template <> const inline VectorPoint<1>              gaussPointsCartPartialLeg<2,2>    = constructGaussPointsCartPartialLeg<2,2>();
template <> const inline VectorPoint<1>              gaussWeigthsCart<2,2>             = constructGaussWeigthsCart<2,2>();

template <> const inline std::vector<VectorPoint<2>> gaussPointsCartFullLegs<2,3>      = constructGaussPointsCartFullLegs<2,3>();
template <> const inline VectorPoint<1>              gaussPointsCartPartialLeg<2,3>    = constructGaussPointsCartPartialLeg<2,3>();
template <> const inline VectorPoint<1>              gaussWeigthsCart<2,3>             = constructGaussWeigthsCart<2,3>();

template <> const inline VectorPoint<2>              gaussPointsCartPartialLeg<3,0>    = constructGaussPointsCartPartialLeg<3,0>();
template <> const inline VectorPoint<1>              gaussWeigthsCart<3,0>             = constructGaussWeigthsCart<3,0>();

template <> const inline std::vector<VectorPoint<3>> gaussPointsCartFullLegs<3,1>      = constructGaussPointsCartFullLegs<3,1>();
template <> const inline VectorPoint<2>              gaussPointsCartPartialLeg<3,1>    = constructGaussPointsCartPartialLeg<3,1>();
template <> const inline VectorPoint<1>              gaussWeigthsCart<3,1>             = constructGaussWeigthsCart<3,1>();

template <> const inline std::vector<VectorPoint<3>> gaussPointsCartFullLegs<3,2>      = constructGaussPointsCartFullLegs<3,2>();
template <> const inline VectorPoint<2>              gaussPointsCartPartialLeg<3,2>    = constructGaussPointsCartPartialLeg<3,2>();
template <> const inline VectorPoint<1>              gaussWeigthsCart<3,2>             = constructGaussWeigthsCart<3,2>();

template <> const inline std::vector<VectorPoint<3>> gaussPointsCartFullLegs<3,3>      = constructGaussPointsCartFullLegs<3,3>();
template <> const inline VectorPoint<2>              gaussPointsCartPartialLeg<3,3>    = constructGaussPointsCartPartialLeg<3,3>();
template <> const inline VectorPoint<1>              gaussWeigthsCart<3,3>             = constructGaussWeigthsCart<3,3>();

// Below are routines to construct Gauss integration in domains which are split

template <int integrOrder> void generateGaussPointsForSplitElement1D(
    const int numberOfSplits,
    Eigen::Matrix<Real, 1, Eigen::Dynamic> & gaussPointsCart,
    Eigen::Matrix<Real, 1, Eigen::Dynamic> & gaussWeightsCart)
{
    //resize carts to number of splits
    gaussPointsCart.resize(Eigen::NoChange, ngausspoint(1, integrOrder) * numberOfSplits);
    gaussWeightsCart.resize(Eigen::NoChange, ngausspoint(1, integrOrder) * numberOfSplits);

    //set of gauss points of starting split element
    Eigen::Matrix<Real, 1, ngausspoint(1, integrOrder)> firstGaussPoint;
    Eigen::Matrix<Real, 1, ngausspoint(1, integrOrder)> gaussWeight;

    firstGaussPoint = (gaussPoints<1, integrOrder>.array() / numberOfSplits).matrix();
    gaussWeight = (gaussWeigths<1, integrOrder>.array() / numberOfSplits).matrix();

    for (int i = 0; i < numberOfSplits; i++) {
        gaussPointsCart.block<1, ngausspoint(1, integrOrder)>(0, i * ngausspoint(1, integrOrder)) = (firstGaussPoint.array() + (i * (1.0 / numberOfSplits))).matrix();
        gaussWeightsCart.block<1, ngausspoint(1, integrOrder)>(0, i * ngausspoint(1, integrOrder)) = gaussWeight;
    };
}


template <int integrOrder> void generateGaussPointsForSplitElement2D(
    const int numberOfSplits,
    Eigen::Matrix<Real, 2, Eigen::Dynamic> & gaussPointsCart,
    Eigen::Matrix<Real, 1, Eigen::Dynamic> & gaussWeightsCart)
{
    //resize carts to number of splits
    gaussPointsCart.resize(Eigen::NoChange, ngausspoint(2, integrOrder) * std::pow(numberOfSplits, 2));
    gaussWeightsCart.resize(Eigen::NoChange, ngausspoint(2, integrOrder) * std::pow(numberOfSplits, 2));
   
    //set of gauss points of starting split element
    Eigen::Matrix<Real, 2, ngausspoint(2, integrOrder)> firstGaussPoint;
    Eigen::Matrix<Real, 2, ngausspoint(2, integrOrder)> flipGaussPoint;
    Eigen::Matrix<Real, 1, ngausspoint(2, integrOrder)> gaussWeight;

    gaussWeight = (gaussWeigths<2, integrOrder>.array() / std::pow(numberOfSplits,2)); //gaussWeight no change

    int elementIndex = 0;
    for (int i = numberOfSplits; i > 0; i--) {
        //reset first gausspoint of column becoz y coordinate gets accumulated
        firstGaussPoint = (gaussPoints<2, integrOrder>.array() / numberOfSplits).matrix();
        flipGaussPoint = ((1.0 / numberOfSplits) - firstGaussPoint.array()).matrix();
       
        //modify x coordinate for each column
        firstGaussPoint.row(0) = (firstGaussPoint.row(0).array() + ((numberOfSplits - i) * (1.0 / numberOfSplits))).matrix();
        flipGaussPoint.row(0) = (flipGaussPoint.row(0).array() + ((numberOfSplits - i) * (1.0 / numberOfSplits))).matrix();

        for (int j = 0; j < i * 2 - 1; j++) {
            if (j % 2 == 0) {
                gaussPointsCart.block<2, ngausspoint(2, integrOrder)>(0, (elementIndex + j) * ngausspoint(2, integrOrder)) = firstGaussPoint;
                firstGaussPoint.row(1) = (firstGaussPoint.row(1).array() + (1.0 / numberOfSplits)).matrix();
            }
            else {
                gaussPointsCart.block<2, ngausspoint(2, integrOrder)>(0, (elementIndex + j) * ngausspoint(2, integrOrder)) = flipGaussPoint;
                flipGaussPoint.row(1) = (flipGaussPoint.row(1).array() + (1.0 / numberOfSplits)).matrix();
            }
            gaussWeightsCart.block<1, ngausspoint(2, integrOrder)>(0, (elementIndex + j) * ngausspoint(2, integrOrder)) = gaussWeight;
        };
        //add previous number of elements after each column
        elementIndex += i * 2 - 1;
    };
}


template <int integrOrder> void recursionForSplits(
    const int numberOfTriangles,
    const std::vector<int>& numberOfSplits,
    std::vector<int>& indexVec,
    std::vector< Eigen::Matrix<Real, 2, Eigen::Dynamic> >& gaussPointsCartFull,
    Eigen::Matrix<Real, 1, Eigen::Dynamic>& gaussPointsCartPartial,
    Eigen::Matrix<Real, 1, Eigen::Dynamic>& gaussWeightsCart)
{
    if (numberOfTriangles != 0) {
        for (int i = 0; i < ngausspoint(2, integrOrder) * std::pow(numberOfSplits[numberOfSplits.size() - numberOfTriangles - 1],2); i++) {
            indexVec[numberOfTriangles - 1] = i;             // indices are stored i, j, n from the back of indexVec
            recursionForSplits<integrOrder>(numberOfTriangles - 1, numberOfSplits, indexVec, gaussPointsCartFull, gaussPointsCartPartial, gaussWeightsCart);
        };
    }
    else {
        for (int i = 0; i < ngausspoint(1, integrOrder) * numberOfSplits[numberOfSplits.size() - 1]; i++) {
            int colIndex = 0;
            int times = 1;
            for (int j = 0; j < indexVec.size(); j++) {
                colIndex += indexVec[j] * times;
                times *= ngausspoint(2,integrOrder) * std::pow(numberOfSplits[numberOfSplits.size() - 2 - j],2);
            };
            colIndex = (colIndex * ngausspoint(1, integrOrder)) * numberOfSplits[numberOfSplits.size() - 1] + i;

            Eigen::Matrix<Real, 2, Eigen::Dynamic> gaussPointsCart2D;
            Eigen::Matrix<Real, 1, Eigen::Dynamic> gaussWeightsCartTemp;
            for (int j = 0; j < indexVec.size(); j++) {
                generateGaussPointsForSplitElement2D<integrOrder>(numberOfSplits[j], gaussPointsCart2D, gaussWeightsCartTemp);
                gaussPointsCartFull[j].col(colIndex) = gaussPointsCart2D.col(indexVec[indexVec.size() - 1 - j]);
                gaussWeightsCart[colIndex] *= gaussWeightsCartTemp(indexVec[indexVec.size() - 1 - j]);
            };
            //add last columns
            Eigen::Matrix<Real, 1, Eigen::Dynamic> gaussPointsCart1D;
            generateGaussPointsForSplitElement1D<integrOrder>(numberOfSplits[numberOfSplits.size() - 1], gaussPointsCart1D, gaussWeightsCartTemp);
            gaussPointsCartPartial.col(colIndex) = gaussPointsCart1D.col(i);
            gaussWeightsCart[colIndex] *= gaussWeightsCartTemp(i);
        }
    }
}


template <int integrOrder> void generateGaussPointsForIntegrationSplitTriangles(
    const std::vector<int>& numberOfSplits,
    std::vector< Eigen::Matrix<Real, 2, Eigen::Dynamic> >& gaussPointsCartFull,
    Eigen::Matrix<Real, 1, Eigen::Dynamic>& gaussPointsCartPartial,
    Eigen::Matrix<Real, 1, Eigen::Dynamic>& gaussWeightsCart)
{
    int numberOfTriangles = numberOfSplits.size() - 1;
    gaussPointsCartFull.clear();
    gaussPointsCartFull.resize(numberOfTriangles);

    int totalGaussPoints = 1;
    for (int i = 0; i < numberOfTriangles; i++) {
        totalGaussPoints *= ngausspoint(2, integrOrder) * std::pow(numberOfSplits[i], 2);
    }
    totalGaussPoints *= ngausspoint(1, integrOrder) * numberOfSplits[numberOfSplits.size() - 1];
   
    for (int i = 0; i < numberOfTriangles; i++) {
        gaussPointsCartFull[i].resize(Eigen::NoChange, totalGaussPoints);
    }
    gaussPointsCartPartial.resize(Eigen::NoChange, totalGaussPoints);
    gaussWeightsCart.resize(Eigen::NoChange, totalGaussPoints);

    std::vector<int> indexVec;
    indexVec.resize(numberOfTriangles);

    //start with gaussWeights = 1
    for (int i = 0; i < gaussWeightsCart.size(); i++) {
        gaussWeightsCart[i] = 1.0;
    };

    recursionForSplits<integrOrder>(numberOfTriangles, numberOfSplits, indexVec, gaussPointsCartFull, gaussPointsCartPartial, gaussWeightsCart);
}

} // namespace GeometryCore 

} // namespace Tortoise 


#endif /* GaussQuadratureCartConst_hpp */
