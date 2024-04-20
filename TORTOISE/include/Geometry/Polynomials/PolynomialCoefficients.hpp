////
////  PolynomialCoefficients.hpp
////  TORTOISE
////
////  Created by Marco Battiato on 20/9/22.
////
//
//#ifndef PolynomialCoefficients_hpp
//#define PolynomialCoefficients_hpp
//
//#include <Generics/Utilities/UsefulMathFunctions.hpp>
//#include <Geometry/GeometryCore/Point.hpp>
//
//namespace Tortoise {
//
//// Calculates the number of coefficients in a NDim variable polynomial of degree MaxDegree
//constexpr int polynomialCoeffN(const int NDim, const int MaxDegree){
//    return Utilities::binomial(NDim + MaxDegree, MaxDegree);
//}
//
//// Define as Polynomial the Eigen vector of appropriate length to store all the monomials' coefficients
//template<int NDim, int MaxDegree> requires (MaxDegree>=0 && NDim>0)
//using PolynomialCoeff  = Eigen::Matrix<Real, 1, polynomialCoeffN(NDim, MaxDegree)>;
//
//
//// We will need to extract the structure of the coefficients of all the monomials at a certain degree
//template<int NDim, int MaxDegree, typename Derived> auto monomialCoeffDegree0 (const Eigen::MatrixBase<Derived>& polyn) {
//    assert(polyn.size() == polynomialCoeffN(NDim, MaxDegree) && MaxDegree>=0);
//    return polyn(0,0);
//};
//template<int NDim, int MaxDegree, typename Derived> auto monomialCoeffDegree1 (const Eigen::MatrixBase<Derived>& polyn) {
//    assert(polyn.size() == polynomialCoeffN(NDim, MaxDegree) && MaxDegree>=1);
//    return polyn.template block<1,polynomialCoeffN(NDim,1)-polynomialCoeffN(NDim,0)>(0, polynomialCoeffN(NDim,0));
//};
//template<int NDim, int MaxDegree, typename Derived> auto monomialCoeffDegree2 (const Eigen::MatrixBase<Derived>& polyn) {
//    assert(polyn.size() == polynomialCoeffN(NDim, MaxDegree) && MaxDegree>=2);
//    const int startPositionInPolyn = polynomialCoeffN(NDim,1);
//    const int numberCoeffToExtract = polynomialCoeffN(NDim,2) -polynomialCoeffN(NDim,1);
//    return upperTriangular<NDim>(polyn.block(0,startPositionInPolyn,1,numberCoeffToExtract), 0);
//};
//template<int NDim, int MaxDegree, typename Derived> auto monomialCoeffDegree3 (const Eigen::MatrixBase<Derived>& polyn, const int which) {
//    assert(polyn.size() == polynomialCoeffN(NDim, MaxDegree) && MaxDegree>=3);
//    assert(which<NDim);
//    const int startPositionInPolyn = polynomialCoeffN(NDim,2) + which*(which*which-3*(which+NDim*(which-2)-NDim*NDim)+2 )/6;
//    const int numberCoeffToExtract = (NDim-which+1)*(NDim-which)/2;
//    return upperTriangular<NDim>(polyn.block(0,startPositionInPolyn,1,numberCoeffToExtract), which);
//};
//
//
//// Calculates the result of the polynomial calculated at either a single point (in which case it returns a number) or several points (in which case, it returns a vector with the results)
//template<int NDim, int MaxDegree, typename DerivedPol, typename DerivedPoint>
//auto apply(const Eigen::MatrixBase<DerivedPol>& polyn, const Eigen::MatrixBase<DerivedPoint>& point){
//    assert(polyn.size() == polynomialCoeffN(NDim, MaxDegree));
//    assert(point.rows() == NDim);
//    const auto nPoints = point.cols();
//    Eigen::Matrix<Real, 1, Eigen::Dynamic> toreturn(Eigen::Matrix<Real, 1, Eigen::Dynamic>::Zero(1, nPoints));
//    if constexpr(MaxDegree >=4) {
//        assert(false && "Polynomial Degree > 3 are not implemented");
//    }
//    if constexpr(MaxDegree >=3) {
//        for (int i = 0; i<NDim; ++i) { toreturn += (point.row(i).transpose().array() * ((point.transpose() * monomialCoeffDegree3<NDim,MaxDegree>(polyn,i)* point).diagonal().array())).matrix(); }
//    }
//    if constexpr(MaxDegree >=2) {
//        toreturn += (point.transpose() * monomialCoeffDegree2<NDim,MaxDegree>(polyn) * point).diagonal();
//    }
//    if constexpr(MaxDegree >=1) {
//        toreturn += monomialCoeffDegree1<NDim,MaxDegree>(polyn) * point;
//    }
//    if constexpr(MaxDegree >=0) {
//        toreturn += Eigen::Matrix<Real, 1, Eigen::Dynamic>::Constant(1, nPoints, monomialCoeffDegree0<NDim,MaxDegree>(polyn));
//    }
//    return toreturn;
//}
//
//
//// Same function as above, but optimised for evaluation at a single point
//template<int NDim, int MaxDegree, typename DerivedPol, typename DerivedPoint>
//Real applySingle(const Eigen::MatrixBase<DerivedPol>& polyn, const Eigen::MatrixBase<DerivedPoint>& point)
//requires (Eigen::MatrixBase<DerivedPoint>::ColsAtCompileTime == 1) {
//    assert(polyn.size() == polynomialCoeffN(NDim, MaxDegree));
//    assert(point.rows() == NDim);
//    Real toreturn = 0.;
//    if constexpr(MaxDegree >=4) {
//        assert(false && "Polynomial Degree > 3 are not implemented");
//    }
//    if constexpr(MaxDegree >=3) {
//        for (int i = 0; i<NDim; ++i) { toreturn += point(i) * (point.dot(monomialCoeffDegree3<NDim,MaxDegree>(polyn,i)* point)); }
//    }
//    if constexpr(MaxDegree >=2) {
//        toreturn += point.dot(monomialCoeffDegree2<NDim,MaxDegree>(polyn) * point);
//    }
//    if constexpr(MaxDegree >=1) {
//        toreturn += monomialCoeffDegree1<NDim,MaxDegree>(polyn).dot(point);
//    }
//    if constexpr(MaxDegree >=0) {
//        toreturn += monomialCoeffDegree0<NDim,MaxDegree>(polyn);
//    }
//    return toreturn;
//}
//
//
//} // namespace Tortoise
//
//#endif /* PolynomialCoefficients_hpp */
