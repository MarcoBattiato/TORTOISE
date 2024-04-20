////
////  Polynomial.hpp
////  TORTOISE
////
////  Created by Marco Battiato on 20/9/22.
////
//
//#ifndef Polynomial_hpp
//#define Polynomial_hpp
//
//
//#include "Point.hpp"
//
//#include <Generics/Features/MathFieldSpace.hpp>
//#include <Generics/Features/VectorSpace.hpp>
//#include <Generics/Utilities/StringUtilities.hpp>
//#include <Generics/Utilities/UsefulMathFunctions.hpp>
//
//#include <Geometry/Polynomials/PolynomialCoefficients.hpp>
//
//#include <iostream>
//#include <vector>
//#include <string>
//#include <functional>
//
//namespace Tortoise {
//
//namespace GeometryCore {
//
//// Structures containing all the data of a polynomial
//template<int NDim, int MaxDegree> requires (MaxDegree>=0 && NDim>0) struct Polynomial;
//
//
////*******************************
//// Technicalities (Skip to line 86)
////*******************************
//
//// Calculates the result of the polynomial calculated at either a single point (in which case it returns a number) or several points (in which case, it returns a vector with the results)
//
//namespace PolynomialImplementation {
//
//extern inline const std::vector<std::string> varNames = {"x","y","z","u","v","w"};
//
//template<int NDim, int MaxDegree, typename DerivedPoint>
//auto apply(const Polynomial<NDim,MaxDegree>& polyn, const Eigen::MatrixBase<DerivedPoint>& point){
//    assert(point.rows() == NDim);
//    const auto nPoints = point.cols();
//    Eigen::Matrix<Real, 1, Eigen::Dynamic> toreturn(Eigen::Matrix<Real, 1, Eigen::Dynamic>::Zero(1, nPoints));
//    if constexpr(MaxDegree >=4) {
//        assert(false && "Polynomial Degree > 3 are not implemented");
//    }
//    if constexpr(MaxDegree >=3) {
//        for (int i = 0; i<NDim; ++i) { toreturn += (point.row(i).transpose().array()* (point.transpose() * polyn.coeff3[i] * point).diagonal().array()).matrix(); }
//    }
//    if constexpr(MaxDegree >=2) {
//        toreturn += (point.transpose() * polyn.coeff2 * point).diagonal();
//    }
//    if constexpr(MaxDegree >=1) {
//        toreturn += polyn.coeff1 * point;
//    }
//    if constexpr(MaxDegree >=0) {
//        toreturn += Eigen::Matrix<Real, 1, Eigen::Dynamic>::Constant(1, nPoints, polyn.coeff0);
//    }
//    return toreturn;
//}
//// Same function as above, but optimised for evaluation at a single point
//template<int NDim, int MaxDegree, typename DerivedPoint>
//Real applySingle(const Polynomial<NDim,MaxDegree>& polyn, const Eigen::MatrixBase<DerivedPoint>& point)
//requires (Eigen::MatrixBase<DerivedPoint>::ColsAtCompileTime == 1) {
//    assert(point.rows() == NDim);
//    Real toreturn = 0.;
//    if constexpr(MaxDegree >=4) {
//        assert(false && "Polynomial Degree > 3 are not implemented");
//    }
//    if constexpr(MaxDegree >=3) {
//        for (int i = 0; i<NDim; ++i) { toreturn += point(i) * (point.dot(polyn.coeff3[i]* point)); }
//    }
//    if constexpr(MaxDegree >=2) {
//        toreturn += point.dot(polyn.coeff2 * point);
//    }
//    if constexpr(MaxDegree >=1) {
//        toreturn += polyn.coeff1.dot(point);
//    }
//    if constexpr(MaxDegree >=0) {
//        toreturn += polyn.coeff0;
//    }
//    return toreturn;
//}
//
//template<int NDim, int MaxDegree, typename DerivedP, typename DerivedV> auto
//findPolynomialCoefficients(const Eigen::MatrixBase<DerivedP>& points, const Eigen::MatrixBase<DerivedV>& values){
//    assert(points.cols() == polynomialCoeffN(NDim, MaxDegree));
//    assert(values.size() == polynomialCoeffN(NDim, MaxDegree));
//    
//    constexpr int nCoef = polynomialCoeffN(NDim, MaxDegree);
//    Eigen::Matrix<Real,nCoef,nCoef> matr;
//    for (int i=0; i<nCoef; ++i){
//        PolynomialCoeff<NDim,MaxDegree> polCoeff(PolynomialCoeff<NDim,MaxDegree>::Zero());
//        polCoeff(0,i) = 1.;
//        for (int j = 0; j < nCoef; ++j) {
//            matr(i,j) = applySingle<NDim,MaxDegree>(polCoeff,points.col(j));
//        }
//    }
//    
//    return matr.solve(values).eval();
//}
//
//template<int NDim, int MaxDegree, typename DerivedP> auto
//findPolynomialCoefficients(const Eigen::MatrixBase<DerivedP>& points, const std::function<Real(Point<NDim>)>& funct){
//    assert(points.cols() == polynomialCoeffN(NDim, MaxDegree));
//   
//    constexpr int nCoef = polynomialCoeffN(NDim, MaxDegree);
//    Eigen::Matrix<Real,nCoef,nCoef> matr;
//    Eigen::Matrix<Real,nCoef,1> values;
//    for (int i=0; i<nCoef; ++i){
//        values(i,0) = funct(points.col(i));
//        PolynomialCoeff<NDim,MaxDegree> polCoeff(PolynomialCoeff<NDim,MaxDegree>::Zero());
//        polCoeff(0,i) = 1.;
//        for (int j = 0; j < nCoef; ++j) {
//            matr(j,i) = applySingle<NDim,MaxDegree>(polCoeff,points.col(j));
//        }
//    }
//    
//    return matr.colPivHouseholderQr().solve(values).transpose().eval();
//}
//
//} // namespace PolynomialImplementation
//
//
////*******************************
//// Main class definition
////*******************************
//
//template<int NDim, int MaxDegree> VectorPoint<NDim> referenceInteporlationPoints() requires (NDim <= 3){
//    VectorPoint<NDim> points;
//    points.resize(Eigen::NoChange, polynomialCoeffN(NDim, MaxDegree));
//    if constexpr (MaxDegree == 0) {
//        points.col(0) = Point<NDim>::Zero();
//        return points;
//    }
//    int counter = 0;
//    for (int i = 0; i < MaxDegree+1; ++i){
//        if constexpr (NDim>1){
//            for (int j = 0; j < MaxDegree + 1 - i; ++j){
//                if constexpr (NDim>2){
//                    for (int k = 0; k < MaxDegree + 1 - i - j; ++k){
//                        points(0,counter) = static_cast<Real>(i)/static_cast<Real>(MaxDegree);
//                        points(1,counter) = static_cast<Real>(j)/static_cast<Real>(MaxDegree);
//                        points(2,counter) = static_cast<Real>(k)/static_cast<Real>(MaxDegree);
//                        ++counter;
//                    }
//                } else {
//                    points(0,counter) = static_cast<Real>(i)/static_cast<Real>(MaxDegree);
//                    points(1,counter) = static_cast<Real>(j)/static_cast<Real>(MaxDegree);
//                    ++counter;
//                }
//            }
//        } else {
//            points(0,counter) = static_cast<Real>(i)/static_cast<Real>(MaxDegree);
//            ++counter;
//        }
//    }
//    return points;
//}
//
//
////=======================================================
//// Polynomials of order 0
////===================
//template<int NDim> class Polynomial<NDim,0> :
//public Features::GroupSpace<Polynomial<NDim,0>>,
//public Features::AsymmetricMathFieldSpaceByValue<Polynomial<NDim,0>,Real>
//{
//public:
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//// Main properties
////~~~~~~~~~~~~~~~~~~
//    Real coeff0;
//        
//    static constexpr int numberCoefficients = polynomialCoeffN(NDim, 0);
//    
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//// Main methods
////~~~~~~~~~~~~~~~~~~
//public:
//    //=======================================================
//    // Constructors
//    //===================
//    Polynomial(): coeff0(0.){};
//    Polynomial(Real value): coeff0(value) {};
//    template<typename Derived> Polynomial(const Eigen::MatrixBase<Derived>& polynCoeff): coeff0(monomialCoeffDegree0<NDim,0>(polynCoeff)){};
//    template<typename DerivedP, typename DerivedV> Polynomial(const Eigen::MatrixBase<DerivedP>& points, const Eigen::MatrixBase<DerivedV>& values): coeff0(values(0,0)) {}
//    template<typename DerivedP> Polynomial(const Eigen::MatrixBase<DerivedP>& points, const std::function<Real(Point<NDim>)>& funct): coeff0(funct(points.col(0))) {}
//    
//    // Generators
//    static Polynomial<NDim,0> random(){ return Polynomial<NDim,0>(Eigen::Matrix<Real,1,1>::Random()(0)); }
//        
//    //********************************
//    //* Arithmetic
//    //********************************
//    Polynomial<NDim,0>& operator+=(const Polynomial<NDim,0>& other){coeff0 += other.coeff0; return *this;};
//    Polynomial<NDim,0>& operator-=(const Polynomial<NDim,0>& other){coeff0 -= other.coeff0; return *this;};
//
//    Polynomial<NDim,0> operator-() const {return Polynomial<NDim,0>(-coeff0);};
//        
//    Polynomial<NDim,0>& operator=(const Real scalar){coeff0 = scalar; return *this;};
//    Polynomial<NDim,0>& operator+=(const Real scalar){coeff0 += scalar; return *this;};
//    Polynomial<NDim,0>& operator-=(const Real scalar){coeff0 -= scalar; return *this;};
//    Polynomial<NDim,0>& operator*=(const Real scalar){coeff0 *= scalar; return *this;};
//    Polynomial<NDim,0>& operator/=(const Real scalar){coeff0 /= scalar; return *this;};
//
//    //********************************
//    //* Evaluation
//    //********************************
//
//    template<typename DerivedPoint> auto operator()(const Eigen::MatrixBase<DerivedPoint>& point) const {return PolynomialImplementation::apply<NDim,0>(*this,point);}
//    template<typename DerivedPoint> auto eval(const Eigen::MatrixBase<DerivedPoint>& point) const {return PolynomialImplementation::apply<NDim,0>(*this,point);}
//    template<typename DerivedPoint> Real evalSingle(const Eigen::MatrixBase<DerivedPoint>& point) const {return PolynomialImplementation::applySingle<NDim,0>(*this,point);};
//    
//    // Substitute the expression X = A y + B into the polynomial
//    // where X is the vector of the variables of the polynomial, A and B are vectors
//    // and y the new scalar variable
//    template<typename DerivedA, typename DerivedB> Polynomial<1,0> substitute(const Eigen::MatrixBase<DerivedA>& a, const Eigen::MatrixBase<DerivedB>& b) const { return Polynomial<1,0>(coeff0);}
//
//};
////********************************
////* I/O
////********************************
//template<int NDim> std::ostream &operator<<(std::ostream &t_os, const Polynomial<NDim,0>& t_poly){
//    t_os << std::setprecision(4) << "P(" + PolynomialImplementation::varNames[0];
//    for (int i = 1; i < NDim; ++i) {t_os << "," + PolynomialImplementation::varNames[i]; }
//    t_os << ") = " << t_poly.coeff0;
//    return t_os;
//}
//
//
////=======================================================
//// Polynomials of order 1
////===================
//template<int NDim> class Polynomial<NDim,1> :
//public Features::GroupSpace<Polynomial<NDim,1>>,
//public Features::AsymmetricMathFieldSpaceByValue<Polynomial<NDim,1>,Real>
//{
//public:
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//// Main properties
////~~~~~~~~~~~~~~~~~~
//    Real coeff0;
//    Eigen::Matrix<Real, 1, NDim> coeff1;
//    
//    static constexpr int numberCoefficients = polynomialCoeffN(NDim, 1);
//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//// Main methods
////~~~~~~~~~~~~~~~~~~
//public:
//    //=======================================================
//    // Constructors
//    //===================
//    Polynomial(): coeff0(0.), coeff1(Eigen::Matrix<Real, 1, NDim>::Zero()){};
//    template<typename Derived> Polynomial(Real value, const Eigen::MatrixBase<Derived>& coeff1val): coeff0(value), coeff1(coeff1val) {};
//    template<typename Derived> Polynomial(Real value, Eigen::MatrixBase<Derived>&& coeff1val): coeff0(value), coeff1(std::move(coeff1val)) {};
//    template<typename Derived> Polynomial(const Eigen::MatrixBase<Derived>& polynCoeff): coeff0(monomialCoeffDegree0<NDim,1>(polynCoeff)), coeff1(monomialCoeffDegree1<NDim,1>(polynCoeff)){};
//    template<typename DerivedP, typename DerivedV> Polynomial(const Eigen::MatrixBase<DerivedP>& points, const Eigen::MatrixBase<DerivedV>& values): Polynomial(PolynomialImplementation::findPolynomialCoefficients<NDim,1>(points, values)) {}
//    template<typename DerivedP> Polynomial(const Eigen::MatrixBase<DerivedP>& points, const std::function<Real(Point<NDim>)>& funct): Polynomial(PolynomialImplementation::findPolynomialCoefficients<NDim,1>(points, funct)) {}
//
//    
//    
//    // Generators
//    static Polynomial<NDim,1> random(){ return Polynomial<NDim,1>(Eigen::Matrix<Real,1,1>::Random()(0), Eigen::Matrix<Real, 1, NDim>::Random()); }
//        
//    //********************************
//    //* Arithmetic
//    //********************************
//    Polynomial<NDim,1>& operator+=(const Polynomial<NDim,1>& other){coeff0 += other.coeff0; coeff1 += other.coeff1; return *this;};
//    Polynomial<NDim,1>& operator-=(const Polynomial<NDim,1>& other){coeff0 -= other.coeff0; coeff1 -= other.coeff1; return *this;};
//
//    Polynomial<NDim,1> operator-() const {return Polynomial<NDim,1>(-coeff0, -coeff1);};
//        
//    Polynomial<NDim,1>& operator=(const Real scalar){coeff0 = scalar; coeff1 = Eigen::Matrix<Real, 1, NDim>::Zero(); return *this;};
//    Polynomial<NDim,1>& operator+=(const Real scalar){coeff0 += scalar; return *this;};
//    Polynomial<NDim,1>& operator-=(const Real scalar){coeff0 -= scalar; return *this;};
//    Polynomial<NDim,1>& operator*=(const Real scalar){coeff0 *= scalar; coeff1 *= scalar; return *this;};
//    Polynomial<NDim,1>& operator/=(const Real scalar){coeff0 /= scalar; coeff1 /= scalar; return *this;};
//
//    //********************************
//    //* Evaluation
//    //********************************
//
//    template<typename DerivedPoint> auto operator()(const Eigen::MatrixBase<DerivedPoint>& point) const {return PolynomialImplementation::apply<NDim,1>(*this,point);}
//    template<typename DerivedPoint> auto eval(const Eigen::MatrixBase<DerivedPoint>& point) const {return PolynomialImplementation::apply<NDim,1>(*this,point);}
//    template<typename DerivedPoint> Real evalSingle(const Eigen::MatrixBase<DerivedPoint>& point) const {return PolynomialImplementation::applySingle<NDim,1>(*this,point);};
//    
//    // Substitute the expression X = A y + B into the polynomial
//    // where X is the vector of the variables of the polynomial, A and B are vectors
//    // and y the new scalar variable
//    template<typename DerivedA, typename DerivedB> Polynomial<1,1> substitute(const Eigen::MatrixBase<DerivedA>& a, const Eigen::MatrixBase<DerivedB>& b) const { return Polynomial<1,1>(coeff0 + coeff1*b, coeff1*a);}
//    
//    //********************************
//    //* Technicalities
//    //********************************
//    // Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
//    // See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
//        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//
//};
////********************************
////* I/O
////********************************
//template<int NDim> std::ostream &operator<<(std::ostream &t_os, const Polynomial<NDim,1>& t_poly){
//    t_os << std::setprecision(4) << "P(" + PolynomialImplementation::varNames[0];
//    for (int i = 1; i < NDim; ++i) {t_os << "," + PolynomialImplementation::varNames[i]; }
//    t_os << ") = " << t_poly.coeff0;
//    for (int i = 0; i <NDim; ++i){
//        if (t_poly.coeff1(0,i) >= 0.) { t_os << " +" << t_poly.coeff1(0,i) << " " << PolynomialImplementation::varNames[i];}
//        else { t_os << " " << t_poly.coeff1(0,i)<< " " << PolynomialImplementation::varNames[i];}
//    }
//    return t_os;
//}
//
//
////=======================================================
//// Polynomials of order 2
////===================
//template<int NDim> class Polynomial<NDim,2> :
//public Features::GroupSpace<Polynomial<NDim,2>>,
//public Features::AsymmetricMathFieldSpaceByValue<Polynomial<NDim,2>,Real>
//{
//public:
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//// Main properties
////~~~~~~~~~~~~~~~~~~
//    Real coeff0;
//    Eigen::Matrix<Real, 1, NDim> coeff1;
//    Eigen::Matrix<Real, NDim, NDim> coeff2;
//    
//    static constexpr int numberCoefficients = polynomialCoeffN(NDim, 2);
//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//// Main methods
////~~~~~~~~~~~~~~~~~~
//public:
//    //=======================================================
//    // Constructors
//    //===================
//    Polynomial(): coeff0(0.), coeff1(Eigen::Matrix<Real, 1, NDim>::Zero()), coeff2 (Eigen::Matrix<Real, NDim, NDim>::Zero()){};
//    template<typename Derived1, typename Derived2>
//        Polynomial(Real value, const Eigen::MatrixBase<Derived1>& coeff1val, const Eigen::MatrixBase<Derived2>& coeff2val): coeff0(value), coeff1(coeff1val), coeff2(coeff2val) {};
//    template<typename Derived1, typename Derived2>
//        Polynomial(Real value, Eigen::MatrixBase<Derived1>&& coeff1val, Eigen::MatrixBase<Derived2>&& coeff2val): coeff0(value), coeff1(std::move(coeff1val)), coeff2(std::move(coeff2val)) {};
//    template<typename Derived> Polynomial(const Eigen::MatrixBase<Derived>& polynCoeff): coeff0(monomialCoeffDegree0<NDim,2>(polynCoeff)), coeff1(monomialCoeffDegree1<NDim,2>(polynCoeff)), coeff2(monomialCoeffDegree2<NDim,2>(polynCoeff)){};
//    template<typename DerivedP, typename DerivedV> Polynomial(const Eigen::MatrixBase<DerivedP>& points, const Eigen::MatrixBase<DerivedV>& values): Polynomial(PolynomialImplementation::findPolynomialCoefficients<NDim,2>(points, values)) {}
//    template<typename DerivedP> Polynomial(const Eigen::MatrixBase<DerivedP>& points, const std::function<Real(Point<NDim>)>& funct): Polynomial(PolynomialImplementation::findPolynomialCoefficients<NDim,2>(points, funct)) {}
//
//    
//    // Generators
//    static Polynomial<NDim,2> random(){
//        Polynomial<NDim,2> pol;
//        PolynomialCoeff<NDim,2> coeff = PolynomialCoeff<NDim,2>::Random();
//        pol.coeff0 = monomialCoeffDegree0<NDim,2>(coeff);
//        pol.coeff1 = monomialCoeffDegree1<NDim,2>(coeff);
//        pol.coeff2 = monomialCoeffDegree2<NDim,2>(coeff);
//        return pol;
//    }
//        
//    //********************************
//    //* Arithmetic
//    //********************************
//    Polynomial<NDim,2>& operator+=(const Polynomial<NDim,2>& other){coeff0 += other.coeff0; coeff1 += other.coeff1; coeff2 += other.coeff2; return *this;};
//    Polynomial<NDim,2>& operator-=(const Polynomial<NDim,2>& other){coeff0 -= other.coeff0; coeff1 -= other.coeff1; coeff2 -= other.coeff2; return *this;};
//
//    Polynomial<NDim,2> operator-() const {return Polynomial<NDim,2>(-coeff0, -coeff1, -coeff2);};
//        
//    Polynomial<NDim,2>& operator=(const Real scalar){coeff0 = scalar; coeff1 = Eigen::Matrix<Real, 1, NDim>::Zero(); coeff2 = Eigen::Matrix<Real, NDim, NDim>::Zero(); return *this;};
//    Polynomial<NDim,2>& operator+=(const Real scalar){coeff0 += scalar; return *this;};
//    Polynomial<NDim,2>& operator-=(const Real scalar){coeff0 -= scalar; return *this;};
//    Polynomial<NDim,2>& operator*=(const Real scalar){coeff0 *= scalar; coeff1 *= scalar; coeff2 *= scalar; return *this;};
//    Polynomial<NDim,2>& operator/=(const Real scalar){coeff0 /= scalar; coeff1 /= scalar; coeff2 /= scalar; return *this;};
//
//    //********************************
//    //* Evaluation
//    //********************************
//
//    template<typename DerivedPoint> auto operator()(const Eigen::MatrixBase<DerivedPoint>& point) const {return PolynomialImplementation::apply<NDim,2>(*this,point);}
//    template<typename DerivedPoint> auto eval(const Eigen::MatrixBase<DerivedPoint>& point) const {return PolynomialImplementation::apply<NDim,2>(*this,point);}
//    template<typename DerivedPoint> Real evalSingle(const Eigen::MatrixBase<DerivedPoint>& point) const {return PolynomialImplementation::applySingle<NDim,2>(*this,point);};
//    
//    // Substitute the expression X = A y + B into the polynomial
//    // where X is the vector of the variables of the polynomial, A and B are vectors
//    // and y the new scalar variable
//    template<typename DerivedA, typename DerivedB> Polynomial<1,2> substitute(const Eigen::MatrixBase<DerivedA>& a, const Eigen::MatrixBase<DerivedB>& b) const {
//        return Polynomial<1,2>(coeff0 + coeff1*b + b.transpose()*coeff2*b, coeff1*a + a.transpose() * coeff2 * b + b.transpose() *coeff2 * a, a.transpose() * coeff2 * a);}
//    
//    //********************************
//    //* Technicalities
//    //********************************
//    // Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
//    // See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
//        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//
//};
////********************************
////* I/O
////********************************
//template<int NDim> std::ostream &operator<<(std::ostream &t_os, const Polynomial<NDim,2>& t_poly){
//    t_os << std::setprecision(4) << "P(" + PolynomialImplementation::varNames[0];
//    for (int i = 1; i < NDim; ++i) {t_os << "," + PolynomialImplementation::varNames[i]; }
//    t_os << ") = " << t_poly.coeff0;
//    for (int i = 0; i <NDim; ++i){
//        if (t_poly.coeff1(0,i) >= 0.) { t_os << " +" << t_poly.coeff1(0,i) << " " << PolynomialImplementation::varNames[i];}
//        else { t_os << " " << t_poly.coeff1(0,i)<< " " << PolynomialImplementation::varNames[i];}
//    }
//    for (int i = 0; i <NDim; ++i){
//        if (t_poly.coeff2(i,i) >= 0.) { t_os << " +" << t_poly.coeff2(i,i) << " " << PolynomialImplementation::varNames[i] << "^2";}
//        else { t_os << " " << t_poly.coeff2(i,i) << " " << PolynomialImplementation::varNames[i] << "^2";}
//        for (int j = i+1; j <NDim; ++j){
//            if (t_poly.coeff2(i,j) >= 0.) { t_os << " +" << t_poly.coeff2(i,j) << " " << PolynomialImplementation::varNames[i] << " " << PolynomialImplementation::varNames[j];}
//            else { t_os << " " << t_poly.coeff2(i,j)<< " " << PolynomialImplementation::varNames[i] << " " << PolynomialImplementation::varNames[j];}
//        }
//    }
//    return t_os;
//}
//
//
////=======================================================
//// Polynomials of order 3
////===================
//template<int NDim> class Polynomial<NDim,3> :
//public Features::GroupSpace<Polynomial<NDim,3>>,
//public Features::AsymmetricMathFieldSpaceByValue<Polynomial<NDim,3>,Real>
//{
//public:
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//// Main properties
////~~~~~~~~~~~~~~~~~~
//    Real coeff0;
//    Eigen::Matrix<Real, 1, NDim> coeff1;
//    Eigen::Matrix<Real, NDim, NDim> coeff2;
//    std::array<Eigen::Matrix<Real, NDim, NDim>, NDim> coeff3;
//    
//    static constexpr int numberCoefficients = polynomialCoeffN(NDim, 3);
//
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//// Main methods
////~~~~~~~~~~~~~~~~~~
//public:
//    //=======================================================
//    // Constructors
//    //===================
//    Polynomial(): coeff0(0.), coeff1(Eigen::Matrix<Real, 1, NDim>::Zero()), coeff2 (Eigen::Matrix<Real, NDim, NDim>::Zero()){
//        for (int i = 0; i < NDim; ++i) { coeff3[i] = Eigen::Matrix<Real, NDim, NDim>::Zero(); }
//    };
//    template<typename Derived1, typename Derived2>
//        Polynomial(Real value, const Eigen::MatrixBase<Derived1>& coeff1val, const Eigen::MatrixBase<Derived2>& coeff2val, const std::array<Eigen::Matrix<Real, NDim, NDim>, NDim>& coeff3val): coeff0(value), coeff1(coeff1val), coeff2(coeff2val), coeff3(coeff3val) {};
//    template<typename Derived1, typename Derived2>
//        Polynomial(Real value, Eigen::MatrixBase<Derived1>&& coeff1val, Eigen::MatrixBase<Derived2>&& coeff2val, std::array<Eigen::Matrix<Real, NDim, NDim>, NDim>&& coeff3val): coeff0(value), coeff1(std::move(coeff1val)), coeff2(std::move(coeff2val)), coeff3(std::move(coeff3val)) {};
//    template<typename Derived> Polynomial(const Eigen::MatrixBase<Derived>& polynCoeff): coeff0(monomialCoeffDegree0<NDim,3>(polynCoeff)), coeff1(monomialCoeffDegree1<NDim,3>(polynCoeff)), coeff2(monomialCoeffDegree2<NDim,3>(polynCoeff)){
//        for (int i=0; i<NDim; ++i) {
//            coeff3[i] = monomialCoeffDegree3<NDim,3>(polynCoeff,i);
//        }
//    };
//    template<typename DerivedP, typename DerivedV> Polynomial(const Eigen::MatrixBase<DerivedP>& points, const Eigen::MatrixBase<DerivedV>& values): Polynomial(PolynomialImplementation::findPolynomialCoefficients<NDim,3>(points, values)) {}
//    template<typename DerivedP> Polynomial(const Eigen::MatrixBase<DerivedP>& points, const std::function<Real(Point<NDim>)>& funct): Polynomial(PolynomialImplementation::findPolynomialCoefficients<NDim,3>(points, funct)) {}
//
//    // Generators
//    static Polynomial<NDim,3> random(){
//        Polynomial<NDim,3> pol;
//        PolynomialCoeff<NDim,3> coeff = PolynomialCoeff<NDim,3>::Random();
//        pol.coeff0 = monomialCoeffDegree0<NDim,3>(coeff);
//        pol.coeff1 = monomialCoeffDegree1<NDim,3>(coeff);
//        pol.coeff2 = monomialCoeffDegree2<NDim,3>(coeff);
//        for (int i = 0; i < NDim; ++i) { pol.coeff3[i] = monomialCoeffDegree3<NDim,3>(coeff,i); }
//        return pol;
//    }
//        
//    //********************************
//    //* Arithmetic
//    //********************************
//    Polynomial<NDim,3>& operator+=(const Polynomial<NDim,3>& other){coeff0 += other.coeff0; coeff1 += other.coeff1; coeff2 += other.coeff2;
//        for (int i = 0; i < NDim; ++i) {coeff3[i] += other.coeff3[i]; } return *this;};
//    Polynomial<NDim,3>& operator-=(const Polynomial<NDim,3>& other){coeff0 -= other.coeff0; coeff1 -= other.coeff1; coeff2 -= other.coeff2;
//        for (int i = 0; i < NDim; ++i) {coeff3[i] -= other.coeff3[i]; } return *this;};
//
//    Polynomial<NDim,3> operator-() const {
//        std::array<Eigen::Matrix<Real, NDim, NDim>, NDim> temp;
//        for (int i = 0; i < NDim; ++i) { temp[i] = -coeff3[i]; }
//        return Polynomial<NDim,3>(-coeff0, -coeff1, -coeff2, temp);};
//        
//    Polynomial<NDim,3>& operator=(const Real scalar){coeff0 = scalar; coeff1 = Eigen::Matrix<Real, 1, NDim>::Zero();
//        coeff2 = Eigen::Matrix<Real, NDim, NDim>::Zero(); for (int i = 0; i < NDim; ++i) { coeff3[i] = Eigen::Matrix<Real, NDim, NDim>::Zero(); } return *this;};
//    Polynomial<NDim,3>& operator+=(const Real scalar){coeff0 += scalar; return *this;};
//    Polynomial<NDim,3>& operator-=(const Real scalar){coeff0 -= scalar; return *this;};
//    Polynomial<NDim,3>& operator*=(const Real scalar){coeff0 *= scalar; coeff1 *= scalar; coeff2 *= scalar; for (int i = 0; i < NDim; ++i) { coeff3[i] *= scalar;} return *this;};
//    Polynomial<NDim,3>& operator/=(const Real scalar){coeff0 /= scalar; coeff1 /= scalar; coeff2 /= scalar; for (int i = 0; i < NDim; ++i) { coeff3[i] /= scalar;} return *this;};
//
//    //********************************
//    //* Evaluation
//    //********************************
//
//    template<typename DerivedPoint> auto operator()(const Eigen::MatrixBase<DerivedPoint>& point) const {return PolynomialImplementation::apply<NDim,3>(*this,point);}
//    template<typename DerivedPoint> auto eval(const Eigen::MatrixBase<DerivedPoint>& point) const {return PolynomialImplementation::apply<NDim,3>(*this,point);}
//    template<typename DerivedPoint> Real evalSingle(const Eigen::MatrixBase<DerivedPoint>& point) const {return PolynomialImplementation::applySingle<NDim,3>(*this,point);};
//    
//    // Substitute the expression X = A y + B into the polynomial
//    // where X is the vector of the variables of the polynomial, A and B are vectors
//    // and y the new scalar variable
//    template<typename DerivedA, typename DerivedB> Polynomial<1,3> substitute(const Eigen::MatrixBase<DerivedA>& a, const Eigen::MatrixBase<DerivedB>& b) const {
//        Eigen::Matrix<Real, 1, 1> coeff0Res = Eigen::Matrix<Real, 1, 1>::Constant(coeff0) + coeff1*b + b.transpose()*coeff2*b;
//        Eigen::Matrix<Real, 1, 1> coeff1Res = coeff1*a + a.transpose() * coeff2 * b + b.transpose() *coeff2 * a;
//        Eigen::Matrix<Real, 1, 1> coeff2Res = a.transpose() * coeff2 * a;
//        Eigen::Matrix<Real, 1, 1> coeff3Res ;
//        for (int i=0; i<NDim; ++i) {
//            Eigen::Matrix<Real, 1, 1> temp1 = b.transpose() * coeff3[i] * b;
//            Eigen::Matrix<Real, 1, 1> temp2 = (b.transpose() *coeff3[i] * a + a.transpose() *coeff3[i] * b);
//            Eigen::Matrix<Real, 1, 1> temp3 = a.transpose() * coeff3[i] * a;
//            coeff0Res += b(i,0) * temp1;
//            coeff1Res += a(i,0) * temp1 + b(i,0) * temp2;
//            coeff2Res += b(i,0) *temp3 + a(i,0) * temp2;
//            coeff2Res += a(i,0) *temp3;
//        }
//        
//        return Polynomial<1,3>(coeff0Res(0,0), coeff1Res, coeff2Res, std::array<Eigen::Matrix<Real, 1, 1>,1>{coeff3Res}) ;}
//    
//    //********************************
//    //* Technicalities
//    //********************************
//    // Required to allow for correct creation of objects containing fixed-size vectorizable Eigen types, as they require 128-bit alignment
//    // See: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
//        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//
//};
////********************************
////* I/O
////********************************
//template<int NDim> std::ostream &operator<<(std::ostream &t_os, const Polynomial<NDim,3>& t_poly){
//    t_os << std::setprecision(4) << "P(" + PolynomialImplementation::varNames[0];
//    for (int i = 1; i < NDim; ++i) {t_os << "," + PolynomialImplementation::varNames[i]; }
//    t_os << ") = " << t_poly.coeff0;
//    for (int i = 0; i <NDim; ++i){
//        if (t_poly.coeff1(0,i) >= 0.) { t_os << " +" << t_poly.coeff1(0,i) << " " << PolynomialImplementation::varNames[i];}
//        else { t_os << " " << t_poly.coeff1(0,i)<< " " << PolynomialImplementation::varNames[i];}
//    }
//    for (int i = 0; i <NDim; ++i){
//        if (t_poly.coeff2(i,i) >= 0.) { t_os << " +" << t_poly.coeff2(i,i) << " " << PolynomialImplementation::varNames[i] << "^2";}
//        else { t_os << " " << t_poly.coeff2(i,i) << " " << PolynomialImplementation::varNames[i] << "^2";}
//        for (int j = i+1; j <NDim; ++j){
//            if (t_poly.coeff2(i,j) >= 0.) { t_os << " +" << t_poly.coeff2(i,j) << " " << PolynomialImplementation::varNames[i] << " " << PolynomialImplementation::varNames[j];}
//            else { t_os << " " << t_poly.coeff2(i,j)<< " " << PolynomialImplementation::varNames[i] << " " << PolynomialImplementation::varNames[j];}
//        }
//    }
//    for (int i = 0; i <NDim; ++i){
//        if (t_poly.coeff3[i](i,i) >= 0.) { t_os << " +" << t_poly.coeff3[i](i,i) << " " << PolynomialImplementation::varNames[i] << "^3";}
//        else { t_os << " " << t_poly.coeff3[i](i,i) << " " << PolynomialImplementation::varNames[i] << "^3";}
//        for (int j = i+1; j <NDim; ++j){
//            if (t_poly.coeff3[i](j,j) >= 0.) { t_os << " +" << t_poly.coeff3[i](i,j) << " " << PolynomialImplementation::varNames[i] << " " << PolynomialImplementation::varNames[j] << "^2";}
//            else { t_os << " " << t_poly.coeff3[i](j,j)<< " " << PolynomialImplementation::varNames[i] << " " << PolynomialImplementation::varNames[j] << "^2";}
//            for (int k = j+1; k <NDim; ++k){
//                if (t_poly.coeff3[i](j,k) >= 0.) { t_os << " +" << t_poly.coeff3[i](i,k) << " " << PolynomialImplementation::varNames[i] << " " << PolynomialImplementation::varNames[j] << " " << PolynomialImplementation::varNames[k];}
//                else { t_os << " " << t_poly.coeff3[i](j,k)<< " " << PolynomialImplementation::varNames[i] << " " << PolynomialImplementation::varNames[j] << " " << PolynomialImplementation::varNames[k];}
//            }
//        }
//    }
//    return t_os;
//}
//
//} // namespace GeometryCore
//
//} // namespace Tortoise
//
//#endif /* Polynomial_hpp */
