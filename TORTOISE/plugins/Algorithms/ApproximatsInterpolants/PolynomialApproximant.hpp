////
////  PolynomialApproximant.hpp
////  TORTOISE
////
////  Created by Marco Battiato on 26/5/23.
////
//
//#ifndef PolynomialApproximant_hpp
//#define PolynomialApproximant_hpp
//
//#include <Eigen/Dense>
//
//template <int order> class PolynomialApproximant;
//
//template <> class PolynomialApproximant<0> {
//// Properties
//    const int                                   nDim;
//    double                                      coeff;
//
//// Methods
//public:
//    int nDimensions() const {
//        return nDim;
//    }
//    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& x){
//        assert(x.rows() == nDim);
//        assert(x.cols() == 1);
//        return coeff;
//    }
//    template <typename Derived> auto nextOrderX(const Eigen::MatrixBase<Derived>& x){
//        return x.array().matrix(); // This is a trick to allow to return x as an expression
//    }
//};
//
//template <int order> class PolynomialApproximant : PolynomialApproximant<order-1> {
//// Properties
//    Eigen::Matrix<double, 1, Eigen::Dynamic>    coeff;
//
//// Methods
//public:
//    int nDimensions() const {
//            return PolynomialApproximant<order-1>::nDimensions();
//        }
//    template <typename Derived>
//    auto operator()(const Eigen::MatrixBase<Derived>& x){
//        return PolynomialApproximant<order-1>::operator()(x) + (coeff * PolynomialApproximant<order-1>::nextOrderX(x)).array();
//    }
//    template <typename Derived> auto nextOrderX(const Eigen::MatrixBase<Derived>& x){
//        return (PolynomialApproximant<order-1>::nextOrderX(x) * x.transpose()).reshaped(Eigen::AutoSize, 1);
//    }
//};
//
//
////class PolynomialApproximant3 {
////
////// Properties
////    const int                                   nDim;
////    double                                      a0;
////    Eigen::Matrix<double, 1, Eigen::Dynamic>    a1;
////    Eigen::Matrix<double, 1, Eigen::Dynamic>    a2;
////    Eigen::Matrix<double, 1, Eigen::Dynamic>    a3;
////
////// Methods
////public:
//////~~~~ Construction/Assignment
////    PolynomialApproximant(int nDim): nDim(nDim) {};
////
//////~~~~ Evaluation
////    template <typename Derived>
////    auto operator()(const Eigen::MatrixBase<Derived>& x){
////
////        assert(x.rows() == nDim);
////        assert(x.cols() == 1);
////
////        return a0 + (a1 * x).array() + (a2 * (x * x.transpose()).reshaped(Eigen::AutoSize, 1)).array() + (a3 * ((x * x.transpose()).reshaped(Eigen::AutoSize, 1) * x.transpose()).reshaped(Eigen::AutoSize, 1)).array();
////
////    }
////
//////~~~~ I/O
////
////};
//
//#endif /* PolynomialApproximant_hpp */
//
//
//
////template <int order> struct PolynomialApproximant;
////
////template <> struct PolynomialApproximant<0> {
////// Properties
////    const int                                   nDim;
////    double                                      coeff;
////// Methods
////public:
////// ~~~ Constructors
////    PolynomialApproximant<0>(int nDim):nDim(nDim){}
////// ~~~ Fitting
////        template <typename Derived1, typename Derived2>
////        void fitData(const Eigen::MatrixBase<Derived1>& x, const Eigen::MatrixBase<Derived2>& y){
////        }
////
////// ~~~ Informations
////    int nDimensions() const {
////        return nDim;
////    }
////    int nCoefficients() const {
////        return 1;
////    }
////// ~~~ Evaluations
////    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& x) const {
////        assert(x.rows() == nDim);
////        assert(x.cols() == 1);
////        return coeff;
////    }
////
////// ~~~ Technicalities (Not to be used by final user)
////    template <typename Derived> auto nextOrderX(const Eigen::MatrixBase<Derived>& x) const {
////        return x.array().matrix(); // This is a trick to allow to return x as an expression
////    }
////    template <typename Derived, typename Derived2>
////    void applyAllBasisFunctions(const Eigen::MatrixBase<Derived>& x, Eigen::MatrixBase<Derived2>& psyx) const {
////        psyx.col(0) = Eigen::Matrix<double,Eigen::Dynamic,1>::Constant(x.cols(), 1, 1.);
////    }
////    template <typename Derived> void extractCoeffValues(const Eigen::MatrixBase<Derived>& coeffToExtract){
////        coeff = coeffToExtract(0,0);
////    }
////};
////
////template <int order> struct PolynomialApproximant : PolynomialApproximant<order-1> {
////// Properties
////    Eigen::Matrix<double, 1, Eigen::Dynamic>    coeff;
////// Methods
////public:
////// ~~~ Constructors
////    PolynomialApproximant(int nDim):PolynomialApproximant<order-1>(nDim){
////        coeff.resize(1, Tortoise::Utilities::intPower<order>(nDim));
////    }
////// ~~~ Fitting
////    template <typename Derived1, typename Derived2>
////    void fitData(const Eigen::MatrixBase<Derived1>& x, const Eigen::MatrixBase<Derived2>& y){
////        assert(x.rows() == nDimensions()); assert(x.cols() == y.cols()); assert(y.rows() == 1);
////
////        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> psyx;
////        psyx.resize(x.cols(), nCoefficients());
////        applyAllBasisFunctions(x,psyx);
////
////        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> psyxTy = psyx.transpose()*y.transpose();
////        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> psyTPsy = psyx.transpose()*psyx;
////
////        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> solution = psyTPsy.llt().solve(psyxTy);
////
////        extractCoeffValues(solution);
////
////        cout << solution << "\n\n";
////    }
////// ~~~ Informations
////    int nDimensions() const {
////            return PolynomialApproximant<order-1>::nDimensions();
////        }
////    int nCoefficients() const {
////        if constexpr (order == 1) return 1 + nDimensions();
////        return Tortoise::Utilities::intPower<1+order>(nDimensions())/(nDimensions()-1);
////    }
////// ~~~ Evaluations
////    template <typename Derived> auto operator()(const Eigen::MatrixBase<Derived>& x) const {
////        return PolynomialApproximant<order-1>::operator()(x) + (coeff * PolynomialApproximant<order-1>::nextOrderX(x)).array();
////    }
////
////// ~~~ Technicalities (Not to be used by final user)
////    template <typename Derived> auto nextOrderX(const Eigen::MatrixBase<Derived>& x) const {
////        return (PolynomialApproximant<order-1>::nextOrderX(x) * x.transpose()).reshaped(Eigen::AutoSize, 1);
////    }
////    template <typename Derived, typename Derived2>
////    void applyAllBasisFunctions(const Eigen::MatrixBase<Derived>& x, Eigen::MatrixBase<Derived2>& psyx) const {
////        const int start = PolynomialApproximant<order-1>::nCoefficients();
////        const int size  = nCoefficients() - PolynomialApproximant<order-1>::nCoefficients();
////        for (int i = 0; i < x.cols(); ++i){
////            psyx.block(i, start, 1, size) = PolynomialApproximant<order-1>::nextOrderX(x.col(i)).transpose();
////        }
////        PolynomialApproximant<order-1>::applyAllBasisFunctions(x, psyx);
////    }
////    template <typename Derived> void extractCoeffValues(const Eigen::MatrixBase<Derived>& coeffToExtract){
////        coeff.transpose() = coeffToExtract.block(PolynomialApproximant<order-1>::nCoefficients(), 0, Tortoise::Utilities::intPower<order>(nDimensions()), 1);
////        PolynomialApproximant<order-1>::extractCoeffValues(coeffToExtract);
////    }
////};
////
////
////int main(int argc, const char * argv[]) {
////
////    const int                nDim = 3;
////    PolynomialApproximant<2> bla(nDim);
////
////    bla.PolynomialApproximant<0>::coeff = 0.;
////    bla.PolynomialApproximant<1>::coeff = Eigen::Matrix<double, 1, Eigen::Dynamic>::Constant(1, nDim, 0.);
////    bla.PolynomialApproximant<2>::coeff = Eigen::Matrix<double, 1, Eigen::Dynamic>::Constant(1, nDim*nDim, 0.);
//////    bla.PolynomialApproximant<3>::coeff = Eigen::Matrix<double, 1, Eigen::Dynamic>::Constant(1, nDim*nDim*nDim, 0.);
////
////    const int nData = 1;
////    Eigen::Matrix<double, Eigen::Dynamic, nData>    x = Eigen::Matrix<double, Eigen::Dynamic, nData>::Random(nDim, nData);
////    cout << x << "\n\n";
////
////    bla.PolynomialApproximant<2>::coeff(0,1) = 1.;
////
////    for (int i=0; i<nData; ++i) cout << bla(x.col(i)) << " ";
////    cout << "\n\n";
////
//////    Eigen::Matrix<double, 1, nData>                 y = Eigen::Matrix<double, 1, nData>::Random(1, nData);
//////
//////    cout << x << "\n\n";
//////    cout << y << "\n\n";
//////
//////    bla.fitData(x, y);
//////
////////    for (int i=0; i<nData; ++i) cout << bla(x.col(i)) << " ";
////////    cout << "\n";
//////
//////    cout << bla.PolynomialApproximant<0>::coeff << "\n";
//////    cout << bla.PolynomialApproximant<1>::coeff << "\n";
//////    cout << bla.PolynomialApproximant<2>::coeff << "\n";
////
////    return 0;
////}
