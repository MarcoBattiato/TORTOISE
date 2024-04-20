// Copyright © 2019, Marco Battiato <marco.battiato@ntu.edu.sg; battiato.marco@gmail.com>, All rights reserved.
//
// Licensed under the GNU GENERAL PUBLIC LICENSE Version 3 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   https://www.gnu.org/licenses/
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
// implied. See the License for the specific language governing
// permissions and limitations under the License.
//
//
// DISCLAIMER: This is a version under active development and testing.
// Not all features have been sufficiently tested, and several features
// are only partially implemented. Moreover both the core and the interface
// may change. You are discouraged from using this version to publish
// results without the supervision of the developer.
//
// If you want, you are welcome to act as a beta tester. In that case
// please contact Marco Battiato at:
// marco.battiato@ntu.edu.sg or battiato.marco@gmail.com
//
// Check if a newer, stable and tested version has, in the meanwhile,
// been made available at:
// https://github.com/MarcoBattiato/TORTOISE
//
//
//  DormandPrince.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 20/10/20.
//
// It calculates the time propagation with adaptive time step Dormand Prince 853
// It needs to work on a container that describes a function of time into a vector space
// The easiest way to ensure that the container has all the necessary methods is to inherit UnivariateRelationContainer
// and ensuring that containerArg in UnivariateRelationContainer is a container of Real type values like for instance double
// and that containerVal is a container of vectors on the containerArg.
// The easiest way to ensure that is that the type containerVal::value_type inherits from VectorSpace
//
// IMPORTANT: remeber to inehrit everything as >> public << otherwise its members might be considered
// as private (if you are defining a class and not a struct)
//
// If the user does not want to inherit from the suggested base classes, the user must ensure that the following type/methods are defined
//  ContainerType must define public:
//      typename    ->  ContainerType::ArgType
//      typename    ->  ContainerType::ValType
//      constructor ->  ContainerType(const ContainerType::ArgType (&), const ContainerType::ValType (&));
//      accessor    ->  ContainerType::ValType     back()
//      accessor    ->  ContainerType::ArgType     timesback();
//      add element ->  (void) emplace_back(const ContainerType::ArgType (&), const ContainerType::ValType (&));
//
//   ContainerType::ArgType must be a continuous scalar type like double and must have a constructor from double
//
//   ContainerType::ValType must define
//      copy constr ->  ContainerType::ValType(const ContainerType::ValType (&));
//      +, - between ContainerType::ValType
//      * by the scalar ContainerType::ArgType
//
//
// =====================================
// Input description
// -> timeContainer: DP853 will start the propagation from the last entry. If the timeContainer is empty if will stop. The calculated dense time steps will be appendend
//                      to this container. This makes this parameter both input and output
// -> finalTime:     DP853 will propagate until the time reaches finalTime. Notice that this is not the time on top of the time of the last entry in timeContainer. If
//                      the time of the last entry in the timeContainer is larger than finalTime, DP853 will do nothing
// -> denseTimeStep: time step of the dense output (see description of the numerical method)
// -> timeOp:        operator that gives the time derivative.
//                      It must be > timeOp(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> ContainerType::ValType
// -> sqrErrorForm:  give the square of the normalised global error given a local error and the y status.
//                      It must be > sqrErrorForm(ContainerType::ValType localErrY, ContainerType::ValType currentY) -> ContainerType::ArgType
// -> outputFunct:   used to produce user defined output. It will attached to the standard output of DP853. It is suggested to keep the user-defined output short
//                      and within a single line for readibility. Also it is better not to add newlines, again to improve readability.
//                      It must be in the form > outputFunct(ContainerType::ArgType currentTime, ContainerType::ValType currentY) -> void
//
// The function returns a timeContainer with all the calculated adaptive steps. Usually these are not important, but are returned for testing reasons
//
// =====================================
// Example
// =====================================
//
//#include "DormandPrince.hpp"
//#include "UnivariateRelationContainer.hpp"
//#include <Eigen/Dense>
//#include <iostream>
//#include <vector>
//
//using namespace Tortoise;
//using std::cout;
// // The following class does little more than UnivariateRelationContainer itself. However it is defined so the user has an example
// // that can be expanded.
//class TimePos2 : public UnivariateRelationContainer<std::vector<double>,std::vector<Eigen::Matrix<double, 2, 1>>> {
//public:
//    TimePos2(const double& time_p, const ValType& vector_p): UnivariateRelationContainer(time_p,vector_p) {};
//    TimePos2(const TimePos2& other): UnivariateRelationContainer(other){};
//    TimePos2(TimePos2&& other): UnivariateRelationContainer(other){};
//    friend std::ostream &operator<<(std::ostream &os, const TimePos2 & funct);
//};
//
//std::ostream &operator<<(std::ostream &os, const TimePos2 & timpos) {
//    std::string output;
//    for (int i=0; i <timpos.size(); ++i){ output+= std::to_string(timpos.arg(i));output+= "  ";output+= std::to_string(timpos[i][0]);output+= "  ";output+= std::to_string(timpos[i][1]);output+= "\n";}
//    return os << output;
//}
//
//Eigen::Matrix<double, 2, 1> timeOp2(const double time, const Eigen::Matrix<double, 2, 1>& stat){return {cos(time),stat[0]};}
//double errorForm2(const Eigen::Matrix<double, 2, 1>& err, const Eigen::Matrix<double, 2, 1>& stat) { return 10000.*(err[0]*err[0]+err[1]*err[1]);}
//void outputFunct2(const double time, const Eigen::Matrix<double, 2, 1>& stat){ std::cout << " t = " << time << " x0 = " << stat[0] << " x1 = " << stat[1] << " ";}
//
//int main(int argc, const char * argv[]) {
//    TimePos2 solution2(0.0,{0.0,0.0});
//    auto blall2 = dormandPrince853<TimePos2>(solution2,1.0,0.1,timeOp2, errorForm2, outputFunct2);
//    std::cout << "\n\n";
//    std::cout << solution2;
//    std::cout << "\n\n";
//    std::cout << sin(1.)<< " "<< 1.0-cos(1.)<< "\n\n\n";
//    return 0;
//}

// TODO
// Remove the output of adaptive step results to increase speed
// Add function to constrain solution at every time step


#ifndef DormandPrince_hpp
#define DormandPrince_hpp

#include <functional>
#include <array>
#include <vector>
#include <iostream>
#include <math.h>

namespace Tortoise {

namespace Algorithms {

template <typename ContainerType, typename timeOpFunctType, typename errorFunctType, typename outputFunctType> ContainerType dormandPrince853
 (ContainerType&                                timeContainer,      // input/output
  const typename ContainerType::ArgType         finalTime,          // end time of propagation
  const typename ContainerType::ArgType         denseTimeStep,      // step for dense grid output
  timeOpFunctType                               timeOp,             // operator that gives the derivative
  errorFunctType                                sqrErrorForm,       // squared normalised error
  outputFunctType                               outputFunct         // Defines some output operations
);


// =======================================
// =======================================
// IMPLEMENTATION
// =======================================
// =======================================

template <typename ContainerType, typename timeOpFunctType, typename errorFunctType, typename outputFunctType> ContainerType dormandPrince853
 (ContainerType&                            timeContainer,                                  // input/output
  const typename ContainerType::ArgType     finalTime,                                      //
  const typename ContainerType::ArgType     denseTimeStep,
  timeOpFunctType                           timeOp,
  errorFunctType                            sqrErrorForm,
  outputFunctType                           outputFunct   // Defines some output operations
  ){
      // Constants for Runge-Kutta
      const std::vector<std::vector<double>> dopr853_a = {
/* 1 */   {},
/* 2 */   {5.26001519587677318785587544488e-2},
/* 3 */   {1.97250569845378994544595329183e-2, 5.91751709536136983633785987549e-2} ,
/* 4 */   {2.95875854768068491816892993775e-2,                                  0,  8.87627564304205475450678981324e-2} ,
/* 5 */   {2.41365134159266685502369798665e-1,                                  0, -8.84549479328286085344864962717e-1, 9.24834003261792003115737966543e-1},
/* 6 */   {3.7037037037037037037037037037e-2,                                   0,                                   0, 1.70828608729473871279604482173e-1, 1.25467687566822425016691814123e-1},
/* 7 */   {3.7109375e-2, 0, 0, 1.70252211019544039314978060272e-1, 6.02165389804559606850219397283e-2, -1.7578125e-2},
/* 8 */   {3.70920001185047927108779319836e-2, 0, 0, 1.70383925712239993810214054705e-1, 1.07262030446373284651809199168e-1, -1.53194377486244017527936158236e-2, 8.27378916381402288758473766002e-3},
/* 9 */   {6.24110958716075717114429577812e-1, 0, 0, -3.36089262944694129406857109825e+0,
          -8.68219346841726006818189891453e-1, 2.75920996994467083049415600797e+1, 2.01540675504778934086186788979e+1, -4.34898841810699588477366255144e+1},
/* 10 */  {4.77662536438264365890433908527e-1, 0, 0, -2.48811461997166764192642586468e+0, -5.90290826836842996371446475743e-1, 2.12300514481811942347288949897e+1,
           1.52792336328824235832596922938e+1, -3.32882109689848629194453265587e+1, -2.03312017085086261358222928593e-2},
/* 11 */  {-9.3714243008598732571704021658e-1, 0, 0,
           5.18637242884406370830023853209e+0, 1.09143734899672957818500254654e+0, -8.14978701074692612513997267357e+0, -1.85200656599969598641566180701e+1,
           2.27394870993505042818970056734e+1, 2.49360555267965238987089396762e+0, -3.0467644718982195003823669022e+0},
/* 12 */  {2.27331014751653820792359768449e+0, 0, 0,
           -1.05344954667372501984066689879e+1, -2.00087205822486249909675718444e+0, -1.79589318631187989172765950534e+1, 2.79488845294199600508499808837e+1,
           -2.85899827713502369474065508674e+0, -8.87285693353062954433549289258e+0, 1.23605671757943030647266201528e+1, 6.43392746015763530355970484046e-1},
/* 13 */   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
/* 14 */   {5.61675022830479523392909219681e-2, 0, 0, 0, 0, 0, 2.53500210216624811088794765333e-1, -2.46239037470802489917441475441e-1, -1.24191423263816360469010140626e-1, 1.5329179827876569731206322685e-1, 8.20105229563468988491666602057e-3, 7.56789766054569976138603589584e-3, -8.298e-3},
/* 15 */   {3.18346481635021405060768473261e-2, 0, 0, 0, 0, 2.83009096723667755288322961402e-2, 5.35419883074385676223797384372e-2, -5.49237485713909884646569340306e-2, 0, 0, -1.08347328697249322858509316994e-4, 3.82571090835658412954920192323e-4, -3.40465008687404560802977114492e-4, 1.41312443674632500278074618366e-1},
/* 16 */   {-4.28896301583791923408573538692e-1, 0, 0, 0, 0, -4.69762141536116384314449447206e+0, 7.68342119606259904184240953878e+0, 4.06898981839711007970213554331e+0,
            3.56727187455281109270669543021e-1, 0, 0, 0, -1.39902416515901462129418009734e-3, 2.9475147891527723389556272149e+0, -9.15095847217987001081870187138e+0}};
      const std::vector<double> dopr853_b = {5.42937341165687622380535766363e-2, 0, 0, 0, 0, 4.45031289275240888144113950566e+0,1.89151789931450038304281599044e+0
      , -5.8012039600105847814672114227e+0, 3.1116436695781989440891606237e-1, -1.52160949662516078556178806805e-1, 2.01365400804030348374776537501e-1,
       4.47106157277725905176885569043e-2};
      const std::vector<double> dopr853_c = {0.0, 0.526001519587677318785587544488e-01, 0.789002279381515978178381316732e-01, 0.118350341907227396726757197510e+00,  0.281649658092772603273242802490e+00, 0.333333333333333333333333333333e+00, 0.25e+00, 0.307692307692307692307692307692e+00, 0.651282051282051282051282051282e+00, 0.6e+00, 0.857142857142857142857142857142e+00, 0.0};
      const std::vector<double> dopr853_cdense = {0.1e+00, 0.2e+00, 0.777777777777777777777777777778e+00};
      const std::vector<double> dopr853_bhh = {0.244094488188976377952755905512e+00, 0.733846688281611857341361741547e+00, 0.220588235294117647058823529412e-01};
      const std::vector<double> dopr853_er = {0.1312004499419488073250102996e-01, 0, 0, 0, 0, -0.1225156446376204440720569753e+01,
      -0.4957589496572501915214079952e+00, 0.1664377182454986536961530415e+01, -0.3503288487499736816886487290e+00, 0.3341791187130174790297318841e+00,
       0.8192320648511571246570742613e-01, -0.2235530786388629525884427845e-01};

      const double d41 = -0.84289382761090128651353491142e+1,d46 = 0.56671495351937776962531783590e+0,d47 = -0.30689499459498916912797304727e+1,d48 = 0.23846676565120698287728149680e+1,d49 = 0.21170345824450282767155149946e+1,d410 = -0.87139158377797299206789907490e+0,d411 = 0.22404374302607882758541771650e+1,d412 = 0.63157877876946881815570249290e+0,d413 = -0.88990336451333310820698117400e-1,d414 = 0.18148505520854727256656404962e+2,d415 = -0.91946323924783554000451984436e+1,d416 = -0.44360363875948939664310572000e+1;

      const double d51 = 0.10427508642579134603413151009e+2,d56 = 0.24228349177525818288430175319e+3,d57 = 0.16520045171727028198505394887e+3,d58 = -0.37454675472269020279518312152e+3,d59 = -0.22113666853125306036270938578e+2,d510 = 0.77334326684722638389603898808e+1,d511 = -0.30674084731089398182061213626e+2,d512 = -0.93321305264302278729567221706e+1,d513 = 0.15697238121770843886131091075e+2,d514 = -0.31139403219565177677282850411e+2,d515 = -0.93529243588444783865713862664e+1,d516 = 0.35816841486394083752465898540e+2;

      const double d61 = 0.19985053242002433820987653617e+2,d66 = -0.38703730874935176555105901742e+3,d67 = -0.18917813819516756882830838328e+3,d68 = 0.52780815920542364900561016686e+3,d69 = -0.11573902539959630126141871134e+2,d610 = 0.68812326946963000169666922661e+1,d611 = -0.10006050966910838403183860980e+1,d612 = 0.77771377980534432092869265740e+0,d613 = -0.27782057523535084065932004339e+1,d614 = -0.60196695231264120758267380846e+2,d615 = 0.84320405506677161018159903784e+2,d616 = 0.11992291136182789328035130030e+2;

      const double d71 = -0.25693933462703749003312586129e+2,d76 = -0.15418974869023643374053993627e+3,d77 = -0.23152937917604549567536039109e+3,d78 = 0.35763911791061412378285349910e+3,d79 = 0.93405324183624310003907691704e+2,d710 = -0.37458323136451633156875139351e+2,d711 = 0.10409964950896230045147246184e+3,d712 = 0.29840293426660503123344363579e+2,d713 = -0.43533456590011143754432175058e+2,d714 = 0.96324553959188282948394950600e+2,d715 = -0.39177261675615439165231486172e+2,d716 = -0.14972683625798562581422125276e+3;

      if(timeContainer.empty()) { return timeContainer;}                             // If the container is empty there is no usable starting point so it quits

      typedef typename ContainerType::ArgType  RealType;       // Short name for the Real type (for instance double)
      typedef typename ContainerType::ValType  vecType;        // Short name for the vector type (for instance

      int                       steps           = 0;
      RealType                  step_siz        = denseTimeStep;        // As initial guess for the time step uses the dense step
      RealType                  step_siz_old;
      RealType                  step_siz_min    = std::pow(10, -10);

      RealType                  curr_time       = timeContainer.argback();
      RealType                  dense_time      = curr_time + denseTimeStep;

      ContainerType             adaptiveSteps (timeContainer.argback(), timeContainer.back());

      vecType                   next_t (timeContainer.back());  // This initialisation is needed!!!!

      vecType                   y_err3 (next_t);    // The initialisation is not necessary, but it is chosen since vecType might not have a default constructor, but should have a copy constructor
      vecType                   y_err5 (next_t);    // The initialisation is not necessary, but it is chosen since vecType might not have a default constructor, but should have a copy constructor
      vecType                   y_rk8  (next_t);    // The initialisation is not necessary, but it is chosen since vecType might not have a default constructor, but should have a copy constructor
      vecType                   temp_y (next_t);    // The initialisation is not necessary, but it is chosen since vecType might not have a default constructor, but should have a copy constructor
      vecType                   y_next (next_t);    // The initialisation is not necessary, but it is chosen since vecType might not have a default constructor, but should have a copy constructor
      vecType                   dydx_new (next_t);  // The initialisation is not necessary, but it is chosen since vecType might not have a default constructor, but should have a copy constructor

      std::vector<vecType>    temp(12, next_t);       // The initialisation is not necessary, but it is chosen since vecType might not have a default constructor, but should have a copy constructor
      std::vector<vecType>     rcont(8, next_t);      // The initialisation is not necessary, but it is chosen since vecType might not have a default constructor, but should have a copy constructor

      const RealType safe_fac = 0.9, minscale = 0.333333, maxscale = 6.0, beta=0.04, alpha = 1.0/8.0-beta*0.2;
      RealType scale;

      RealType                  err_old = 1.0e-4;
      bool                      reject = false;
      bool                      notfailed = true;

      while (curr_time < finalTime && notfailed){
          std::cout << "** Time: " << std::fixed << curr_time << " timestep: " << step_siz;

          for(int j=0; j<12; j++){
              temp_y = next_t;
              for(int k=0; k<j; k++){ temp_y += dopr853_a[j][k] * step_siz * temp[k]; }
              temp[j] = timeOp(curr_time + dopr853_c[j] * step_siz, temp_y);
//              temp[j] =timeOp(curr_time+dopr853_c[j]*step_siz, toExecute(init_time+dopr853_c[j]*step_siz, temp_y));
          }
          
//          for(int j=0; j<12; j++){ std::cout << temp[j] << "|";}
          
          y_rk8 *= 0.0;             // y_rk8 needs to be initialised to 0
          for(int j=0; j<12; j++){ if(dopr853_b[j] != 0) { y_rk8 += dopr853_b[j]*temp[j];}  }
          y_err3 = y_rk8 - dopr853_bhh[0] * temp[0] - dopr853_bhh[1] * temp[8] - dopr853_bhh[2] * temp[11];

          y_next = next_t + (step_siz * y_rk8);
//          y_next= toExecute(curr_time+step_siz, next_t+(step_siz*y_rk8));

          y_err5 *= 0.0;            // y_err5 needs to be initialised to 0
          for(int j=0; j<12; j++){  if(dopr853_er[j] != 0)  {  y_err5 += dopr853_er[j]*temp[j]; } }

//          std::cout << " " << y_next << " " << y_err3 << " " << y_err5 << " ";
          
          // Evaluation of the error
          double errbnd5 = sqrErrorForm(y_err5,y_next);
          double deno = errbnd5 + 0.01 * sqrErrorForm(y_err3,y_next);
          if(deno<=0.0){ deno=1.0;}
          double err = errbnd5 / std::sqrt(deno);

          std::cout << " error: " << err ;

          // Checking error
          if(err <= 1.){
              // Case when error is acceptable

              // Adjustment of step size
              step_siz_old = step_siz;                                                  // The current step size is saved as old
              // The step size will then be modified by multiplying it by a scale factor
              if(err == 0){scale = maxscale;}                                           // If the error is 0, the scale is the maximum one
              else{
                  scale = (safe_fac * std::pow(err_old,beta) * std::pow(err,-alpha));   // This formula predicts how much the scale should be, but it adds an inertia (it uses the old error as well)
                  scale = std::min(std::max(scale, minscale), maxscale);                // The scale is however capped between min and max allowed values
              }
              if(reject) { step_siz *= std::min(1.0, scale);}                           // If the previous step was rejected it will not increase the step size
              else { step_siz *= scale;}
              err_old = std::max(err, 1.0e-4);
              reject = false;                                                           // As the current step was accepted, it will store this info

              // output
              std::cout << " ACCEPTED <<< ";
              curr_time += step_siz_old;
              outputFunct(curr_time, y_next);        // user-defined output
              std::cout << "\n";

              if(curr_time >= dense_time){
                  // Preparation for Dense Output
                  dydx_new = timeOp(curr_time, y_next);
                  rcont[0] = next_t;
                  rcont[1] = y_next - next_t;
                  rcont[2] = step_siz_old * temp[0] - rcont[1];
                  rcont[3] = rcont[1] - step_siz_old * dydx_new - rcont[2];
                  rcont[4] = d41 * temp[0] + d46 * temp[5] + d47 * temp[6] + d48 * temp[7] + d49 * temp[8] + d410 * temp[9]+ d411 * temp[10]+ d412 * temp[11];
                  rcont[5] = d51 * temp[0] + d56 * temp[5] + d57 * temp[6] + d58 * temp[7] + d59 * temp[8] + d510 * temp[9]+ d511 * temp[10]+ d512 * temp[11];
                  rcont[6] = d61 * temp[0] + d66 * temp[5] + d67 * temp[6] + d68 * temp[7] + d69 * temp[8] + d610 * temp[9]+ d611 * temp[10]+ d612 * temp[11];
                  rcont[7] = d71 * temp[0] + d76 * temp[5] + d77 * temp[6] + d78 * temp[7] + d79 * temp[8] + d710 * temp[9]+ d711 * temp[10]+ d712 * temp[11];

                  temp[9] = timeOp(curr_time - step_siz_old + dopr853_cdense[0] * step_siz_old,
                                 next_t + step_siz_old * ( dopr853_a[13][0] * temp[0] + dopr853_a[13][6] * temp[6] + dopr853_a[13][7] * temp[7] + dopr853_a[13][8] * temp[8] + dopr853_a[13][9] * temp[9] + dopr853_a[13][10] * temp[10] + dopr853_a[13][11] * temp[11] + dopr853_a[13][12] * dydx_new) );


                  temp[10] = timeOp(curr_time - step_siz_old + dopr853_cdense[1] * step_siz_old,
                                  next_t + step_siz_old * ( dopr853_a[14][0] * temp[0] + dopr853_a[14][5] * temp[5] + dopr853_a[14][6] * temp[6] + dopr853_a[14][7] * temp[7] + dopr853_a[14][10] * temp[10] + dopr853_a[14][11] * temp[11] + dopr853_a[14][12] * dydx_new + dopr853_a[14][13] * temp[9]));

                  temp[11] = timeOp(curr_time - step_siz_old + dopr853_cdense[2] * step_siz_old,
                                  next_t + step_siz_old * ( dopr853_a[15][0] * temp[0] + dopr853_a[15][5] * temp[5] + dopr853_a[15][6] * temp[6] + dopr853_a[15][7] * temp[7] + dopr853_a[15][8] * temp[8] + dopr853_a[15][12] * dydx_new + dopr853_a[15][13] * temp[9] + dopr853_a[15][14] * temp[10]));

                  rcont[4] = step_siz_old * ( rcont[4] + d413 * dydx_new + d414 * temp[9] + d415 * temp[10] + d416 * temp[11] );
                  rcont[5] = step_siz_old * ( rcont[5] + d513 * dydx_new + d514 * temp[9] + d515 * temp[10] + d516 * temp[11] );
                  rcont[6] = step_siz_old * ( rcont[6] + d613 * dydx_new + d614 * temp[9] + d615 * temp[10] + d616 * temp[11] );
                  rcont[7] = step_siz_old * ( rcont[7] + d713 * dydx_new + d714 * temp[9] + d715 * temp[10] + d716 * temp[11] );

                  while(curr_time > dense_time && dense_time <= finalTime){

                      double s = ( dense_time - ( curr_time - step_siz_old ) ) / step_siz_old;
                      double s1 = 1.0 - s;
                      timeContainer.emplace_back(dense_time, rcont[0]+s*(rcont[1]+s1*(rcont[2]+s*(rcont[3]+s1*(rcont[4]+s*(rcont[5]+s1*(rcont[6]+s*rcont[7])))))) );
                      dense_time += denseTimeStep;
                  }
                  next_t=y_next;
                  adaptiveSteps.emplace_back(curr_time,y_next);
              }
              else{
                  next_t=y_next;
                  adaptiveSteps.emplace_back(curr_time,y_next);
              }
              steps++;
          }
          else{
              scale = std::max(safe_fac*std::pow(err,-alpha), minscale);
              step_siz *= scale;
              reject = true;

              // output
              std::cout << " REJECTED !!! ";
              outputFunct(curr_time, y_next);        // user-defined output
              std::cout << "\n";
              if(step_siz < step_siz_min){
                  std::cout << "The step size has become too small to be handled" << std::endl;
                  notfailed=false;
              }
          }
      }
      return adaptiveSteps;
  }

} // namespace Algorithms 

} // namespace Tortoise

#endif /* DormandPrince_hpp */
