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
//  Algorithms.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 4/7/20.
//

#ifndef NelderMead_hpp
#define NelderMead_hpp

#include <array>
#include <functional>
#include <assert.h>
#include <iostream>
#include <math.h>

// Minimization of a given function using Nelder-Mead algorithm
//
// The method constructs a simplex
// Implemented following http://www.scholarpedia.org/article/Nelder-Mead_algorithm#The_Nelder-Mead_simplex_algorithm
//
// After each iteration a convergence test is executed
// Two tests are performed
// - Domain convergence: the sizes of the simplex across each direction are compared with the provided tolerances
// - Function-value convergence: this test is executed in two different ways. If the minimum value is not known then the spread of the function values across the vertices of the simplex is compared with the given tolerance. Notice that this does not guarantee that the function is at a minimum (see later). If the minimum value is known (for whatever reason, for instance when minimising the squared error) and is 0 (if not the function can be shifted to have a 0 minimum), then the user should specify minimizeToZero = true. In that case the values of the function at the simplex's vertices are compared with 0.
// The method stops if either of the convergence test are passed (i.e. as long as one is passed, the method stops).
// If the user wants to use a single test, this can be easily obtained by setting the tolerance for the other test to 0.0 (that test will always fail, and the method will stop only when the wanted condition is satisfied)
//
// The method implemented is greedy minimization by default (called by for instance minimizeNelderMead<12>(...) or equivalently  minimizeNelderMead<12,false>(...) )
// Sometimes the simplex becomes extremely small prematurely slowing down the convergence
// The problem can be mitigated by using greedy expansion (called by for instance minimizeNelderMead<12,true>(...) )
//
// However the best way to obtain good results is to have an initial guess and initial steps reasonably close to the correct solution
//
// The function returns
// 0 : failed convergence
// 1 : domain convergence
// 2 : function-value convergence
//
// It saves the solution (best performing vertex of the simplex) in solution
// In solutionVar it saves the variance (resolved for every direction
// The minimum value obtained is saved at valueAtSolution and its variance in valueVarAtSolution


// EXAMPLE
//
//#include "Algorithms.hpp"
//#include "Clock.hpp"
//
//#include <iostream>
//
//int main(int argc, const char * argv[]) {
//
//    std::cout << "Start\n";
//
//    Clock clock;
//
//    const int nDims = 350;
//    std::function<double(std::array<double, nDims>)> fun = [nDims](std::array<double, nDims> x){
//        double toreturn = 0.0;
//        for (int i=0; i<nDims; ++i) {
//            toreturn += (x[i]-static_cast<double>(i))*(x[i]-static_cast<double>(i));
//        }
//        return toreturn;};
//
//    std::array<double, nDims>  solution;
//    std::array<double, nDims>  solutionVar;
//    double valueAtSolution;
//    double valueVarAtSolution;
//    int executedIterations;
//
//    std::array<double, nDims> start, startStepSizes, variablesTolerances;
//    for (int i=0; i<nDims; ++i) { start[i]=static_cast<double>(i)-0.9; startStepSizes[i]=5.0; variablesTolerances[i]= 0.000;}
//
//    clock.Tic();
//
//    int convergencetype = minimizeNelderMead<nDims,true>(fun, start, startStepSizes, variablesTolerances, 0.1, true, 500000, solution, solutionVar, valueAtSolution, valueVarAtSolution, executedIterations);
//    std::cout << "Time taken " << clock.Toc() << " s\n";
//    std::cout << "Convergence of type " << convergencetype << " in " << executedIterations<< " iterations\n";
//    visualizeNelderMeadResult<nDims>(solution, solutionVar, valueAtSolution, valueVarAtSolution);
//
//    return 0;
//}

namespace Tortoise {

template <int nParameters, bool greedyExpansion = false> int minimizeNelderMead(
        const std::function<double(std::array<double, nParameters>)> &fun,  // Function to minimise
        const std::array<double, nParameters> &start,                       // Starting point
        const std::array<double, nParameters> &startStepSizes,              // Step sizes on each vaariable for first iteration
        const std::array<double, nParameters> &variablesTolerances,         // Tolerances required on simplex to have domain convergence
        double functionValueTolerance,                                      // Tolerance on the function value to have function-value convergence
        bool minimizeToZero,                                                // Set to true if the known minimum is zero
        int maxIterations,                                                  // Max number of iterations
        std::array<double, nParameters> &solution,                          // OUTPUT: solution
        std::array<double, nParameters> &solutionVar,                       // OUTPUT: variance of the solution
        double &valueAtSolution,                                            // OUTPUT: function value at solution
        double &valueVarAtSolution,                                         // OUTPUT: variance of function value at solution
        int &executedIterations );                                          // OUTPUT: number of iterations that have been executed


// Plots to terminal with the following format
// ====== SOLUTION =======
// [0.686798 ± 0.075531, -0.0112915 ± 0.0181274]  -> 0.00272425 ± 0.00130491
template <int nParameters> void visualizeNelderMeadResult(std::array<double, nParameters> &solution, std::array<double, nParameters> &solutionVar, double &valueAtSolution, double &valueVarAtSolution) ;








// =======================================
// =======================================
// IMPLEMENTATION
// =======================================
// =======================================


// Functions showing the status of the data structure (Used for debug)
template <int nParameters> void viewPointNelderMead(std::array<double, nParameters> const &point, double valueAtPoint){
    std::cout << "[";
    for (int dim =0; dim <nParameters-1; ++dim){
        std::cout << point[dim] << ", ";
    }
    std::cout << point[nParameters-1] << "]  -> " << valueAtPoint << "\n";
}
template <int nParameters> void viewStatusNelderMead(std::array<std::array<double, nParameters>, nParameters+1> const &simplex, std::array<double, nParameters+1> const &valuesAtVertices, int worst, int secondWorst, int best, double worstValue, double secondWorstValue, double bestValue){
    std::cout << "====== STATUS ========\n";
    for (int vertex=0; vertex < nParameters+1; ++vertex){
        std::cout << vertex << ": ";
        viewPointNelderMead<nParameters>(simplex[vertex], valuesAtVertices[vertex]);
    }
    std::cout << "Worst        : "<< worst << " -> " << worstValue <<"\n";
    std::cout << "Second Worst : "<< secondWorst << " -> " << secondWorstValue <<"\n";
    std::cout << "Best         : "<< best << " -> " << bestValue <<"\n";
    std::cout << "======================\n";
}


// Convergence test of a single step. Returns
// 0 : not yet converged
// 1 : domain convergence
// 2 : function-value convergence
// To speed up the convergence test, a simplified version has been implemented. Even when the present convergence test is passed, a proper convergence test might be off by a factor of at most 2.
template <int nParameters> int convergenceTestNelderMead (const std::array<std::array<double, nParameters>, nParameters+1> &simplex, const std::array<double, nParameters+1> &valuesAtVertices, const std::array<double, nParameters> &variablesTolerances, double functionValueTolerance, bool minimizeToZero) {
    bool domainConvergence = true;
    bool functionConvergence = true;
    if (minimizeToZero) {
        for (int vertex=0; vertex < nParameters+1; ++vertex) {
            if (valuesAtVertices[vertex] > functionValueTolerance ) {
                functionConvergence = false; break;
            }
        }
    } else {
        for (int vertex=1; vertex < nParameters+1; ++vertex) {
            if (fabs(valuesAtVertices[0]-valuesAtVertices[vertex]) > functionValueTolerance ) {
                functionConvergence = false; break;
            }
        }
    }
    for (int vertex=1; vertex < nParameters+1; ++vertex) {
        for (int dim = 0; dim < nParameters; ++dim) {
            if (fabs(simplex[0][dim]-simplex[vertex][dim]) > variablesTolerances[dim] ) {
                domainConvergence = false; break;
            }
        }
    }
    if (domainConvergence) return 1;
    if (functionConvergence) return 2;
    return 0;
}

template <int nParameters> void fullOrderingNelderMead(const std::array<double, nParameters+1> &valuesAtVertices, int &worst, int &secondWorst, int &best,  double &worstValue, double &secondWorstValue, double &bestValue){
    if (valuesAtVertices[0] > valuesAtVertices[1]){
        worstValue = valuesAtVertices[0];        worst = 0;
        secondWorstValue = valuesAtVertices[1];  secondWorst = 1;
        bestValue = valuesAtVertices[1];         best = 1;
    } else {
        worstValue = valuesAtVertices[1];        worst = 1;
        secondWorstValue = valuesAtVertices[0];  secondWorst = 0;
        bestValue = valuesAtVertices[0];         best = 0;
    }
    for (int vert=2; vert<nParameters+1; ++vert) {
        if (valuesAtVertices[vert] >= worstValue) {
            secondWorstValue = worstValue;              secondWorst = worst;
            worstValue = valuesAtVertices[vert];        worst = vert;
        } else if (valuesAtVertices[vert] >= secondWorstValue) {
            secondWorstValue = valuesAtVertices[vert];  secondWorst = vert;
        } else if (valuesAtVertices[vert] < bestValue) {
            bestValue = valuesAtVertices[vert];         best = vert;
        }
    }
}


template <int nParameters> void visualizeNelderMeadResult(std::array<double, nParameters> &solution, std::array<double, nParameters> &solutionVar, double &valueAtSolution, double &valueVarAtSolution) {
    std::cout << "====== SOLUTION =======\n";
    std::cout << "[";
    for (int dim =0; dim <nParameters-1; ++dim){
        std::cout << solution[dim] << " ± " << solutionVar[dim] << ", ";
    }
    std::cout << solution[nParameters-1] << " ± " << solutionVar[nParameters-1] << "]  ->  " << valueAtSolution << " ± " << valueVarAtSolution << "\n";
}

template <int nParameters, bool greedyExpansion> int minimizeNelderMead(
        const std::function<double(std::array<double, nParameters>)> &fun,  // Function to minimise
        const std::array<double, nParameters> &start,                       // Starting point
        const std::array<double, nParameters> &startStepSizes,              // Step sizes on each vaariable for first iteration
        const std::array<double, nParameters> &variablesTolerances,         // Tolerances required on simplex to have domain convergence
        double functionValueTolerance,                                      // Tolerance on the function value to have function-value convergence
        bool minimizeToZero,                                                // Set to true if the known minimum is zero
        int maxIterations,                                                  // Max number of iterations
        std::array<double, nParameters> &solution,                          // OUTPUT: solution
        std::array<double, nParameters> &solutionVar,                       // OUTPUT: variance of the solution
        double &valueAtSolution,                                            // OUTPUT: function value at solution
        double &valueVarAtSolution,                                         // OUTPUT: variance of function value at solution
        int &executedIterations) {
    
    // Implemented following http://www.scholarpedia.org/article/Nelder-Mead_algorithm#The_Nelder-Mead_simplex_algorithm
    
    assert(nParameters>1);
    
    // Parameters for simplex transformation
    const double alpha = 1.0, beta = 0.5, gamma = 2.0, delta = 0.5;
    
    // Variables for iterative method
    std::array<std::array<double, nParameters>, nParameters+1>  simplex;
    std::array<double, nParameters+1>                           valuesAtVertices;
    int                                                         worst, secondWorst, best;
    double                                                      worstValue, secondWorstValue, bestValue;
    std::array<double, nParameters>                             centroid, reflectPoint, expandPoint, contractPoint;
    double                                                      reflectPointValue, expandPointValue, contractPointValue;
    int                                                         convergenceStatus;
    
    // Constructs the initial simplex as a right-angled simplex at start
    // and with sides parallel to axis and startStepSizes long
    simplex[0] = start;
    for (int side=0; side<nParameters; ++side){
        simplex[side+1] = start;
        simplex[side+1][side] += startStepSizes[side];
    }
    
    //=====================
    // ORDERING:
    //=====================
    for (int vert=0; vert<nParameters+1; ++vert){ valuesAtVertices[vert] = fun(simplex[vert]);} // Evaluates the function at each vertex of the simplex
    fullOrderingNelderMead<nParameters>(valuesAtVertices, worst, secondWorst, best, worstValue, secondWorstValue, bestValue);

    executedIterations=0;
    while (maxIterations--) {                   // Iterates at most maxIterations times
        
        if ( (convergenceStatus = convergenceTestNelderMead<nParameters>(simplex, valuesAtVertices, variablesTolerances,functionValueTolerance,minimizeToZero)) ) {
            break;
        }
        
        executedIterations++;
        
//        std::cout << bestValue << "\n";
//        std::cout << "CONVERGENCE RESULT: " << convergenceStatus << "\n";
        
        //=====================
        // CENTROID:
        //=====================
        for (int dim=0; dim < nParameters; ++dim){ centroid[dim] = - simplex[worst][dim]; }
        for (int vert=0; vert<nParameters+1; ++vert){
            for (int dim=0; dim < nParameters; ++dim){
                centroid[dim] += simplex[vert][dim];
            }
        }
        for (int dim=0; dim < nParameters; ++dim){ centroid[dim] /= nParameters; }
                
//************
//        std::cout << " >>> BEGIN ITERATION <<<\n";
//        viewStatusNelderMead<nParameters>(simplex, valuesAtVertices, worst, secondWorst, best, worstValue, secondWorstValue, bestValue);
//        std::cout << "CENTROID  "; viewPointNelderMead<nParameters>(centroid,0.0);
//************
        
        
        //=====================
        // TRANSFORMATION:
        //=====================
        
        // >>>> REFLECT:
        // >>>>>>>>>>>>
        for (int dim=0; dim<nParameters; ++dim){ reflectPoint[dim] = centroid[dim]+alpha*(centroid[dim] - simplex[worst][dim]);}  // Calculates reflection point
        reflectPointValue = fun(reflectPoint);                                                                                      // Calculates value at reflection point
        
//************
//        std::cout << "REFLECT "; viewPointNelderMead<nParameters>(reflectPoint,reflectPointValue);
//************
        
        if (bestValue <= reflectPointValue && reflectPointValue < secondWorstValue) {
            simplex[worst] = reflectPoint;                                                                                          // Accepts reflection point
            valuesAtVertices[worst] = reflectPointValue;
            // -------
            // Optimised Ordering
            //  reflectPointValue >= bestValue   => No need to update best
            //  reflectPointValue < secondWorstValue  =>  secondWorstValue becomes worstValue & we need to find the new secondWorst
            worst = secondWorst;    worstValue = secondWorstValue;          // secondWorstValue becomes worstValue
            secondWorst = 0;     secondWorstValue = valuesAtVertices[0];
            for (int vert=1; vert<nParameters+1; ++vert) {
                if(valuesAtVertices[vert]> secondWorstValue && vert != worst) {secondWorst = vert; secondWorstValue = valuesAtVertices[vert]; }
            }
            // -------

//************
//            std::cout << "ACCEPTING REFLECT \n";
//            viewStatusNelderMead<nParameters>(simplex, valuesAtVertices, worst, secondWorst, best, worstValue, secondWorstValue, bestValue);
//************
            
            continue;                                                                                                               // terminates iteration
        }
        
        
        // >>>> EXPAND:
        // >>>>>>>>>>>>
        if (reflectPointValue < bestValue) {
            for (int dim=0; dim < nParameters; ++dim){ expandPoint[dim] = centroid[dim]+gamma*(reflectPoint[dim] - centroid[dim]);} // Calculates expansion point
            expandPointValue = fun(expandPoint);                                                                                    // Calculates value at expansion point
            
//************
//            std::cout << "EXPAND "; viewPointNelderMead<nParameters>(expandPoint,expandPointValue);
//************
            if (greedyExpansion) {
                if (expandPointValue < bestValue) {
                simplex[worst] = expandPoint;                                                                                          // Accepts reflection point
                valuesAtVertices[worst] = expandPointValue;
                } else {
                    simplex[worst] = reflectPoint;                                                                                          // Accepts reflection point
                    valuesAtVertices[worst] = reflectPointValue;
                }
            } else {
            if (expandPointValue < reflectPointValue) {
                simplex[worst] = expandPoint;                                                                                          // Accepts reflection point
                valuesAtVertices[worst] = expandPointValue;
            } else {
                simplex[worst] = reflectPoint;                                                                                          // Accepts reflection point
                valuesAtVertices[worst] = reflectPointValue;
            }
            }
            // -------
            // Optimised Ordering
            // The new point is for sure the best
            // The secondWorstValue becomes worstValue & we need to find the new secondWorst
            // -------
            best = worst;           bestValue = valuesAtVertices[worst];
            worst = secondWorst;    worstValue = secondWorstValue;          // secondWorstValue becomes worstValue
            secondWorst = 0;        secondWorstValue = valuesAtVertices[0];
            for (int vert=1; vert<nParameters+1; ++vert) {
                if(valuesAtVertices[vert]> secondWorstValue && vert != worst) {secondWorst = vert; secondWorstValue = valuesAtVertices[vert]; }
            }
            
//************
//            std::cout << "ACCEPTING EXPAND or REFLECT \n";
//            viewStatusNelderMead<nParameters>(simplex, valuesAtVertices, worst, secondWorst, best, worstValue, secondWorstValue, bestValue);
//************
            
            continue;                                                                                                               // terminates iteration
        }
        
        // >>>> CONTRACT:
        // >>>>>>>>>>>>
        if (reflectPointValue >= secondWorstValue) {
            // >>>> CONTRACT OUTSIDE:
            if (reflectPointValue < worstValue) {
                for (int dim=0; dim < nParameters; ++dim){ contractPoint[dim] = centroid[dim]+beta*(reflectPoint[dim] - centroid[dim]);} // Calculates contraction point
                contractPointValue = fun(contractPoint);                                                                                 // Calculates value at contraction point
                
//************
//            std::cout << "CONTRACT OUTSIDE"; viewPointNelderMead<nParameters>(contractPoint,contractPointValue);
//************

                if (contractPointValue <= reflectPointValue) {
                    simplex[worst] = contractPoint;                                                                                      // Accepts contraction point
                    valuesAtVertices[worst] = contractPointValue;
                    // -------
                    // Ordering: in this case we have no information on the ordering
                    fullOrderingNelderMead<nParameters>(valuesAtVertices, worst, secondWorst, best, worstValue, secondWorstValue, bestValue);
                    // -------
                    
//************
//                    std::cout << "ACCEPTING CONTRACT OUTSIDE \n";
//                    viewStatusNelderMead<nParameters>(simplex, valuesAtVertices, worst, secondWorst, best, worstValue, secondWorstValue, bestValue);
//************
                    
                    continue;                                                                                                            // terminates iteration
                }
            } else {
            // >>>> CONTRACT INSIDE:
                for (int dim=0; dim < nParameters; ++dim){ contractPoint[dim] = centroid[dim]+beta*(simplex[worst][dim] - centroid[dim]);} // Calculates contraction point
                contractPointValue = fun(contractPoint);                                                                                   // Calculates value at contraction point
                
//************
//                std::cout << "CONTRACT INSIDE"; viewPointNelderMead<nParameters>(contractPoint,contractPointValue);
//************

                if (contractPointValue <= worstValue) {
                    simplex[worst] = contractPoint;                                                                                        // Accepts contraction point
                    valuesAtVertices[worst] = contractPointValue;
                    
                    // -------
                    // Ordering: in this case we have no information on the ordering
                    fullOrderingNelderMead<nParameters>(valuesAtVertices, worst, secondWorst, best, worstValue, secondWorstValue, bestValue);
                    // -------
                                        
//************
//                    std::cout << "ACCEPTING CONTRACT INSIDE \n";
//                    viewStatusNelderMead<nParameters>(simplex, valuesAtVertices, worst, secondWorst, best, worstValue, secondWorstValue, bestValue);
//************
                    continue;                                                                                                              // terminates iteration
                }
            }
        }
        
        // >>>> SHRINK:
        // >>>>>>>>>>>>
        for (int vert=0; vert<nParameters+1; ++vert){
            if (vert != best) {
                for (int dim=0; dim < nParameters; ++dim){
                    simplex[vert][dim] = simplex[best][dim]+delta*(simplex[vert][dim]-simplex[best][dim]);
                }
                valuesAtVertices[vert] = fun(simplex[vert]);
            }
        }
//************
//        std::cout << "SHRINK\n";
//************

        // -------
        // Ordering: in this case we have no information on the ordering
        fullOrderingNelderMead<nParameters>(valuesAtVertices, worst, secondWorst, best, worstValue, secondWorstValue, bestValue);
        // -------
        
//************
//        std::cout << "ACCEPTING SHRINK \n";
//        viewStatusNelderMead<nParameters>(simplex, valuesAtVertices, worst, secondWorst, best, worstValue, secondWorstValue, bestValue);
//************

    }
    
    
    // >>>> Returning results
    solution = simplex[best];
    for (int dim = 0; dim < nParameters; ++dim) {
        solutionVar[dim]=fabs(simplex[0][dim]-simplex[1][dim]);
        for (int vertex=2; vertex < nParameters+1; ++vertex) {
            if (fabs(simplex[0][dim]-simplex[vertex][dim]) > solutionVar[dim] ) {
                solutionVar[dim] = fabs(simplex[0][dim]-simplex[vertex][dim]);
            }
        }
    }
    valueAtSolution = bestValue;
    valueVarAtSolution = worstValue-bestValue;
       
    return convergenceStatus;
}


} // namespace Tortoise

#endif /* NelderMead_hpp */



