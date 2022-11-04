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
//  DataMatrix.hpp
//  TORTOISE
//
//  Created by Marco Battiato on 11/3/22.
//

#ifndef DataMatrix_hpp
#define DataMatrix_hpp

#include <Geometry/Structured/FunctionSpace/Function.hpp>
#include <Generics/Containers/StringMap.hpp>


#include <stdio.h>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>      // std::setw

namespace Tortoise {

class DataColumn: public std::vector<double>{
public:
    using std::vector<double>::operator=;
    template <int NDim> DataColumn& operator=(const Function<NDim>& other){
        auto mat = other.toNodalFromLinForm();
        this->std::vector<double>::operator=(std::vector<double>(mat.data(), mat.data() + mat.rows() * mat.cols()));
        return *this;
    };
};

class DataMatrix: public StringMap<DataColumn>{
public:
    DataMatrix(const std::vector<std::string>& colNames) {
        for (auto col: colNames){ emplace_back(col, DataColumn()); }
    }

    DataMatrix& addColumn(const std::string& colName){
        emplace_back(colName, DataColumn());
        return *this;
    }
    
    // Sort will crash of the columns do not have the same lengths
    void sortBy(const std::string& colName){
        int idOfMain = findID(colName);
        assert(idOfMain>=0);
        std::vector<std::vector<double>> data;
        data.resize(operator[](colName).size());
        for (int i=0; i<operator[](colName).size(); ++i){
            data[i].resize(size());
            data[i][0] = operator[](colName)[i];
            int cur = 1;
            for (int j=0; j<size(); ++j){
                if (j!=idOfMain) { data[i][cur] = operator[](j)[i]; ++cur; }
            }
        }
        std::sort(data.begin(), data.end(), [](const std::vector<double>& a, const std::vector<double>& b) { return a[0] < b[0]; });
        for (int i=0; i<operator[](colName).size(); ++i){
            operator[](colName)[i] = data[i][0];
            int cur = 1;
            for (int j=0; j<size(); ++j){
                if (j!=idOfMain) { operator[](j)[i] = data[i][cur]; ++cur; }
            }
        }
    }
    
    void printToFile(const std::string& fileName){
        std::ofstream myfile;
        myfile.open (fileName);
        int  currOutputLines = 0;
        bool done = true;    // false if still lines to print to file
        myfile << "# ";
        for (int i=0; i<this->size(); ++i){
            std::string title(this->id(i)); if (title.size()>10) { title.resize(10);}
            myfile << std::setw(10) << title << " ";
            if ((*this)[i].size()>currOutputLines) { done = false; }
        }
        myfile << "\n";
        
        while (!done) {
            done = true;
            for (int i=0; i<this->size(); ++i){
                if ((*this)[i].size()>currOutputLines) {
                    myfile << std::setw(10) << (*this)[i][currOutputLines] << " ";
                } else {
                    myfile << std::setw(10) << " " << " ";
                }
                if ((*this)[i].size()>currOutputLines+1) { done = false; }
            }
            myfile << "\n";
            ++currOutputLines;
        }
        myfile.close();
    }
};

} // namespace Tortoise


#endif /* DataMatrix_hpp */
