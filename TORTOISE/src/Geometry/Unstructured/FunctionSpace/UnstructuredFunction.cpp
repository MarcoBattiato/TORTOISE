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
//  UnstructuredFunction.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 7/8/21.
//
//  INCOMPLETE

#include <Geometry/Unstructured/FunctionSpace/UnstructuredFunction.hpp>

namespace Tortoise {

UnstructuredFunction1D::UnstructuredFunction1D(const UnstructuredMesh1D& t_mesh): mesh(&t_mesh), vec(DataVector::Zero(t_mesh.dataVectorDim)){};
//UnstructuredFunction1D::UnstructuredFunction1D(const UnstructuredMesh1D& mesh, const real value): mesh(&t_mesh), vec(DataVector::Zero(t_mesh.dataVectorDim)){
//    Eigen::Map<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), 2, mesh->numberElements);
//    elementview.row(0).array() = t_value;
//};
//UnstructuredFunction1D(const UnstructuredMesh1D& mesh, const DataVector& nodalValues): mesh(&t_mesh), vec(t_nodalValues){
//    assert(vec.size() == t_mesh.dataVectorDim);
//    Eigen::Map<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), 2, mesh->numberElements);
//    elementview.block(1,0,mesh->nVertElement-1,mesh->numberElements).rowwise() -= elementview.row(0);
//};
//UnstructuredFunction1D(const UnstructuredMesh1D& mesh, const std::function<real(real)>& f): mesh(&t_mesh), vec(t_mesh.dataVectorDim){
//    for (int elem = 0; elem < mesh->numberElements ; ++elem) {
//        for (int vert = 0; vert < 2; ++vert){
//            vec(vert + mesh->nVertElement * elem) = t_f(mesh->operator()(elem,vert));
//        }
//    }
//    Eigen::Map<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>> elementview(vec.data(), 2, mesh->numberElements);
//    elementview.block(1,0,1,mesh->numberElements).rowwise() -= elementview.row(0);
//}
//
//
//template<int NDim> Function<NDim>::Function(const Mesh<NDim>& t_mesh, const std::function<real(Point<NDim>)>& t_f)
//
//
//explicit UnstructuredFunction1D(const UnstructuredMesh1D* mesh);                                                  // Constructs a constant 0 function
//UnstructuredFunction1D(const UnstructuredMesh1D* mesh, const real value);                                // Constructs a constant function
//UnstructuredFunction1D(const UnstructuredMesh1D* mesh, const DataVector& nodalValues);                            // Constructs using the values at nodes
//UnstructuredFunction1D(const UnstructuredMesh1D* mesh, const std::function<real(real)>& f);       // Constructs using the function passed
//UnstructuredFunction1D(const std::function<real(real)>& f, const UnstructuredFunction1D& other);// Constructs a function of a function f(other)
//UnstructuredFunction1D(const std::function<real(std::vector<real>)>& f, const std::vector<UnstructuredFunction1D>& t_others);// Constructs a function of a function f(other[0], other[1], ...)
//UnstructuredFunction1D(const UnstructuredFunction1D& t_other);                                                        // Copy constructor
//UnstructuredFunction1D(UnstructuredFunction1D&& t_other);                                                             // Move constructor


} // namespace Tortoise
