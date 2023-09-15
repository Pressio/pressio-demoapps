/*
//@HEADER
// ************************************************************************
//
// functor_fill_stencil.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef PRESSIODEMOAPPS_IMPL_GRADIENTS_2D_HPP_
#define PRESSIODEMOAPPS_IMPL_GRADIENTS_2D_HPP_

namespace pressiodemoapps{ namespace impl{

template<class ScalarType, class CellConnectivityVectorLike, class StateType>
ScalarType fd_face_normal_gradient_for_cell_centered_state_2d(const StateType & f,
							      const CellConnectivityVectorLike & cc,
							      char axis,
							      GradFdMode mode,
							      const ScalarType & h,
							      int nDofPerCell = 1,
							      int dofShift = 0)
{

  assert(axis == 'x' || axis == 'y');
  const int graphNumCols = cc.size();
  assert(graphNumCols ==5 || graphNumCols==9);

  const int i_05  = cc[0]; // this always exists
  const int i_p15 = (axis=='x') ? cc[3] : (axis=='y') ? cc[2] : -1;
  const int i_m15 = (axis=='x') ? cc[1] : (axis=='y') ? cc[4] : -1;
  const int i_p30 = (graphNumCols==9) ? ((axis=='x') ? cc[7] : ((axis=='y') ? cc[6] : -1)) : -1;
  const int i_m30 = (graphNumCols==9) ? ((axis=='x') ? cc[5] : ((axis=='y') ? cc[8] : -1)) : -1;

  const auto f_05  = f(i_05*nDofPerCell + dofShift);
  const auto f_p15 = (i_p15 != -1) ? f(i_p15*nDofPerCell + dofShift) : 0;
  const auto f_m15 = (i_m15 != -1) ? f(i_m15*nDofPerCell + dofShift) : 0;
  const auto f_p30 = (i_p30 != -1) ? f(i_p30*nDofPerCell + dofShift) : 0;
  const auto f_m30 = (i_m30 != -1) ? f(i_m30*nDofPerCell + dofShift) : 0;

  // https://web.media.mit.edu/~crtaylor/calculator.html
  switch(mode){
  case GradFdMode::ForwardTwoPt:  return (-f_05 + f_p15)/h;
  case GradFdMode::BackwardTwoPt: return ( f_05 - f_m15)/h;
  case GradFdMode::CenterThreePt: return (-f_m15 + f_p15)/(2*h);

  case GradFdMode::ForwardThreePt:  return (-2.*f_05 + 3.*f_p15 - 1.*f_p30)/h;
  case GradFdMode::BackwardThreePt: return ( 2.*f_05 - 3.*f_m15 + 1.*f_m30)/h;
  case GradFdMode::CenterFivePt:    return (-f_m30 - 8.*f_m15 + 8.*f_p15 - f_p30)/(12.*h);

  default:
    return 0;
  }
}

template<class ScT, std::size_t N>
struct Face{
  static constexpr std::size_t N_ = N;
  std::array<ScT, 3> centerCoords = {};
  std::array<ScT, N> normalGradValue = {};
};

struct MyKey{
  int parentCellGID;
  FacePosition pos;

  bool operator==(const MyKey & p) const {
    return parentCellGID == p.parentCellGID && pos == p.pos;
  }
};

struct hash_fn{
  std::size_t operator() (const MyKey & k) const{
    const std::size_t h1 = std::hash<int>()(k.parentCellGID);
    const std::size_t h2 = std::hash<int>()(static_cast<int>(k.pos));
    return h1 ^ h2;
  }
};

template<
  class MeshType,
  class FaceType = Face<typename MeshType::scalar_type, 1> >
class GradientInternal{
  using sc_t = typename MeshType::scalar_type;

public:
  GradientInternal(const MeshType & mesh,
		   BoundaryFacesGradientScheme schemeAtBoundaryFaces)
    : m_meshObj(mesh),
      m_schemeAtBoundaryFaces(schemeAtBoundaryFaces)
  {
    initialize(mesh);
  }

  template<class StateType>
  void compute(const StateType & state){
    if (m_schemeAtBoundaryFaces == BoundaryFacesGradientScheme::OneSidedFdAutoStencil){
      this->computeOneSidedFdAutoStencil(state);
    }
    else{
      throw std::runtime_error("invalid choice");
    }
  }

  const auto & queryFace(int cellGID, FacePosition fp) const{
    const MyKey key{cellGID, fp};
    std::cout << m_data.count(key) << '\n';
    assert(m_data.count(key) == 1);
    auto it = m_data.find(key);
    return it->second;
  }

private:
  template<class StateType>
  void computeOneSidedFdAutoStencil(const StateType & state)
  {
    constexpr std::size_t nDofPerCell = FaceType::N_;
    const int ss = m_meshObj.get().stencilSize();
    const auto SkewRight  = (ss == 3) ? GradFdMode::ForwardTwoPt
      : GradFdMode::ForwardThreePt;
    const auto SkewLeft   = (ss == 3) ? GradFdMode::BackwardTwoPt
      : GradFdMode::BackwardThreePt;

    const auto dx  = m_meshObj.get().dx();
    const auto dy  = m_meshObj.get().dy();

    const auto & rowsForLoop = m_meshObj.get().graphRowsOfCellsStrictlyOnBd();
    const auto & G = m_meshObj.get().graph();
    for (auto rowInd : rowsForLoop){
      const bool bdL = m_meshObj.get().cellHasLeftFaceOnBoundary2d(rowInd);
      const bool bdF = m_meshObj.get().cellHasFrontFaceOnBoundary2d(rowInd);
      const bool bdR = m_meshObj.get().cellHasRightFaceOnBoundary2d(rowInd);
      const bool bdB = m_meshObj.get().cellHasBackFaceOnBoundary2d(rowInd);
      //std::cout << "row = " << rowInd << " " << G(rowInd,0) << std::endl;

      const auto currCellConnec = G.row(rowInd);
      const int cellGID = currCellConnec(0);

      if (bdL){
	MyKey key{cellGID, FacePosition::Left};
	assert(m_data.count(key) == 1);
	auto & faceObj= m_data[key];
	for (int j=0; j<nDofPerCell; ++j){
	  faceObj.normalGradValue[j] =
	    fd_face_normal_gradient_for_cell_centered_state_2d(state, currCellConnec,
							       'x', SkewRight, dx,
							       nDofPerCell, j);
	}
      }

      // if (bdF){
      // 	MyKey key{cellGID, FacePosition::Front};
      // 	assert(m_data.count(key) == 1);
      // 	auto & faceObj= m_data[key];
      // 	for (int j=0; j<nDofPerCell; ++j){
      // 	  faceObj.normalGradValue[j] =
      // 	    fd_face_normal_gradient_for_cell_centered_state_2d(state, currCellConnec,
      // 							       'y', SkewLeft, dy,
      // 							       nDofPerCell, j);
      // 	}
      // }

      if (bdR){
	MyKey key{cellGID, FacePosition::Right};
	assert(m_data.count(key) == 1);
	auto & faceObj= m_data[key];
	for (int j=0; j<nDofPerCell; ++j){
	  faceObj.normalGradValue[j] =
	    fd_face_normal_gradient_for_cell_centered_state_2d(state, currCellConnec,
							       'x', SkewLeft, dx,
							       nDofPerCell, j);
	}
      }

      // if (bdB){
      // 	MyKey key{cellGID, FacePosition::Back};
      // 	assert(m_data.count(key) == 1);
      // 	auto & faceObj= m_data[key];
      // 	for (int j=0; j<nDofPerCell; ++j){
      // 	  faceObj.normalGradValue[j] =
      // 	    fd_face_normal_gradient_for_cell_centered_state_2d(state, currCellConnec,
      // 							       'y', SkewRight, dy,
      // 							       nDofPerCell, j);
      // 	}
      // }
    }
  }

private:
  void initialize(const MeshType & mesh){
    const auto & rowsForLoop = mesh.graphRowsOfCellsStrictlyOnBd();
    const auto & G = mesh.graph();
    const auto & x = mesh.viewX();
    const auto & y = mesh.viewY();
    const auto & z = mesh.viewZ();
    const auto half   = static_cast<sc_t>(1)/static_cast<sc_t>(2);
    const auto dxHalf = mesh.dx()*half;
    const auto dyHalf = mesh.dy()*half;

    for (auto rowInd : rowsForLoop){
      const bool bdL = mesh.cellHasLeftFaceOnBoundary2d(rowInd);
      const bool bdF = mesh.cellHasFrontFaceOnBoundary2d(rowInd);
      const bool bdR = mesh.cellHasRightFaceOnBoundary2d(rowInd);
      const bool bdB = mesh.cellHasBackFaceOnBoundary2d(rowInd);
      //std::cout << "row = " << rowInd << " " << G(rowInd,0) << std::endl;

      const auto graphRow = G.row(rowInd);
      const int cellGID   = graphRow(0);
      const auto cellX    = x(cellGID);
      const auto cellY    = y(cellGID);
      const auto cellZ    = z(cellGID);

      if (bdL) {
	m_data[MyKey{cellGID, FacePosition::Left}] = FaceType{{cellX-dxHalf, cellY, cellZ}};
      }

      // if (bdF) {
      // 	m_data[MyKey{cellGID, FacePosition::Front}] = FaceType{{cellX, cellY+dyHalf, cellZ}};
      // }

      if (bdR) {
	m_data[MyKey{cellGID, FacePosition::Right}] = FaceType{{cellX+dxHalf, cellY, cellZ}};
      }
      // if (bdB) {
      // 	m_data[MyKey{cellGID, FacePosition::Back}] = FaceType{{cellX, cellY-dyHalf, cellZ}};
      // }
    }
  }

private:
  BoundaryFacesGradientScheme m_schemeAtBoundaryFaces;
  std::reference_wrapper<const MeshType> m_meshObj;
  std::unordered_map<MyKey, FaceType, hash_fn> m_data = {};
};

}}
#endif
