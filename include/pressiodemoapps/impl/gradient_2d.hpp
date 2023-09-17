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

/*
  FRIZZI: Sept 15, 2023
  the code below has been done fairly quickly so not much thought has gone
  into how to best do things, but this is a starting point and note that all
  this is intentionally just impl details, the actual public API is kept separate
*/

template<class IntType, class ScalarType, class CellConnectivity, class FType>
ScalarType face_normal_gradient_for_cell_centered_function_2d(const FType & f,
							      const CellConnectivity & cc,
							      char axis,
							      GradFdMode mode,
							      const ScalarType & h,
							      int numDofPerCell,
							      int dofShift)
{

  assert(axis == 'x' || axis == 'y');
  const int graphNumCols = cc.size();
  assert(graphNumCols ==5 || graphNumCols==9 || graphNumCols==13);

  // fd rules taken from https://web.media.mit.edu/~crtaylor/calculator.html
  // manually verified only the 2pt one

  const IntType i_05  = cc[0]; // this always exists
  const IntType i_p15 = (axis=='x') ? cc[3] : (axis=='y') ? cc[2] : -1;
  const IntType i_m15 = (axis=='x') ? cc[1] : (axis=='y') ? cc[4] : -1;
  const IntType i_p30 = (graphNumCols>=9) ? ((axis=='x') ? cc[7] : ((axis=='y') ? cc[6] : -1)) : -1;
  const IntType i_m30 = (graphNumCols>=9) ? ((axis=='x') ? cc[5] : ((axis=='y') ? cc[8] : -1)) : -1;

  const auto f_05  = f(i_05*numDofPerCell + dofShift);
  const auto f_p15 = (i_p15 != -1) ? f(i_p15*numDofPerCell + dofShift) : 0;
  const auto f_m15 = (i_m15 != -1) ? f(i_m15*numDofPerCell + dofShift) : 0;
  const auto f_p30 = (i_p30 != -1) ? f(i_p30*numDofPerCell + dofShift) : 0;
  const auto f_m30 = (i_m30 != -1) ? f(i_m30*numDofPerCell + dofShift) : 0;

  switch(mode){
  case GradFdMode::ForwardTwoPt:  return (-f_05 + f_p15)/h;
  case GradFdMode::BackwardTwoPt: return ( f_05 - f_m15)/h;
  case GradFdMode::CenterThreePt: return (-f_m15 + f_p15)/(2*h);

  case GradFdMode::ForwardThreePt:  return (-2*f_05 + 3*f_p15 - 1*f_p30)/h;
  case GradFdMode::BackwardThreePt: return ( 2*f_05 - 3*f_m15 + 1*f_m30)/h;
  case GradFdMode::CenterFivePt:    return (-f_m30 - 8*f_m15 + 8*f_p15 - f_p30)/(12*h);

  default: return 0;
  }
}

constexpr int _faceNormalX = 1;
constexpr int _faceNormalY = 2;

template<class ScT, class IndexType, std::size_t N>
struct Face{
  std::array<ScT, 3> centerCoordinates = {};
  std::array<ScT, N> normalGradient = {};
  IndexType parentCellGraphRow = {};
  int normalDirection = {}; //1 for x, 2 for y

  static constexpr std::size_t N_ = N;
};

template<class ScT, class IndexType>
struct Face<ScT, IndexType, 1>{
  std::array<ScT, 3> centerCoordinates = {};
  ScT normalGradient = {};
  IndexType parentCellGraphRow = {};
  int normalDirection = {}; //1 for x, 2 for y

  static constexpr std::size_t N_ = 1;
};

template<class IndexType>
struct MyKey{
  IndexType parentCellGID;
  FacePosition pos;

  bool operator==(const MyKey & p) const {
    return parentCellGID == p.parentCellGID && pos == p.pos;
  }
};

struct HashFnc{
  template<class IndexType>
  std::size_t operator() (const MyKey<IndexType> & k) const{
    const std::size_t h1 = std::hash<IndexType>()(k.parentCellGID);
    const std::size_t h2 = std::hash<int>()(static_cast<int>(k.pos));
    return h1 ^ h2;
  }
};


template<class MeshType, class FaceType>
class GradientEvaluatorInternal
{
  using sc_t = typename MeshType::scalar_type;
  using index_t	= typename MeshType::index_t;
  using key_t = MyKey<index_t>;

public:
  // constructor for when we only want the normal grad at domain boundaries' faces
  GradientEvaluatorInternal(const MeshType & mesh,
			    BoundaryFacesNormalGradientScheme schemeAtBoundaryFaces)
    : m_schemeAtBoundaryFaces(schemeAtBoundaryFaces),
      m_meshObj(mesh)
  {
    initializeForStoringNormalGradsAtBoundaryFaces(mesh);
  }

  const auto & queryFace(index_t cellGID, FacePosition fp) const{
    const key_t key{cellGID, fp};
    assert(m_data.count(key) == 1);
    auto it = m_data.find(key);
    return it->second;
  }

  template<class FieldType>
  void compute(const FieldType & field, int numDofPerCell)
  {
    if (m_schemeAtBoundaryFaces == BoundaryFacesNormalGradientScheme::OneSidedFdAutoStencil){
      this->normalGradBoundaryFacesOneSidedFdAutoStencil(field, numDofPerCell);
    }
    else{
      throw std::runtime_error("invalid choice for BoundaryFacesNormalGradientScheme");
    }
  }

private:
  template<class FieldType>
  void normalGradBoundaryFacesOneSidedFdAutoStencil(const FieldType & field,
						    int numDofPerCell)
  {
    /*
      compute normal gradients using ond-sided FD where the width of
      stencil used is determined based on stencil in the mesh object
      if the mesh stencil size == 3, we use a two point one-sided FD
      if the mesh stencil size == 5, we use a three point one-sided FD
      if the mesh stencil size == 7, we use a three point one-sided FD
        - we could make this bigger later
     */

    // determine the FD stencil using mesh stencil size
    const int ss = m_meshObj.get().stencilSize();
    const auto RightSided  = (ss == 3) ? GradFdMode::ForwardTwoPt
                                        : GradFdMode::ForwardThreePt;
    const auto LeftSided   = (ss == 3) ? GradFdMode::BackwardTwoPt
				        : GradFdMode::BackwardThreePt;

    const auto dx  = m_meshObj.get().dx();
    const auto dy  = m_meshObj.get().dy();
    const auto & G = m_meshObj.get().graph();

    for (auto it=m_data.begin(); it!=m_data.end(); ++it){
      const auto fpos = it->first.pos;
      auto & face = it->second;
      const int parentCellRowInd = face.parentCellGraphRow;
      const auto currCellConnec  = G.row(parentCellRowInd);

      GradFdMode fdMode =
	(fpos == FacePosition::Left) ? RightSided
	   : (fpos == FacePosition::Back) ? RightSided
	        : (fpos == FacePosition::Right) ? LeftSided
	             : LeftSided;

      const char axis   = (face.normalDirection == _faceNormalX) ? 'x' : 'y';
      const auto h      = (face.normalDirection == _faceNormalX) ? dx : dy;

      if constexpr (FaceType::N_ == 1){
	face.normalGradient =
	  face_normal_gradient_for_cell_centered_function_2d<index_t>(field, currCellConnec, axis,
							     fdMode, h, numDofPerCell, 0);
      }
      else{
	for (int j=0; j<numDofPerCell; ++j){
	  face.normalGradient[j] =
	    face_normal_gradient_for_cell_centered_function_2d<index_t>(field, currCellConnec, axis,
									fdMode, h, numDofPerCell, j);
	}
      }
    }
  }

private:
  void initializeForStoringNormalGradsAtBoundaryFaces(const MeshType & mesh){
    const auto & G = mesh.graph();
    const auto & x = mesh.viewX();
    const auto & y = mesh.viewY();
    const auto & z = mesh.viewZ();
    constexpr auto half   = static_cast<sc_t>(1)/static_cast<sc_t>(2);
    const auto dxHalf = mesh.dx()*half;
    const auto dyHalf = mesh.dy()*half;

    for (auto rowInd : mesh.graphRowsOfCellsStrictlyOnBd())
    {
      const bool bL = mesh.cellHasLeftFaceOnBoundary2d(rowInd);
      const bool bF = mesh.cellHasFrontFaceOnBoundary2d(rowInd);
      const bool bR = mesh.cellHasRightFaceOnBoundary2d(rowInd);
      const bool bB = mesh.cellHasBackFaceOnBoundary2d(rowInd);

      const int cellGID = G.row(rowInd)[0];
      const auto cellX  = x(cellGID);
      const auto cellY  = y(cellGID);
      const auto cellZ  = z(cellGID);
      const auto faceCZ = cellZ;

      if (bL){
	const auto faceCX = cellX-dxHalf;
	const auto faceCY = cellY;
	m_data[key_t{cellGID, FacePosition::Left}] =
	  FaceType{{faceCX, faceCY, faceCZ}, {}, rowInd, _faceNormalX};
      }
      if (bF){
	const auto faceCX = cellX;
	const auto faceCY = cellY+dyHalf;
	m_data[key_t{cellGID, FacePosition::Front}] =
	  FaceType{{faceCX, faceCY, faceCZ}, {}, rowInd, _faceNormalY};
      }
      if (bR){
	const auto faceCX = cellX+dxHalf;
	const auto faceCY = cellY;
	m_data[key_t{cellGID, FacePosition::Right}] =
	  FaceType{{faceCX, faceCY, faceCZ}, {}, rowInd, _faceNormalX};
      }
      if (bB){
	const auto faceCX = cellX;
	const auto faceCY = cellY-dyHalf;
	m_data[key_t{cellGID, FacePosition::Back}] =
	  FaceType{{faceCX, faceCY, faceCZ}, {}, rowInd, _faceNormalY};
      }
    }
  }

private:
  BoundaryFacesNormalGradientScheme m_schemeAtBoundaryFaces;
  std::reference_wrapper<const MeshType> m_meshObj;
  std::unordered_map<key_t, FaceType, HashFnc> m_data = {};
};

}}
#endif
