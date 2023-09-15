/*
//@HEADER
// ************************************************************************
//
// mesh.hpp
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

#ifndef PRESSIODEMOAPPS_GRADIENT_PUBLIC_HPP_
#define PRESSIODEMOAPPS_GRADIENT_PUBLIC_HPP_

#include "mesh.hpp"
#include "schemes_info.hpp"
#include "impl/gradient_2d.hpp"

namespace pressiodemoapps{

constexpr BoundaryFacesNormalGradientScheme defaultGradSchemeForBoundaryFaces =
  BoundaryFacesNormalGradientScheme::OneSidedFdAutoStencil;

template<class MeshType, std::size_t MaxNumDofsPerCell = 1>
class GradientEvaluator
{
public:
  explicit GradientEvaluator(const MeshType & mesh)
    : m_impl(mesh, defaultGradSchemeForBoundaryFaces)
  {
    if (mesh.dimensionality() != 2){
      throw std::runtime_error("gradients currently only supported for 2D");
    }
  }

  // this overload is for when the field has one dof/value per cell
  template<class FieldType>
  void operator()(const FieldType & field){
    constexpr int numDofsPerCell = 1;
    m_impl.compute(field, numDofsPerCell);
  }

  // this overload is for when the field has more than one dof/value per cell
  template<class FieldType>
  void operator()(const FieldType & field, int numDofsPerCell){
    assert(numDofsPerCell > 1);
    m_impl.compute(field, numDofsPerCell);
  }

  /* return a reference to a struct that meets the following API.

     - if MaxNumDofsPerCell == 1

	struct FaceApiExpositionOnly{
	  std::array<scalar_type, 3> centerCoordinates;
	  scalar_type normalGradient;
	  int normalDirection; // == 1 for normal along x, == 2 for normal along y
	};

     - if MaxNumDofsPerCell >= 2

	struct FaceApiExpositionOnly{
	  std::array<scalar_type, 3> centerCoordinates;
	  std::array<scalar_type, MaxNumDofsPerCell> normalGradient;
	  int normalDirection; // == 1 for normal along x, == 2 for normal along y
	};

     NOTE: this is just for exposition-only to show what it offers,
     it is not the actual name of the class.
     You do not need to know the actual type, just use auto.
  */
  const auto & queryFace(int cellGID, FacePosition fp) const{
    return m_impl.queryFace(cellGID, fp);
  }

private:
  using face_t = impl::Face<typename MeshType::scalar_type, MaxNumDofsPerCell>;
  impl::GradientEvaluatorInternal<MeshType, face_t> m_impl;
};

}//end namespace pressiodemoapps
#endif
