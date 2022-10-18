/*
//@HEADER
// ************************************************************************
//
// resize.hpp
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

#ifndef PRESSIODEMOAPPS_RESIZE_FUNC_HPP_
#define PRESSIODEMOAPPS_RESIZE_FUNC_HPP_

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace pressiodemoapps{

template<class value_type, class sizetype>
void resize(std::vector<value_type> & a,
	    sizetype newSize)
{
  a.resize(newSize);
}

template<class value_type, class sizetype>
void resize(Eigen::Matrix<value_type, -1, 1> & a,
	    sizetype newSize)
{
  a.resize(newSize);
}

template<class value_type, class sizetype>
void resize(Eigen::Matrix<value_type, -1, -1, Eigen::RowMajor> & M,
	    sizetype rows,
	    sizetype cols)
{
  M.resize(rows, cols);
}

template<class value_type, class sizetype>
void resize(Eigen::Matrix<value_type, -1, -1, Eigen::ColMajor> & M,
	    sizetype rows,
	    sizetype cols)
{
  M.resize(rows, cols);
}

template<class value_type, class IntType, class sizetype>
void resize(Eigen::SparseMatrix<value_type, Eigen::RowMajor, IntType> & M,
	    sizetype rows,
	    sizetype cols)
{
  M.resize(rows, cols);
}

template<class value_type, class IntType, class sizetype>
void resize(Eigen::SparseMatrix<value_type, Eigen::ColMajor, IntType> & M,
	    sizetype rows,
	    sizetype cols)
{
  M.resize(rows, cols);
}

#ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
template<class T, class sizetype>
typename std::enable_if<
  pressiodemoapps::predicates::is_array_pybind<T>::value
  >::type
resize(T & o, sizetype newSize)
{
  if (o.ndim() != 1 ){
    throw std::runtime_error
      ("Resize overload for ndim=1 called on pybind array with ndim!=1");
  }

  o.resize({newSize}, false);
}

template<class T, class sizetype>
typename std::enable_if<
  pressiodemoapps::predicates::is_array_pybind<T>::value
  >::type
resize(T & o, sizetype rows, sizetype cols)
{
  if (o.ndim() != 2 ){
    throw std::runtime_error
      ("Resize overload for ndim=2 called on pybind array with ndim!=2");
  }
  o.resize({rows, cols}, false);
}
#endif

}
#endif
