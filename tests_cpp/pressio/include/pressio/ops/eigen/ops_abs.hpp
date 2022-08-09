/*
//@HEADER
// ************************************************************************
//
// ops_abs.hpp
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

#ifndef OPS_EIGEN_OPS_ABS_HPP_
#define OPS_EIGEN_OPS_ABS_HPP_

namespace pressio{ namespace ops{

template <class T1, class T2>
::pressio::mpl::enable_if_t<
  /* a bit restrictive but fine for now */
  ::pressio::all_have_traits_and_same_scalar<T1, T2>::value
  && ::pressio::Traits<T1>::package_identifier == PackageIdentifier::Eigen
  && ::pressio::Traits<T2>::package_identifier == PackageIdentifier::Eigen
  && ::pressio::Traits<T1>::rank == 1
  && ::pressio::Traits<T2>::rank == 1
  >
abs(T1 & y, const T2 & x)
{

  using ord_t = decltype( ::pressio::ops::extent(x, 0) );
  for (ord_t i=0; i< ::pressio::ops::extent(x, 0); ++i){
    y(i) = std::abs(x(i));
  }
}

}}//end namespace pressio::ops
#endif  // OPS_EIGEN_OPS_ABS_HPP_