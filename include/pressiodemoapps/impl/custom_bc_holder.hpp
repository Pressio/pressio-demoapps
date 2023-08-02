/*
//@HEADER
// ************************************************************************
//
// euler2d.hpp
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

#ifndef PRESSIODEMOAPPS_IMPL_CUSTOM_BCS_HOLDER_HPP_
#define PRESSIODEMOAPPS_IMPL_CUSTOM_BCS_HOLDER_HPP_

#include "./ghost_relative_locations.hpp"

namespace pressiodemoapps{
namespace impl{

template<class T1, class T2, class T3, class T4>
struct CustomBCsHolder;

template<class T1, class T2, class T3, class T4>
struct CustomBCsHolder
{
  CustomBCsHolder() = delete;

  template<class T1l, class T2l, class T3l, class T4l>
  CustomBCsHolder(T1l && t1l, T2l && t2l, T3l && t3l, T4l && t4l)
    : a_(std::forward<T1l>(t1l)),
      b_(std::forward<T2l>(t2l)),
      c_(std::forward<T3l>(t3l)),
      d_(std::forward<T4l>(t4l)){}

  template<class ...Types>
  void operator() (GhostRelativeLocation rloc, Types && ... args) const noexcept
  {
    if (rloc == GhostRelativeLocation::Left){
      a_(std::forward<Types>(args)...);
    }
    else if (rloc == GhostRelativeLocation::Front){
      b_(std::forward<Types>(args)...);
    }
    else if (rloc == GhostRelativeLocation::Right){
      c_(std::forward<Types>(args)...);
    }
    else{
      d_(std::forward<Types>(args)...);
    }
  }

  T1 a_;
  T2 b_;
  T3 c_;
  T4 d_;
};

template<class T1, class T2, class T3, class T4>
struct CustomBCsHolder<T1 &, T2 &, T3 &, T4 &>
{
  CustomBCsHolder() = delete;

  CustomBCsHolder(T1 & t1l, T2 & t2l, T3 & t3l, T4 & t4l)
    : a_(t1l), b_(t2l), c_(t3l), d_(t4l){}

  template<class ...Types>
  void operator() (GhostRelativeLocation rloc, Types && ... args) const noexcept
  {
    if (rloc == GhostRelativeLocation::Left){
      a_(std::forward<Types>(args)...);
    }
    else if (rloc == GhostRelativeLocation::Front){
      b_(std::forward<Types>(args)...);
    }
    else if (rloc == GhostRelativeLocation::Right){
      c_(std::forward<Types>(args)...);
    }
    else{
      d_(std::forward<Types>(args)...);
    }
  }

  std::reference_wrapper<const T1> a_;
  std::reference_wrapper<const T2> b_;
  std::reference_wrapper<const T3> c_;
  std::reference_wrapper<const T4> d_;
};

}} //end namespace pressiodemoapps::impl
#endif
