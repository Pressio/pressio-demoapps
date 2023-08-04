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

template<class FuncLeftT, class FuncFrontT, class FuncRightT, class FuncBackT>
struct CustomBCsHolder;

template<class FuncLeftT, class FuncFrontT, class FuncRightT, class FuncBackT>
struct CustomBCsHolder
{
  CustomBCsHolder() = delete;

  template<class T1, class T2, class T3, class T4>
  CustomBCsHolder(T1 && a, T2 && b, T3 && c, T4 && d)
    : m_funcLeft(std::forward<T1>(a)),
      m_funcFront(std::forward<T2>(b)),
      m_funcRight(std::forward<T3>(c)),
      m_funcBack(std::forward<T4>(d)){}

  template<class ...Types>
  void operator() (GhostRelativeLocation rloc, Types && ... args) const noexcept
  {
    if (rloc == GhostRelativeLocation::Left){
      m_funcLeft(std::forward<Types>(args)...);
    }
    else if (rloc == GhostRelativeLocation::Front){
      m_funcFront(std::forward<Types>(args)...);
    }
    else if (rloc == GhostRelativeLocation::Right){
      m_funcRight(std::forward<Types>(args)...);
    }
    else{
      m_funcBack(std::forward<Types>(args)...);
    }
  }

private:
  FuncLeftT m_funcLeft;
  FuncFrontT m_funcFront;
  FuncRightT m_funcRight;
  FuncBackT m_funcBack;
};

template<class FuncLeftT, class FuncFrontT, class FuncRightT, class FuncBackT>
struct CustomBCsHolder<FuncLeftT &, FuncFrontT &, FuncRightT &, FuncBackT &>
{
  CustomBCsHolder() = delete;

  CustomBCsHolder(FuncLeftT & a, FuncFrontT & b,
		  FuncRightT & c, FuncBackT & d)
    : m_funcLeft(a), m_funcFront(b), m_funcRight(c), m_funcBack(d){}

  template<class ...Types>
  void operator() (GhostRelativeLocation rloc, Types && ... args) const noexcept
  {
    if (rloc == GhostRelativeLocation::Left){
      m_funcLeft(std::forward<Types>(args)...);
    }
    else if (rloc == GhostRelativeLocation::Front){
      m_funcFront(std::forward<Types>(args)...);
    }
    else if (rloc == GhostRelativeLocation::Right){
      m_funcRight(std::forward<Types>(args)...);
    }
    else{
      m_funcBack(std::forward<Types>(args)...);
    }
  }

private:
  std::reference_wrapper<const FuncLeftT> m_funcLeft;
  std::reference_wrapper<const FuncFrontT> m_funcFront;
  std::reference_wrapper<const FuncRightT> m_funcRight;
  std::reference_wrapper<const FuncBackT> m_funcBack;
};

}} //end namespace pressiodemoapps::impl
#endif
