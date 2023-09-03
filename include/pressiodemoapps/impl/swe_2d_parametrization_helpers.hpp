/*
//@HEADER
// ************************************************************************
//
// swe_2d_prob_class.hpp
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

#ifndef SWE_2D_PARAMETRIZATION_HELPERS_HPP_
#define SWE_2D_PARAMETRIZATION_HELPERS_HPP_

namespace pressiodemoapps{
namespace implswe2d{

constexpr int invalidIndex = std::numeric_limits<int>::max();

// NOTE: physical and IC params must be associated with
// contiguous but separate indexes because they are stored in separate vectors
constexpr int gravity_i  = 0;
constexpr int coriolis_i = 1;

// for ic1, ic2
constexpr int ic1_pulseMag_i  = 0;
constexpr int ic1_pulseX_i    = 1;
constexpr int ic1_pulseY_i    = 2;
constexpr int ic2_pulseMag1_i = 3;
constexpr int ic2_pulseX1_i   = 4;
constexpr int ic2_pulseY1_i   = 5;
constexpr int ic2_pulseMag2_i = 6;
constexpr int ic2_pulseX2_i   = 7;
constexpr int ic2_pulseY2_i   = 8;

// IMPORTANT: we need to fill the vectors with default vaues
// with ALL default paramters **consistently** with the indexing defined above
template<class ScalarType>
const std::vector<ScalarType> defaultPhysicalParams
({
  static_cast<ScalarType>(9.8), //gravity
  static_cast<ScalarType>(-3), //coriolis
});

template<class ScalarType>
const std::vector<ScalarType> defaultInitCondParams
({
  static_cast<ScalarType>(1)/8, //pulseMag
  static_cast<ScalarType>(1),   //pulseX
  static_cast<ScalarType>(1),   //pulseY
  // this is the two gaussians pulses
  static_cast<ScalarType>(1)/10,//pulseMag1
  static_cast<ScalarType>(-2),  //pulseX1
  static_cast<ScalarType>(-2),  //pulseY1
  static_cast<ScalarType>(1)/8, //pulseMag2
  static_cast<ScalarType>(2),   //pulseX2
  static_cast<ScalarType>(2)    //pulseY2
});

const std::vector<std::string> paramNames({
    "gravity", "coriolis",
    "pulseMagnitude", "pulseX", "pulseY",
    "pulseMagnitude1", "pulseX1", "pulseY1",
    "pulseMagnitude2", "pulseX2", "pulseY2"
  });

template<class T = void>
bool is_physical(const std::string & s){
  return (s == "gravity" || s == "coriolis");
}

template<class T = void>
int phys_param_string_to_index(const std::string & s){
  if      (s == "gravity")  { return gravity_i; }
  else if (s == "coriolis") { return coriolis_i; }
  else                      { return invalidIndex; }
}

template<class T = void>
int ic_param_string_to_index(const int icFlag, const std::string & s){
  switch(icFlag){
    case 1:
      if      (s == "pulseMagnitude") { return ic1_pulseMag_i; }
      else if (s == "pulseX")	      { return ic1_pulseX_i; }
      else if (s == "pulseY")	      { return ic1_pulseY_i; }
      else                            { return invalidIndex; }

    case 2:
      if      (s == "pulseMagnitude1")  { return ic2_pulseMag1_i; }
      else if (s == "pulseX1")		{ return ic2_pulseX1_i; }
      else if (s == "pulseY1")		{ return ic2_pulseY1_i; }
      else if (s == "pulseMagnitude2")  { return ic2_pulseMag2_i; }
      else if (s == "pulseX2")		{ return ic2_pulseX2_i; }
      else if (s == "pulseY2")		{ return ic2_pulseY2_i; }
      else				{ return invalidIndex; }

    default: return invalidIndex;
  }
}

template<class T = void>
bool valid_parameter_name(const std::string & pname){
  auto match = [&](const std::string & testS){ return pname == testS; };
  return std::any_of(paramNames.cbegin(), paramNames.cend(), match);
}

template<class ScalarType>
bool contains_valid_parameter_names(const std::unordered_map<std::string, ScalarType> & map){
  return std::all_of(map.cbegin(), map.cend(), [](const auto & it){ return valid_parameter_name(it.first); });
}

template<class ScalarType>
void replace_params_from_map_if_present(std::vector<ScalarType> & physVec,
					std::vector<ScalarType> & icVec,
					const int icFlag,
					const std::unordered_map<std::string, ScalarType> & map)
{

  auto action = [&](auto it){
    if (is_physical(it.first)){ physVec[phys_param_string_to_index<>(it.first)] = it.second; }
    else		      { icVec[ic_param_string_to_index<>(icFlag, it.first)] = it.second;}
  };
  std::for_each(map.cbegin(), map.cend(), action);
}

}}//end namespace
#endif
