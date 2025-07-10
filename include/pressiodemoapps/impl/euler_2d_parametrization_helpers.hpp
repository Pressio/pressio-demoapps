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

#ifndef EULER_2D_PARAMETRIZATION_HELPERS_HPP_
#define EULER_2D_PARAMETRIZATION_HELPERS_HPP_

namespace pressiodemoapps{
namespace impleuler2d{

constexpr int invalidIndex = std::numeric_limits<int>::max();

// NOTE: physical and IC params must be associated with
// contiguous but separate indexes because they are stored in separate vectors
constexpr int gamma_i  = 0;

constexpr int icNormalShock_mach_i          = 0;
constexpr int icCrossShock_density_i        = 1;
constexpr int icCrossShock_inletXVel_i      = 2;
constexpr int icCrossShock_bottomYVel_i     = 3;
constexpr int icRiemann1_topRightPressure_i = 4;
constexpr int icRiemann2_topRightPressure_i = 5;
constexpr int icRiemann2_topRightXVel_i     = 6;
constexpr int icRiemann2_topRightYVel_i     = 7;
constexpr int icRiemann2_topRightDensity_i  = 8;
constexpr int icRiemann2_botLeftPressure_i  = 9;

// IMPORTANT: we need to fill the vectors with default vaues
// with ALL default paramters **consistently** with the indexing defined above
template<class ScalarType>
const std::vector<ScalarType> defaultPhysicalParams
({
  static_cast<ScalarType>(1.4) //gamma
});

template<class ScalarType>
const std::vector<ScalarType> defaultInitCondParams
({
    static_cast<ScalarType>(9)    //normalShockMach
  , static_cast<ScalarType>(0.1)  //cross shock density
  , static_cast<ScalarType>(10.)  //cross shock inlet X valoc
  , static_cast<ScalarType>(1.)   //cross shock bottom Y veloc
  , static_cast<ScalarType>(0.4)  //this is for riemann IC 1
  , static_cast<ScalarType>(1.5)  //riemann IC 2 top right pressure
  , static_cast<ScalarType>(0.0)  //riemann IC 2 top right x-velocity
  , static_cast<ScalarType>(0.0)  //riemann IC 2 top right y-velocity
  , static_cast<ScalarType>(1.5)  //riemann IC 2 top right density
  , static_cast<ScalarType>(0.029)  //riemann IC 2 bottom left pressure
});

const std::vector<std::string> paramNames({
    "gamma",
    "normalShockMach",
    "crossShockDensity", "crossShockInletXVel", "crossShockBottomYVel",
    "riemannTopRightPressure", "riemannTopRightXVel", "riemannTopRightYVel",
    "riemannTopRightDensity", "riemannBotLeftPressure",
  });

template<class T = void>
bool is_physical(const std::string & s){
  return (s == "gamma");
}

template<class T = void>
int phys_param_string_to_index(const std::string & s){
  if      (s == "gamma")  { return gamma_i; }
  else                      { return invalidIndex; }
}

template<class T = void>
int ic_param_string_to_index(const Euler2d probEn, const int icFlag, const std::string & s)
{
  switch(probEn){
  case Euler2d::NormalShock:
    switch(icFlag){
    case 1:
      if      (s == "mach") { return icNormalShock_mach_i; }
      else                  { return invalidIndex; }
    default: return invalidIndex;
    }

  case Euler2d::CrossShock:
    switch(icFlag){
    case 1:
      if      (s == "crossShockDensity")    { return icCrossShock_density_i; }
      else if (s == "crossShockInletXVel")  { return icCrossShock_inletXVel_i; }
      else if (s == "crossShockBottomYVel") { return icCrossShock_bottomYVel_i; }
      else                                  { return invalidIndex; }
    default: return invalidIndex;
    }

  case Euler2d::Riemann:
    switch(icFlag){
    case 1:
      if      (s == "riemannTopRightPressure") { return icRiemann1_topRightPressure_i; }
      else                                     { return invalidIndex; }
    case 2:
      if      (s == "riemannTopRightPressure") { return icRiemann2_topRightPressure_i; }
      else if (s == "riemannTopRightXVel")     { return icRiemann2_topRightXVel_i; }
      else if (s == "riemannTopRightYVel")     { return icRiemann2_topRightYVel_i; }
      else if (s == "riemannTopRightDensity")  { return icRiemann2_topRightDensity_i; }
      else if (s == "riemannBotLeftPressure")  { return icRiemann2_botLeftPressure_i; }
      else                                     { return invalidIndex; }
    default: return invalidIndex;
    }

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
					const Euler2d probEn,
					const int icFlag,
					const std::unordered_map<std::string, ScalarType> & map)
{

  auto action = [&](auto it){
    if (is_physical(it.first)){ physVec[phys_param_string_to_index<>(it.first)] = it.second; }
    else		      { icVec[ic_param_string_to_index<>(probEn, icFlag, it.first)] = it.second;}
  };
  std::for_each(map.cbegin(), map.cend(), action);
}

}}//end namespace
#endif
