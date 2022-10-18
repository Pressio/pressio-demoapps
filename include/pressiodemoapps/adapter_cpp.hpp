/*
//@HEADER
// ************************************************************************
//
// adapter_cpp.hpp
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

#ifndef PRESSIODEMOAPPS_ADAPTER_MIXINS_CPP_HPP_
#define PRESSIODEMOAPPS_ADAPTER_MIXINS_CPP_HPP_

namespace pressiodemoapps{

// this is the class that gets instantiated and passed
// to the user when they create any problem.
// The specific problem class they want becomes the template T.
template<class T>
class PublicProblemEigenMixinCpp : public T
{

public:
  using scalar_type               = typename T::scalar_type;
  using independent_variable_type = typename T::scalar_type;
  using time_type                 = typename T::scalar_type;
  using state_type		  = typename T::state_type;
  using right_hand_side_type      = typename T::velocity_type;
  using jacobian_type             = typename T::jacobian_type;

private:
  using typename T::index_t;
  using eigen_vec_t    = Eigen::Matrix<typename T::scalar_type, Eigen::Dynamic, 1>;
  using eigen_mat_ll_t = Eigen::Matrix<typename T::scalar_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  using eigen_mat_lr_t = Eigen::Matrix<typename T::scalar_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

public:
  template<class ...Args>
  PublicProblemEigenMixinCpp(Args && ... args)
    : T(std::forward<Args>(args)...)
  {
    T::initializeJacobian(m_jacobian);
    ::pressiodemoapps::set_zero(m_jacobian);
    ::pressiodemoapps::resize(m_rhs, T::m_numDofSampleMesh);
    ::pressiodemoapps::set_zero(m_rhs);
  }

  typename T::index_t totalDofSampleMesh()  const{
    return T::m_numDofSampleMesh;
  }

  typename T::index_t totalDofStencilMesh() const{
    return T::m_numDofStencilMesh;
  }

  state_type createState() const{
    state_type tmp = T::initialCondition();
    ::pressiodemoapps::set_zero(tmp);
    return tmp;
  }

  right_hand_side_type createRightHandSide() const{
    namespace pda = ::pressiodemoapps;
    return pda::clone(m_rhs);
  }

  jacobian_type createJacobian() const{
    namespace pda = ::pressiodemoapps;
    return pda::clone(m_jacobian);
  }

  //
  // jacobian action
  //
  eigen_vec_t createApplyJacobianResult(const eigen_vec_t & operand) const{
    eigen_vec_t res(m_jacobian.rows());
    ::pressiodemoapps::set_zero(res);
    return res;
  }

  eigen_mat_ll_t createApplyJacobianResult(const eigen_mat_ll_t & operand) const{
    eigen_mat_ll_t res(m_jacobian.rows(), operand.cols());
    ::pressiodemoapps::set_zero(res);
    return res;
  }

  eigen_mat_lr_t createApplyJacobianResult(const eigen_mat_lr_t & operand) const{
    eigen_mat_lr_t res(m_jacobian.rows(), operand.cols());
    ::pressiodemoapps::set_zero(res);
    return res;
  }

  eigen_vec_t createResultOfJacobianActionOn(const eigen_vec_t & operand) const{
    return createApplyJacobianResult(operand);
  }

  eigen_mat_ll_t createResultOfJacobianActionOn(const eigen_mat_ll_t & operand) const{
    return createApplyJacobianResult(operand);
  }

  eigen_mat_lr_t createResultOfJacobianActionOn(const eigen_mat_lr_t & operand) const{
    return createApplyJacobianResult(operand);
  }

  //
  // evaluation
  //
  void operator()(const state_type & state,
		  const independent_variable_type currentTime,
		  right_hand_side_type & V) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   V, nullptr);
  }

  void operator()(const state_type & state,
		 const independent_variable_type currentTime,
		 right_hand_side_type & V,
		 jacobian_type & jacobian,
		 bool computeJacobian) const
  {
    if (computeJacobian){
      T::velocityAndOptionalJacobian(state, currentTime,
				     V, &jacobian);
    }else{
      T::velocityAndOptionalJacobian(state, currentTime,
				     V, nullptr);
    }
  }

  void rightHandSide(const state_type & state,
		     const independent_variable_type currentTime,
		     right_hand_side_type & V) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   V, nullptr);
  }

  void rightHandSideAndJacobian(const state_type & state,
				const independent_variable_type currentTime,
				right_hand_side_type & V,
				jacobian_type & jacobian) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   V, &jacobian);
  }

  void jacobian(const state_type & state,
		const independent_variable_type currentTime,
		jacobian_type & jacobian) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_rhs, &jacobian);
  }

  void applyJacobian(const state_type & state,
		     const eigen_vec_t & operand,
		     const independent_variable_type currentTime,
		     eigen_vec_t & result) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_rhs, &m_jacobian);
    result = m_jacobian * operand;
  }

  void applyJacobian(const state_type & state,
		     const eigen_mat_ll_t & operand,
		     const independent_variable_type currentTime,
		     eigen_mat_ll_t & result) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_rhs, &m_jacobian);
    result = m_jacobian * operand;
  }

  void applyJacobian(const state_type & state,
		     const eigen_mat_lr_t & operand,
		     const independent_variable_type currentTime,
		     eigen_mat_lr_t & result) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_rhs, &m_jacobian);
    result = m_jacobian * operand;
  }

private:
  mutable jacobian_type m_jacobian;
  mutable right_hand_side_type m_rhs;
};

}
#endif
