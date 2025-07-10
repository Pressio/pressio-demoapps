/*
//@HEADER
// ************************************************************************
//
// adapter_py.hpp
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

#ifndef PRESSIODEMOAPPS_ADAPTER_MIXINS_PYTHON_HPP_
#define PRESSIODEMOAPPS_ADAPTER_MIXINS_PYTHON_HPP_

namespace pressiodemoapps{

// this is the class that gets instantiated and passed
// to the user when they create any problem.
// The specific problem class they want becomes the template T.
template<class T>
class PublicProblemMixinPy : public T
{

public:
  using scalar_type   = typename T::scalar_type;
  using state_type    = typename T::state_type;
  using velocity_type = typename T::velocity_type;
  using jacobian_type = typename T::jacobian_type;

  // need to assert that the parent is using Eigen types
  // because that is required for the bindings to work
  static_assert( ::pressiodemoapps::is_vector_eigen<state_type>::value
		 and ::pressiodemoapps::is_vector_eigen<velocity_type>::value
		 and ::pressiodemoapps::is_sparse_matrix_eigen<jacobian_type>::value,
		 "Parent class does not use Eigen types, this is required for bindings");

private:
  using eigen_vec_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using eigen_ref_vec_t = Eigen::Ref<eigen_vec_t>;
  using eigen_ref_const_vec_t = Eigen::Ref<const eigen_vec_t>;

  using eigen_mat_ll_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  using eigen_ref_mat_ll_t = Eigen::Ref<eigen_mat_ll_t>;
  using eigen_ref_const_mat_ll_t = Eigen::Ref<const eigen_mat_ll_t>;

  using eigen_mat_lr_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using eigen_ref_mat_lr_t = Eigen::Ref<eigen_mat_lr_t>;
  using eigen_ref_const_mat_lr_t = Eigen::Ref<const eigen_mat_lr_t>;

public:
  template<class ...Args>
  PublicProblemMixinPy(Args && ... args)
    : T(std::forward<Args>(args)...)
  {
    T::initializeJacobian(m_jacobian);
    ::pressiodemoapps::set_zero(m_jacobian);
    ::pressiodemoapps::resize(m_velocity, T::m_numDofSampleMesh);
    ::pressiodemoapps::set_zero(m_velocity);
  }

  int numDofPerCell() const {
    return T::numDofPerCellImpl();
  }

  typename T::index_t totalDofSampleMesh()  const{
    return T::m_numDofSampleMesh;
  }

  typename T::index_t totalDofStencilMesh() const{
    return T::m_numDofStencilMesh;
  }

  /* create methods
     For creating objects for Python, we need to use the
     native types from the parent NOT Eigen::Ref<>.
     This is because we want a copy to be made and returned to Python.
     For example by doing:
	velocity_type createVelocity() const

     velocity_type is an Eigen vector, so via the binding this
     will return a numpy array to Python.
   */

  velocity_type createVelocity() const{
    namespace pda = ::pressiodemoapps;
    return pda::clone(m_velocity);
  }

  eigen_vec_t
  createApplyJacobianResultRank1(const eigen_ref_const_vec_t & operand) const{
    eigen_vec_t res(m_jacobian.rows(), 1);
    ::pressiodemoapps::set_zero(res);
    return res;
  }

  eigen_mat_ll_t
  createApplyJacobianResultRank2_ll(const eigen_ref_const_mat_ll_t & operand) const{
    eigen_mat_ll_t res(m_jacobian.rows(), operand.cols());
    ::pressiodemoapps::set_zero(res);
    return res;
  }

  eigen_mat_lr_t
  createApplyJacobianResultRank2_lr(const eigen_ref_const_mat_lr_t & operand) const{
    eigen_mat_lr_t res(m_jacobian.rows(), operand.cols());
    ::pressiodemoapps::set_zero(res);
    return res;
  }

  /* for the arguments, we need to use the proper Eigen::Ref<>.
     This is because we want arguments to "reference" the
     corresponding numpy array in Python not copy it every time.
     So using Ref<> ensures we have view semantics not value semantics.
   */

  void velocity(const Eigen::Ref<const state_type> & state,
		const scalar_type currentTime,
		Eigen::Ref<velocity_type> & V) const
  {
    T::velocityAndOptionalJacobian(state, currentTime, V, nullptr);
  }

  // note that it does not matter these methods are called Rank1, Rank2, etc
  // because this is just the C++ name. We need different names because
  // these do different things but when we bind these we create a single
  // python name "applyJacobian" to exposes all these.
  // It has multiple overloads so it will
  // call the right overload based on the numpy array passed.
  void applyJacobianRank1(const Eigen::Ref<const state_type> & state,
			  const eigen_ref_const_vec_t & operand,
			  const scalar_type currentTime,
			  eigen_ref_vec_t & result) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_velocity, &m_jacobian);
    result = m_jacobian * operand;
  }

  void applyJacobianRank2_ll(const Eigen::Ref<const state_type> & state,
			     const eigen_ref_const_mat_ll_t & operand,
			     const scalar_type currentTime,
			     eigen_ref_mat_ll_t & result) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_velocity, &m_jacobian);
    result = m_jacobian * operand;
  }

  void applyJacobianRank2_lr(const Eigen::Ref<const state_type> & state,
			     const eigen_ref_const_mat_lr_t & operand,
			     const scalar_type currentTime,
			     eigen_ref_mat_lr_t & result) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_velocity, &m_jacobian);
    result = m_jacobian * operand;
  }

private:
  mutable jacobian_type m_jacobian;
  mutable velocity_type m_velocity;
};

}
#endif
