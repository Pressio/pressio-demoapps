
#ifndef PRESSIODEMOAPPS_ADAPTER_MIXINS_HPP_
#define PRESSIODEMOAPPS_ADAPTER_MIXINS_HPP_

namespace pressiodemoapps{

template<class T>
class PublicProblemMixinCpp : public T
{
  using typename T::index_t;

public:
  using scalar_type   = typename T::scalar_type;
  using state_type    = typename T::state_type;
  using velocity_type = typename T::velocity_type;
  using jacobian_type = typename T::jacobian_type;

  template<class ...Args>
  PublicProblemMixinCpp(Args && ... args)
    : T(std::forward<Args>(args)...)
  {
    T::initializeJacobian(m_jacobian);
    ::pressiodemoapps::set_zero(m_jacobian);
    ::pressiodemoapps::resize(m_velocity, T::m_numDofSampleMesh);
    ::pressiodemoapps::set_zero(m_velocity);
  }

  index_t totalDofSampleMesh()  const{ return T::m_numDofSampleMesh;  }
  index_t totalDofStencilMesh() const{ return T::m_numDofStencilMesh; }

  velocity_type createVelocity() const{
    // Eigen has value semantics so this returns a deep copy
    ::pressiodemoapps::set_zero(m_velocity);
    return m_velocity;
  }

  jacobian_type createJacobian() const{
    // Eigen has value semantics so this returns a deep copy
    ::pressiodemoapps::set_zero(m_jacobian);
    return m_jacobian;
  }

  void velocity(const state_type & state,
		const scalar_type currentTime,
		velocity_type & V) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   V, nullptr);
  }

  void velocityAndJacobian(const state_type & state,
			   const scalar_type currentTime,
			   velocity_type & V,
			   jacobian_type & jacobian) const
  {
    T::velocityAndOptionalJacobian(state, currentTime, V, &jacobian);
  }

  void jacobian(const state_type & state,
		const scalar_type currentTime,
		jacobian_type & jacobian) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_velocity, &jacobian);
  }

private:
  mutable jacobian_type m_jacobian;
  mutable velocity_type m_velocity;
};


template<class T>
class PublicProblemMixinPy : public T
{

public:
  using scalar_type   = typename T::scalar_type;
  using state_type    = typename T::state_type;
  using velocity_type = typename T::velocity_type;
  using jacobian_type = typename T::jacobian_type;

private:
  using typename T::index_t;
  using eigen_vec_t = Eigen::Matrix<scalar_type, -1, 1>;
  using eigen_ref_vec_t = Eigen::Ref<eigen_vec_t>;
  using eigen_ref_const_vec_t = Eigen::Ref<const eigen_vec_t>;

  using eigen_mat_ll_t = Eigen::Matrix<scalar_type, -1, -1, Eigen::ColMajor>;
  using eigen_ref_mat_ll_t = Eigen::Ref<eigen_mat_ll_t>;
  using eigen_ref_const_mat_ll_t = Eigen::Ref<const eigen_mat_ll_t>;

  using eigen_mat_lr_t = Eigen::Matrix<scalar_type, -1, -1, Eigen::RowMajor>;
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

  index_t totalDofSampleMesh()  const{ return T::m_numDofSampleMesh;  }
  index_t totalDofStencilMesh() const{ return T::m_numDofStencilMesh; }

  velocity_type createVelocity() const{
    // Eigen has value semantics so this returns a deep copy
    ::pressiodemoapps::set_zero(m_velocity);
    return m_velocity;
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

  void velocity(const Eigen::Ref<const state_type> & state,
		const scalar_type currentTime,
		Eigen::Ref<velocity_type> & V) const
  {
    T::velocityAndOptionalJacobian(state, currentTime, V, nullptr);
  }

  void applyJacobianRank1(const Eigen::Ref<const state_type> & state,
			  const scalar_type currentTime,
			  const eigen_ref_const_vec_t & operand,
			  eigen_ref_vec_t & result) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_velocity, &m_jacobian);
    result = m_jacobian * operand;
  }

  void applyJacobianRank2_ll(const Eigen::Ref<const state_type> & state,
			     const scalar_type currentTime,
			     const eigen_ref_const_mat_ll_t & operand,
			     eigen_ref_mat_ll_t & result) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_velocity, &m_jacobian);
    result = m_jacobian * operand;
  }

  void applyJacobianRank2_lr(const Eigen::Ref<const state_type> & state,
			     const scalar_type currentTime,
			     const eigen_ref_const_mat_lr_t & operand,
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
