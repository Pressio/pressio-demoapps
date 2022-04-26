
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
  using scalar_type   = typename T::scalar_type;
  using state_type    = typename T::state_type;
  using velocity_type = typename T::velocity_type;
  using jacobian_type = typename T::jacobian_type;

private:
  using typename T::index_t;
  using eigen_vec_t    = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
  using eigen_mat_ll_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  using eigen_mat_lr_t = Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

public:
  template<class ...Args>
  PublicProblemEigenMixinCpp(Args && ... args)
    : T(std::forward<Args>(args)...)
  {
    T::initializeJacobian(m_jacobian);
    ::pressiodemoapps::set_zero(m_jacobian);
    ::pressiodemoapps::resize(m_velocity, T::m_numDofSampleMesh);
    ::pressiodemoapps::set_zero(m_velocity);
  }

  typename T::index_t totalDofSampleMesh()  const{
    return T::m_numDofSampleMesh;
  }

  typename T::index_t totalDofStencilMesh() const{
    return T::m_numDofStencilMesh;
  }

  velocity_type createVelocity() const{
    namespace pda = ::pressiodemoapps;
    return pda::clone(m_velocity);
  }

  jacobian_type createJacobian() const{
    namespace pda = ::pressiodemoapps;
    return pda::clone(m_jacobian);
  }

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
    T::velocityAndOptionalJacobian(state, currentTime,
				   V, &jacobian);
  }

  void jacobian(const state_type & state,
		const scalar_type currentTime,
		jacobian_type & jacobian) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_velocity, &jacobian);
  }

  void applyJacobian(const state_type & state,
		     const eigen_vec_t & operand,
		     const scalar_type currentTime,
		     eigen_vec_t & result) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_velocity, &m_jacobian);
    result = m_jacobian * operand;
  }

  void applyJacobian(const state_type & state,
		     const eigen_mat_ll_t & operand,
		     const scalar_type currentTime,
		     eigen_mat_ll_t & result) const
  {
    T::velocityAndOptionalJacobian(state, currentTime,
				   m_velocity, &m_jacobian);
    result = m_jacobian * operand;
  }

  void applyJacobian(const state_type & state,
		     const eigen_mat_lr_t & operand,
		     const scalar_type currentTime,
		     eigen_mat_lr_t & result) const
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
