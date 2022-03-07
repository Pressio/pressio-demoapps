
#ifndef PRESSIODEMOAPPS_ADAPTER_MIXINS_HPP_
#define PRESSIODEMOAPPS_ADAPTER_MIXINS_HPP_

namespace pressiodemoapps{

template<class T>
class PublicProblemMixinCpp : public T
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
  PublicProblemMixinCpp(Args && ... args)
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
		 "Parent class is not based on Eigen");

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
