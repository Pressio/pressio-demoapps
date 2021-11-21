
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
  using eigen_mat_t = Eigen::Matrix<scalar_type, -1, -1>;
  using eigen_ref_vec_t = Eigen::Ref<eigen_vec_t>;
  using eigen_ref_mat_t = Eigen::Ref<eigen_mat_t>;
  using eigen_ref_const_vec_t = Eigen::Ref<const eigen_vec_t>;
  using eigen_ref_const_mat_t = Eigen::Ref<const eigen_mat_t>;

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
    return res;
  }

  eigen_mat_t
  createApplyJacobianResultRank2(const eigen_ref_const_mat_t & operand) const{
    eigen_mat_t res(m_jacobian.rows(), operand.cols());
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

  void applyJacobianRank2(const Eigen::Ref<const state_type> & state,
			  const scalar_type currentTime,
			  const eigen_ref_const_mat_t & operand,
			  eigen_ref_mat_t & result) const
  {

    T::velocityAndOptionalJacobian(state, currentTime,
				   m_velocity, &m_jacobian);
    result = m_jacobian * operand;
  }

private:
  mutable jacobian_type m_jacobian;
  mutable velocity_type m_velocity;
};



// template<class T>
// class ExplicitAdapterMixin : public T
// {
//   using typename T::index_t;

// public:
//   using scalar_type   = typename T::scalar_type;
//   using state_type    = typename T::state_type;
//   using velocity_type = typename T::velocity_type;

//   template<class ...Args>
//   ExplicitAdapterMixin(Args && ... args)
//     : T(std::forward<Args>(args)...)
//   {}

//   index_t totalDofSampleMesh()  const{ return T::m_numDofSampleMesh;  }
//   index_t totalDofStencilMesh() const{ return T::m_numDofStencilMesh; }

//   velocity_type createVelocity() const{
//     velocity_type v(T::m_numDofSampleMesh);
//     ::pressiodemoapps::set_zero(v);
//     return v;
//   }

//   // since this mixin is used both in C++ and bindings,
//   // we need to ensure the right types are used such
//   // for bindings we get reference semantics with Python numpy
//   // see https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html
//   void velocity(
// #ifdef PRESSIODEMOAPPS_ENABLE_BINDINGS
// 		const Eigen::Ref<const state_type> & state,
// 		const scalar_type currentTime,
// 		Eigen::Ref<velocity_type> & V
// #else
// 		const state_type & state,
// 		const scalar_type currentTime,
// 		velocity_type & V
// #endif
// 		) const
//   {
//     T::velocityAndOptionalJacobian(state, currentTime, V, nullptr);
//   }
// };

// template<class T>
// class ImplicitAdapterMixinCpp : public T
// {
//   using typename T::index_t;

// public:
//   using scalar_type   = typename T::scalar_type;
//   using state_type    = typename T::state_type;
//   using velocity_type = typename T::velocity_type;
//   using jacobian_type = typename T::jacobian_type;

//   template<class ...Args>
//   ImplicitAdapterMixinCpp(Args && ... args)
//     : T(std::forward<Args>(args)...)
//   {
//     T::initializeJacobian(m_jacobian);
//     ::pressiodemoapps::set_zero(m_jacobian);
//   }

//   index_t totalDofSampleMesh()  const{ return T::m_numDofSampleMesh;  }
//   index_t totalDofStencilMesh() const{ return T::m_numDofStencilMesh; }

//   velocity_type createVelocity() const{
//     velocity_type v(T::m_numDofSampleMesh);
//     ::pressiodemoapps::set_zero(v);
//     return v;
//   }

//   jacobian_type createJacobian() const{
//     // Eigen has value semantics so this returns a deep copy
//     return m_jacobian;
//   }

//   // here we are in C++, so just use the types we know for arguments
//   void velocity(const state_type & state,
// 		const scalar_type currentTime,
// 		velocity_type & V) const
//   {
//     m_lastRecordedTime = currentTime;
//     T::velocityAndOptionalJacobian(state, currentTime, V, &m_jacobian);
//   }

//   // here we are in C++, so just use the types we know for arguments
//   void jacobian(const state_type & state,
// 		const scalar_type currentTime,
// 		jacobian_type & jacobian) const
//   {
//     if (currentTime != m_lastRecordedTime or m_neverCalled){
//       auto tmpv = createVelocity();
//       T::velocityAndOptionalJacobian(state, currentTime, tmpv, &m_jacobian);
//       m_lastRecordedTime = currentTime;
//       m_neverCalled = false;
//     }

//     jacobian = m_jacobian;
//   }

// private:
//   mutable jacobian_type m_jacobian;

//   mutable bool m_neverCalled = true;

//   /* since here the velocity and Jacobian are fused, we need
//      some logic to ensure that if applyJacobian is called,
//      the jacobian has already bee computed or the applyJacobian
//      kernel will not make any sense.
//      For now, check if the last computed time matches.
//    */
//   mutable scalar_type m_lastRecordedTime = {};
// };


// // see https://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html
// template<class T>
// class ImplicitAdapterMixinPy : public T
// {
//   using typename T::index_t;

// public:
//   using scalar_type   = typename T::scalar_type;
//   using state_type    = typename T::state_type;
//   using velocity_type = typename T::velocity_type;
//   using jacobian_type = typename T::jacobian_type;

//   using eigen_vec_t = Eigen::Matrix<scalar_type, -1, 1>;
//   using eigen_mat_t = Eigen::Matrix<scalar_type, -1, -1>;

//   template<class ...Args>
//   ImplicitAdapterMixinPy(Args && ... args)
//     :T(std::forward<Args>(args)...)
//   {
//     T::initializeJacobian(m_jacobian);
//     ::pressiodemoapps::set_zero(m_jacobian);
//   }

//   index_t totalDofSampleMesh()  const{ return T::m_numDofSampleMesh;  }
//   index_t totalDofStencilMesh() const{ return T::m_numDofStencilMesh; }

//   velocity_type createVelocity() const{
//     velocity_type v(T::m_numDofSampleMesh);
//     ::pressiodemoapps::set_zero(v);
//     return v;
//   }

//   eigen_vec_t
//   createApplyJacobianResultRank1(const Eigen::Ref<const eigen_vec_t> & operand) const
//   {
//     eigen_vec_t res(m_jacobian.rows(), 1);
//     return res;
//   }

//   eigen_mat_t
//   createApplyJacobianResultRank2(const Eigen::Ref<const eigen_mat_t> & operand) const
//   {
//     eigen_mat_t res(m_jacobian.rows(), operand.cols());
//     return res;
//   }

//   void velocity(const Eigen::Ref<const state_type> & state,
// 		const scalar_type currentTime,
// 		Eigen::Ref<velocity_type> & V) const
//   {
//     m_lastRecordedTime = currentTime;
//     T::velocityAndOptionalJacobian(state, currentTime, V, &m_jacobian);
//   }

//   void applyJacobianRank1(const Eigen::Ref<const state_type> & state,
// 			  const scalar_type currentTime,
// 			  const Eigen::Ref< const eigen_vec_t> & operand,
// 			  Eigen::Ref< eigen_vec_t> & result) const
//   {
//     if (currentTime != m_lastRecordedTime or m_neverCalled){
//       auto tmpv = createVelocity();
//       T::velocityAndOptionalJacobian(state, currentTime, tmpv, &m_jacobian);
//       m_lastRecordedTime = currentTime;
//       m_neverCalled = false;
//     }

//     // this relies on jacobian being computed within velocity
//     result = m_jacobian * operand;
//   }

//   void applyJacobianRank2(const Eigen::Ref<const state_type> & state,
// 			  const scalar_type currentTime,
// 			  const Eigen::Ref<const eigen_mat_t> & operand,
// 			  Eigen::Ref<eigen_mat_t> & result) const
//   {

//     if (currentTime != m_lastRecordedTime or m_neverCalled){
//       auto tmpv = createVelocity();
//       T::velocityAndOptionalJacobian(state, currentTime, tmpv, &m_jacobian);
//       m_lastRecordedTime = currentTime;
//       m_neverCalled = false;
//     }

//     // this relies on jacobian being computed within velocity
//     result = m_jacobian * operand;
//   }

// private:
//   mutable jacobian_type m_jacobian;

//   mutable bool m_neverCalled = true;

//   /* since here the velocity and Jacobian are fused, we need
//      some logic to ensure that if applyJacobian is called,
//      the jacobian has already bee computed or the applyJacobian
//      kernel will not make any sense.
//      For now, check if the last computed time matches.
//    */
//   mutable scalar_type m_lastRecordedTime = {};
// };

}
#endif
