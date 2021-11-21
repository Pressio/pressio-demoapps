
#ifndef PRESSIODEMOAPPS_EE_FLUX_FUNCTOR_HPP_
#define PRESSIODEMOAPPS_EE_FLUX_FUNCTOR_HPP_

namespace pressiodemoapps{ namespace ee{ namespace impl{

//
// FLUX VALUES
//

template<class Parent, int ndpc, class scalar_type, class flux_type>
struct ComputeDirectionalFluxValues;

// for 3 dofs, so Euler 1d
template<class Parent, class scalar_type, class flux_type>
struct ComputeDirectionalFluxValues<Parent, 3, scalar_type, flux_type> : Parent{
private:
  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  flux_type & m_fluxL;
  flux_type & m_fluxR;

public:
  template<class ...Args>
  ComputeDirectionalFluxValues(InviscidFluxScheme fluxEnum, scalar_type gamma,
			       flux_type & fluxL, flux_type & fluxR, Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_fluxEnum(fluxEnum), m_gamma(gamma),
      m_fluxL(fluxL), m_fluxR(fluxR)
  {}

  const flux_type & fluxLeft()  const { return m_fluxL; }
  const flux_type & fluxRight() const { return m_fluxR; }

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    Parent::operator()(smPt);
    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();
    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      ee::impl::ee_rusanov_flux_three_dof(m_fluxL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
      ee::impl::ee_rusanov_flux_three_dof(m_fluxR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);
      break;
    }
  }
};

// for 4 dofs, so Euler 2d
template<class Parent, class scalar_type, class flux_type>
struct ComputeDirectionalFluxValues<Parent, 4, scalar_type, flux_type> : Parent {
private:
  const std::array<scalar_type, 2> m_normal = {};
  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  flux_type & m_fluxL;
  flux_type & m_fluxR;

public:
  template<class ...Args>
  ComputeDirectionalFluxValues(InviscidFluxScheme fluxEnum,
			       const std::array<scalar_type, 2> normal,
			       scalar_type gamma, flux_type & fluxL, flux_type & fluxR,
			       Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_normal(normal), m_fluxEnum(fluxEnum), m_gamma(gamma),
      m_fluxL(fluxL), m_fluxR(fluxR)
  {}

  const flux_type & fluxLeft()  const { return m_fluxL; }
  const flux_type & fluxRight() const { return m_fluxR; }

  template<class IndexType>
  void operator()(IndexType smPt){
    Parent::operator()(smPt);
    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();
    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      ee::impl::ee_rusanov_flux_four_dof(m_fluxL, uMinusHalfNeg, uMinusHalfPos,
				     m_normal, m_gamma);
      ee::impl::ee_rusanov_flux_four_dof(m_fluxR, uPlusHalfNeg,  uPlusHalfPos,
				     m_normal, m_gamma);
      break;
    }
  }
};


// for 5 dofs, so Euler 3d
template<class Parent, class scalar_type, class flux_type>
struct ComputeDirectionalFluxValues<Parent, 5, scalar_type, flux_type> : Parent {
private:
  const std::array<scalar_type, 3> m_normal = {};
  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  flux_type & m_fluxL;
  flux_type & m_fluxR;

public:
  template<class ...Args>
  ComputeDirectionalFluxValues(InviscidFluxScheme fluxEnum,
			       const std::array<scalar_type, 3> normal,
			       scalar_type gamma, flux_type & fluxL, flux_type & fluxR,
			       Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_normal(normal), m_fluxEnum(fluxEnum), m_gamma(gamma),
      m_fluxL(fluxL), m_fluxR(fluxR)
  {}

  const flux_type & fluxLeft()  const { return m_fluxL; }
  const flux_type & fluxRight() const { return m_fluxR; }

  template<class IndexType>
  void operator()(IndexType smPt){
    Parent::operator()(smPt);
    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();
    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      ee::impl::ee_rusanov_flux_five_dof(m_fluxL, uMinusHalfNeg, uMinusHalfPos,
					 m_normal, m_gamma);
      ee::impl::ee_rusanov_flux_five_dof(m_fluxR, uPlusHalfNeg,  uPlusHalfPos,
					 m_normal, m_gamma);
      break;
    }
  }
};



//
// FLUX JACOBIANS
//

template<class Parent, int ndpc, class scalar_type, class flux_jac_type>
struct ComputeDirectionalFluxJacobians;

// for 3 dofs, so Euler 1d
template<class Parent, class scalar_type, class flux_jac_type>
struct ComputeDirectionalFluxJacobians<Parent, 3, scalar_type, flux_jac_type> : Parent {
private:
  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;

public:
  template<class ...Args>
  ComputeDirectionalFluxJacobians(InviscidFluxScheme fluxEnum,
				  scalar_type gamma,
				  flux_jac_type & fluxJacLNeg,
				  flux_jac_type & fluxJacLPos,
				  flux_jac_type & fluxJacRNeg,
				  flux_jac_type & fluxJacRPos,
				  Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_fluxEnum(fluxEnum), m_gamma(gamma),
      m_fluxJacLNeg(fluxJacLNeg), m_fluxJacLPos(fluxJacLPos),
      m_fluxJacRNeg(fluxJacRNeg), m_fluxJacRPos(fluxJacRPos)
  {}

  const flux_jac_type & fluxJacLNeg()  const { return m_fluxJacLNeg; }
  const flux_jac_type & fluxJacLPos()  const { return m_fluxJacLPos; }
  const flux_jac_type & fluxJacRNeg()  const { return m_fluxJacRNeg; }
  const flux_jac_type & fluxJacRPos()  const { return m_fluxJacRPos; }

  template<class IndexType>
  void operator()(IndexType smPt){
    Parent::operator()(smPt);
    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();
    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      ee::impl::ee_rusanov_flux_jacobian_three_dof(m_fluxJacLNeg, m_fluxJacLPos,
					      uMinusHalfNeg, uMinusHalfPos,
					      m_gamma);
      ee::impl::ee_rusanov_flux_jacobian_three_dof(m_fluxJacRNeg, m_fluxJacRPos,
					      uPlusHalfNeg, uPlusHalfPos,
					      m_gamma);
      break;
    }
  }
};

// for 4 dofs, so Euler 2d
template<class Parent, class scalar_type, class flux_jac_type>
struct ComputeDirectionalFluxJacobians<Parent, 4, scalar_type, flux_jac_type> : Parent {
private:
  const std::array<scalar_type, 2> m_normal = {};
  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;

public:
  template<class ...Args>
  ComputeDirectionalFluxJacobians(InviscidFluxScheme fluxEnum,
				  const std::array<scalar_type, 2> normal,
				  scalar_type gamma,
				  flux_jac_type & fluxJacLNeg,
				  flux_jac_type & fluxJacLPos,
				  flux_jac_type & fluxJacRNeg,
				  flux_jac_type & fluxJacRPos,
				  Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_normal(normal),
      m_fluxEnum(fluxEnum), m_gamma(gamma),
      m_fluxJacLNeg(fluxJacLNeg), m_fluxJacLPos(fluxJacLPos),
      m_fluxJacRNeg(fluxJacRNeg), m_fluxJacRPos(fluxJacRPos)
  {}

  const flux_jac_type & fluxJacLNeg()  const { return m_fluxJacLNeg; }
  const flux_jac_type & fluxJacLPos()  const { return m_fluxJacLPos; }
  const flux_jac_type & fluxJacRNeg()  const { return m_fluxJacRNeg; }
  const flux_jac_type & fluxJacRPos()  const { return m_fluxJacRPos; }

  template<class IndexType>
  void operator()(IndexType smPt){
    Parent::operator()(smPt);
    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();
    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      ee::impl::ee_rusanov_flux_jacobian_four_dof(m_fluxJacLNeg, m_fluxJacLPos,
					     uMinusHalfNeg, uMinusHalfPos,
					     m_normal, m_gamma);
      ee::impl::ee_rusanov_flux_jacobian_four_dof(m_fluxJacRNeg, m_fluxJacRPos,
					     uPlusHalfNeg, uPlusHalfPos,
					     m_normal, m_gamma);
      break;
    }
  }
};

// for 5 dofs, so Euler 3d
template<class Parent, class scalar_type, class flux_jac_type>
struct ComputeDirectionalFluxJacobians<Parent, 5, scalar_type, flux_jac_type> : Parent {
private:
  const std::array<scalar_type, 3> m_normal = {};
  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;

public:
  template<class ...Args>
  ComputeDirectionalFluxJacobians(InviscidFluxScheme fluxEnum,
				  const std::array<scalar_type, 3> normal,
				  scalar_type gamma,
				  flux_jac_type & fluxJacLNeg,
				  flux_jac_type & fluxJacLPos,
				  flux_jac_type & fluxJacRNeg,
				  flux_jac_type & fluxJacRPos,
				  Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_normal(normal),
      m_fluxEnum(fluxEnum), m_gamma(gamma),
      m_fluxJacLNeg(fluxJacLNeg), m_fluxJacLPos(fluxJacLPos),
      m_fluxJacRNeg(fluxJacRNeg), m_fluxJacRPos(fluxJacRPos)
  {}

  const flux_jac_type & fluxJacLNeg()  const { return m_fluxJacLNeg; }
  const flux_jac_type & fluxJacLPos()  const { return m_fluxJacLPos; }
  const flux_jac_type & fluxJacRNeg()  const { return m_fluxJacRNeg; }
  const flux_jac_type & fluxJacRPos()  const { return m_fluxJacRPos; }

  template<class IndexType>
  void operator()(IndexType smPt){
    Parent::operator()(smPt);
    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();
    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      ee::impl::ee_rusanov_flux_jacobian_five_dof(m_fluxJacLNeg, m_fluxJacLPos,
						  uMinusHalfNeg, uMinusHalfPos,
						  m_normal, m_gamma);
      ee::impl::ee_rusanov_flux_jacobian_five_dof(m_fluxJacRNeg, m_fluxJacRPos,
						  uPlusHalfNeg, uPlusHalfPos,
						  m_normal, m_gamma);
      break;
    }
  }
};


//
// FLUX VALUES AND JACOBIANS
//

template<class Parent, int ndpc, class scalar_type, class flux_type, class flux_jac_type>
struct ComputeDirectionalFluxValuesAndJacobians;

// for 3 dofs, so Euler 1d
template<class Parent, class scalar_type, class flux_type, class flux_jac_type>
struct ComputeDirectionalFluxValuesAndJacobians<
  Parent, 3, scalar_type, flux_type, flux_jac_type> : Parent
{
private:
  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  flux_type & m_fluxL;
  flux_type & m_fluxR;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;

public:
  template<class ...Args>
  ComputeDirectionalFluxValuesAndJacobians(InviscidFluxScheme fluxEnum,
					   scalar_type gamma,
					   flux_type & fluxL,
					   flux_type & fluxR,
					   flux_jac_type & fluxJacLNeg,
					   flux_jac_type & fluxJacLPos,
					   flux_jac_type & fluxJacRNeg,
					   flux_jac_type & fluxJacRPos,
					   Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_fluxEnum(fluxEnum), m_gamma(gamma),
      m_fluxL(fluxL), m_fluxR(fluxR),
      m_fluxJacLNeg(fluxJacLNeg), m_fluxJacLPos(fluxJacLPos),
      m_fluxJacRNeg(fluxJacRNeg), m_fluxJacRPos(fluxJacRPos)
  {}

  const flux_type & fluxLeft()  const { return m_fluxL; }
  const flux_type & fluxRight() const { return m_fluxR; }
  const flux_jac_type & fluxJacLNeg()  const { return m_fluxJacLNeg; }
  const flux_jac_type & fluxJacLPos()  const { return m_fluxJacLPos; }
  const flux_jac_type & fluxJacRNeg()  const { return m_fluxJacRNeg; }
  const flux_jac_type & fluxJacRPos()  const { return m_fluxJacRPos; }

  template<class IndexType>
  void operator()(IndexType smPt){
    Parent::operator()(smPt);
    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      ee::impl::ee_rusanov_flux_three_dof(m_fluxL, uMinusHalfNeg,
				      uMinusHalfPos, m_gamma);
      ee::impl::ee_rusanov_flux_three_dof(m_fluxR, uPlusHalfNeg,
				      uPlusHalfPos,  m_gamma);
      ee::impl::ee_rusanov_flux_jacobian_three_dof(m_fluxJacLNeg, m_fluxJacLPos,
					      uMinusHalfNeg, uMinusHalfPos,
					      m_gamma);
      ee::impl::ee_rusanov_flux_jacobian_three_dof(m_fluxJacRNeg, m_fluxJacRPos,
					      uPlusHalfNeg, uPlusHalfPos,
					      m_gamma);
      break;
    }
  }
};

// for 4 dofs, so Euler 2d
template<class Parent, class scalar_type, class flux_type, class flux_jac_type>
struct ComputeDirectionalFluxValuesAndJacobians<
  Parent, 4, scalar_type, flux_type, flux_jac_type> : Parent
{
private:
  const std::array<scalar_type, 2> m_normal = {};
  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  flux_type & m_fluxL;
  flux_type & m_fluxR;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;

public:
  template<class ...Args>
  ComputeDirectionalFluxValuesAndJacobians(InviscidFluxScheme fluxEnum,
					   const std::array<scalar_type, 2> normal,
					   scalar_type gamma,
					   flux_type & fluxL,
					   flux_type & fluxR,
					   flux_jac_type & fluxJacLNeg,
					   flux_jac_type & fluxJacLPos,
					   flux_jac_type & fluxJacRNeg,
					   flux_jac_type & fluxJacRPos,
					   Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_normal(normal),
      m_fluxEnum(fluxEnum), m_gamma(gamma),
      m_fluxL(fluxL), m_fluxR(fluxR),
      m_fluxJacLNeg(fluxJacLNeg), m_fluxJacLPos(fluxJacLPos),
      m_fluxJacRNeg(fluxJacRNeg), m_fluxJacRPos(fluxJacRPos)
  {}

  const flux_type & fluxLeft()  const { return m_fluxL; }
  const flux_type & fluxRight() const { return m_fluxR; }
  const flux_jac_type & fluxJacLNeg()  const { return m_fluxJacLNeg; }
  const flux_jac_type & fluxJacLPos()  const { return m_fluxJacLPos; }
  const flux_jac_type & fluxJacRNeg()  const { return m_fluxJacRNeg; }
  const flux_jac_type & fluxJacRPos()  const { return m_fluxJacRPos; }

  template<class IndexType>
  void operator()(IndexType smPt){
    Parent::operator()(smPt);
    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      ee::impl::ee_rusanov_flux_four_dof(m_fluxL, uMinusHalfNeg,
				     uMinusHalfPos, m_normal, m_gamma);
      ee::impl::ee_rusanov_flux_four_dof(m_fluxR, uPlusHalfNeg,
				     uPlusHalfPos,  m_normal, m_gamma);
      ee::impl::ee_rusanov_flux_jacobian_four_dof(m_fluxJacLNeg, m_fluxJacLPos,
					     uMinusHalfNeg, uMinusHalfPos,
					     m_normal, m_gamma);
      ee::impl::ee_rusanov_flux_jacobian_four_dof(m_fluxJacRNeg, m_fluxJacRPos,
					     uPlusHalfNeg, uPlusHalfPos,
					     m_normal, m_gamma);
      break;
    }
  }
};

// for 5 dofs, so Euler 3d
template<class Parent, class scalar_type, class flux_type, class flux_jac_type>
struct ComputeDirectionalFluxValuesAndJacobians<
  Parent, 5, scalar_type, flux_type, flux_jac_type> : Parent
{
private:
  const std::array<scalar_type, 3> m_normal = {};
  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  flux_type & m_fluxL;
  flux_type & m_fluxR;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;

public:
  template<class ...Args>
  ComputeDirectionalFluxValuesAndJacobians(InviscidFluxScheme fluxEnum,
					   const std::array<scalar_type, 3> normal,
					   scalar_type gamma,
					   flux_type & fluxL,
					   flux_type & fluxR,
					   flux_jac_type & fluxJacLNeg,
					   flux_jac_type & fluxJacLPos,
					   flux_jac_type & fluxJacRNeg,
					   flux_jac_type & fluxJacRPos,
					   Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_normal(normal),
      m_fluxEnum(fluxEnum), m_gamma(gamma),
      m_fluxL(fluxL), m_fluxR(fluxR),
      m_fluxJacLNeg(fluxJacLNeg), m_fluxJacLPos(fluxJacLPos),
      m_fluxJacRNeg(fluxJacRNeg), m_fluxJacRPos(fluxJacRPos)
  {}

  const flux_type & fluxLeft()  const { return m_fluxL; }
  const flux_type & fluxRight() const { return m_fluxR; }
  const flux_jac_type & fluxJacLNeg()  const { return m_fluxJacLNeg; }
  const flux_jac_type & fluxJacLPos()  const { return m_fluxJacLPos; }
  const flux_jac_type & fluxJacRNeg()  const { return m_fluxJacRNeg; }
  const flux_jac_type & fluxJacRPos()  const { return m_fluxJacRPos; }

  template<class IndexType>
  void operator()(IndexType smPt){
    Parent::operator()(smPt);
    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      ee::impl::ee_rusanov_flux_five_dof(m_fluxL, uMinusHalfNeg,
					 uMinusHalfPos, m_normal, m_gamma);
      ee::impl::ee_rusanov_flux_five_dof(m_fluxR, uPlusHalfNeg,
					 uPlusHalfPos,  m_normal, m_gamma);
      ee::impl::ee_rusanov_flux_jacobian_five_dof(m_fluxJacLNeg, m_fluxJacLPos,
						  uMinusHalfNeg, uMinusHalfPos,
						  m_normal, m_gamma);
      ee::impl::ee_rusanov_flux_jacobian_five_dof(m_fluxJacRNeg, m_fluxJacRPos,
						  uPlusHalfNeg, uPlusHalfPos,
						  m_normal, m_gamma);
      break;
    }
  }
};

}}}
#endif
