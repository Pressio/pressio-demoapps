
#ifndef PRESSIODEMOAPPS_ADVECTION1d_FLUX_FUNCTOR_HPP_
#define PRESSIODEMOAPPS_ADVECTION1d_FLUX_FUNCTOR_HPP_

namespace pressiodemoapps{ namespace impladvection1d{

template<class Parent, class scalar_type, class flux_type>
struct ComputeDirectionalFluxValues : Parent
{
private:
  InviscidFluxScheme m_fluxEnum;
  flux_type & m_fluxL;
  flux_type & m_fluxR;
  scalar_type m_advectionVelocity;

public:
  template<class ...Args>
  ComputeDirectionalFluxValues(InviscidFluxScheme fluxEnum,
			       flux_type & fluxL, flux_type & fluxR,
			       scalar_type advectionVelocity,
			       Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_fluxEnum(fluxEnum), m_fluxL(fluxL), m_fluxR(fluxR),
      m_advectionVelocity(advectionVelocity)
  {}

  const flux_type & fluxLeft()  const { return m_fluxL; }
  const flux_type & fluxRight() const { return m_fluxR; }

  template<class IndexType>
  void operator()(IndexType smPt, int ndpc)
  {
    Parent::operator()(smPt, ndpc);

    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      m_fluxL(0) = uMinusHalfNeg(0)*m_advectionVelocity;
      m_fluxR(0) = uPlusHalfNeg(0)*m_advectionVelocity;
      break;
    }
  }
};

template<class Parent, class scalar_type, class flux_jac_type>
struct ComputeDirectionalFluxJacobians : Parent
{
private:
  InviscidFluxScheme m_fluxEnum;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;
  scalar_type m_advectionVelocity;

public:
  template<class ...Args>
  ComputeDirectionalFluxJacobians(InviscidFluxScheme fluxEnum,
				  flux_jac_type & fluxJacLNeg,
				  flux_jac_type & fluxJacLPos,
				  flux_jac_type & fluxJacRNeg,
				  flux_jac_type & fluxJacRPos,
				  scalar_type advectionVelocity,
				  Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_fluxEnum(fluxEnum),
      m_fluxJacLNeg(fluxJacLNeg), m_fluxJacLPos(fluxJacLPos),
      m_fluxJacRNeg(fluxJacRNeg), m_fluxJacRPos(fluxJacRPos),
      m_advectionVelocity(advectionVelocity)
  {}

  const flux_jac_type & fluxJacLNeg()  const { return m_fluxJacLNeg; }
  const flux_jac_type & fluxJacLPos()  const { return m_fluxJacLPos; }
  const flux_jac_type & fluxJacRNeg()  const { return m_fluxJacRNeg; }
  const flux_jac_type & fluxJacRPos()  const { return m_fluxJacRPos; }

  template<class IndexType>
  void operator()(IndexType smPt, int ndpc)
  {
    Parent::operator()(smPt, ndpc);

    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      m_fluxJacLNeg(0,0) = m_advectionVelocity;
      m_fluxJacLPos(0,0) = 0;
      m_fluxJacRNeg(0,0) = m_advectionVelocity;
      m_fluxJacRPos(0,0) = 0;
      break;
    }
  }
};

template<class Parent, class scalar_type, class flux_type, class flux_jac_type>
struct ComputeDirectionalFluxValuesAndJacobians : Parent
{
private:
  InviscidFluxScheme m_fluxEnum;
  flux_type & m_fluxL;
  flux_type & m_fluxR;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;
  scalar_type m_advectionVelocity;

public:
  template<class ...Args>
  ComputeDirectionalFluxValuesAndJacobians(InviscidFluxScheme fluxEnum,
					   flux_type & fluxL,
					   flux_type & fluxR,
					   flux_jac_type & fluxJacLNeg,
					   flux_jac_type & fluxJacLPos,
					   flux_jac_type & fluxJacRNeg,
					   flux_jac_type & fluxJacRPos,
					   scalar_type advectionVelocity,
					   Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_fluxEnum(fluxEnum),
      m_fluxL(fluxL), m_fluxR(fluxR),
      m_fluxJacLNeg(fluxJacLNeg), m_fluxJacLPos(fluxJacLPos),
      m_fluxJacRNeg(fluxJacRNeg), m_fluxJacRPos(fluxJacRPos),
      m_advectionVelocity(advectionVelocity)
  {}

  const flux_type & fluxLeft()  const { return m_fluxL; }
  const flux_type & fluxRight() const { return m_fluxR; }
  const flux_jac_type & fluxJacLNeg()  const { return m_fluxJacLNeg; }
  const flux_jac_type & fluxJacLPos()  const { return m_fluxJacLPos; }
  const flux_jac_type & fluxJacRNeg()  const { return m_fluxJacRNeg; }
  const flux_jac_type & fluxJacRPos()  const { return m_fluxJacRPos; }

  template<class IndexType>
  void operator()(IndexType smPt, int ndpc)
  {
    Parent::operator()(smPt, ndpc);

    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      m_fluxL(0)    = uMinusHalfNeg(0)*m_advectionVelocity;
      m_fluxR(0)    = uPlusHalfNeg(0)*m_advectionVelocity;
      m_fluxJacLNeg(0,0) = m_advectionVelocity;
      m_fluxJacLPos(0,0) = 0;
      m_fluxJacRNeg(0,0) = m_advectionVelocity;
      m_fluxJacRPos(0,0) = 0;
      break;
    }
  }
};

}}
#endif
