
#ifndef PRESSIODEMOAPPS_ADVECTION1d_FLUX_FUNCTOR_HPP_
#define PRESSIODEMOAPPS_ADVECTION1d_FLUX_FUNCTOR_HPP_

namespace pressiodemoapps{ namespace impladv{

template<class Parent, class scalar_type, class flux_type>
struct ComputeDirectionalFluxValues : Parent
{
private:
  InviscidFluxScheme m_fluxEnum;
  flux_type & m_fluxL;
  flux_type & m_fluxR;

public:
  template<class ...Args>
  ComputeDirectionalFluxValues(InviscidFluxScheme fluxEnum,
			       flux_type & fluxL, flux_type & fluxR,
			       Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_fluxEnum(fluxEnum), m_fluxL(fluxL), m_fluxR(fluxR)
  {}

  const flux_type & fluxLeft()  const { return m_fluxL; }
  const flux_type & fluxRight() const { return m_fluxR; }

  template<class IndexType>
  void operator()(IndexType smPt) {
    Parent::operator()(smPt);
    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      m_fluxL = uMinusHalfNeg;
      m_fluxR = uPlusHalfNeg;
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

public:
  template<class ...Args>
  ComputeDirectionalFluxJacobians(InviscidFluxScheme fluxEnum,
				  flux_jac_type & fluxJacLNeg,
				  flux_jac_type & fluxJacLPos,
				  flux_jac_type & fluxJacRNeg,
				  flux_jac_type & fluxJacRPos,
				  Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_fluxEnum(fluxEnum),
      m_fluxJacLNeg(fluxJacLNeg), m_fluxJacLPos(fluxJacLPos),
      m_fluxJacRNeg(fluxJacRNeg), m_fluxJacRPos(fluxJacRPos)
  {}

  const flux_jac_type & fluxJacLNeg()  const { return m_fluxJacLNeg; }
  const flux_jac_type & fluxJacLPos()  const { return m_fluxJacLPos; }
  const flux_jac_type & fluxJacRNeg()  const { return m_fluxJacRNeg; }
  const flux_jac_type & fluxJacRPos()  const { return m_fluxJacRPos; }

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
      m_fluxJacLNeg = 1;
      m_fluxJacLPos = 0;
      m_fluxJacRNeg = 1;
      m_fluxJacRPos = 0;
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

public:
  template<class ...Args>
  ComputeDirectionalFluxValuesAndJacobians(InviscidFluxScheme fluxEnum,
					   flux_type & fluxL,
					   flux_type & fluxR,
					   flux_jac_type & fluxJacLNeg,
					   flux_jac_type & fluxJacLPos,
					   flux_jac_type & fluxJacRNeg,
					   flux_jac_type & fluxJacRPos,
					   Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_fluxEnum(fluxEnum),
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
  void operator()(IndexType smPt)
  {
    Parent::operator()(smPt);

    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum){
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      m_fluxL = uMinusHalfNeg;
      m_fluxR = uPlusHalfNeg;
      m_fluxJacLNeg = 1;
      m_fluxJacLPos = 0;
      m_fluxJacRNeg = 1;
      m_fluxJacRPos = 0;
      break;
    }
  }
};

}}
#endif
