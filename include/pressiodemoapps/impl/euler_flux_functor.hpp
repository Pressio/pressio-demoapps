
#ifndef PRESSIODEMOAPPS_EE_FLUXE_FUNCTOR_HPP_
#define PRESSIODEMOAPPS_EE_FLUXE_FUNCTOR_HPP_

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<class Parent, int ndpc, class scalar_type, class flux_type>
struct FluxValues : Parent
{
private:
  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  flux_type & m_fluxL;
  flux_type & m_fluxR;

public:
  template<class ...Args>
  FluxValues(InviscidFluxScheme fluxEnum,
	     scalar_type gamma,
	     flux_type & fluxL,
	     flux_type & fluxR,
	     Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_fluxEnum(fluxEnum), m_gamma(gamma),
      m_fluxL(fluxL), m_fluxR(fluxR)
  {}

  const flux_type & fluxLeft()  const { return m_fluxL; }
  const flux_type & fluxRight() const { return m_fluxR; }

  template<class IndexType, int _ndpc = ndpc>
  std::enable_if_t<_ndpc == 3> operator()(IndexType smPt)
  {
    Parent::operator()(smPt);

    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum)
    {
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      ee::impl::eeRusanovFluxThreeDof(m_fluxL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
      ee::impl::eeRusanovFluxThreeDof(m_fluxR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);
      break;
    }
  }
};

template<class Parent, int ndpc, class scalar_type, class flux_jac_type>
struct FluxJacobians : Parent
{
private:
  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;

public:
  template<class ...Args>
  FluxJacobians(InviscidFluxScheme fluxEnum,
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

  template<class IndexType, int _ndpc = ndpc>
  std::enable_if_t<_ndpc == 3> operator()(IndexType smPt)
  {
    Parent::operator()(smPt);

    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum)
    {
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      ee::impl::eeRusanovFluxJacobianThreeDof(m_fluxJacLNeg, m_fluxJacLPos,
					      uMinusHalfNeg, uMinusHalfPos,
					      m_gamma);
      ee::impl::eeRusanovFluxJacobianThreeDof(m_fluxJacRNeg, m_fluxJacRPos,
					      uPlusHalfNeg, uPlusHalfPos,
					      m_gamma);
      break;
    }
  }
};


template<class Parent, int ndpc, class scalar_type, class flux_type, class flux_jac_type>
struct FluxValuesAndJacobians : Parent
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
  FluxValuesAndJacobians(InviscidFluxScheme fluxEnum,
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

  template<class IndexType, int _ndpc = ndpc>
  std::enable_if_t<_ndpc == 3> operator()(IndexType smPt)
  {
    Parent::operator()(smPt);

    const auto & uMinusHalfNeg = Parent::reconstructionLeftNeg();
    const auto & uMinusHalfPos = Parent::reconstructionLeftPos();
    const auto & uPlusHalfNeg  = Parent::reconstructionRightNeg();
    const auto & uPlusHalfPos  = Parent::reconstructionRightPos();

    switch(m_fluxEnum)
    {
    case ::pressiodemoapps::InviscidFluxScheme::Rusanov:
      ee::impl::eeRusanovFluxThreeDof(m_fluxL, uMinusHalfNeg,
				      uMinusHalfPos, m_gamma);
      ee::impl::eeRusanovFluxThreeDof(m_fluxR, uPlusHalfNeg,
				      uPlusHalfPos,  m_gamma);
      ee::impl::eeRusanovFluxJacobianThreeDof(m_fluxJacLNeg, m_fluxJacLPos,
					      uMinusHalfNeg, uMinusHalfPos,
					      m_gamma);
      ee::impl::eeRusanovFluxJacobianThreeDof(m_fluxJacRNeg, m_fluxJacRPos,
					      uPlusHalfNeg, uPlusHalfPos,
					      m_gamma);
      break;
    }
  }
};

}}}
#endif
