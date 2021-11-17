
#ifndef PRESSIODEMOAPPS_EE_FLUXE_FUNCTOR_HPP_
#define PRESSIODEMOAPPS_EE_FLUXE_FUNCTOR_HPP_

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<int ndpc, class scalar_type, class edge_rec_type, class flux_type>
struct FluxFunctorWithoutGradients
{
  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  const edge_rec_type & m_uMinusHalfNeg;
  const edge_rec_type & m_uMinusHalfPos;
  const edge_rec_type & m_uPlusHalfNeg;
  const edge_rec_type & m_uPlusHalfPos;
  flux_type & m_fluxL;
  flux_type & m_fluxR;

  FluxFunctorWithoutGradients(InviscidFluxScheme fluxEnum,
			      scalar_type gamma,
			      const edge_rec_type & uMinusHalfNeg,
			      const edge_rec_type & uMinusHalfPos,
			      const edge_rec_type & uPlusHalfNeg,
			      const edge_rec_type & uPlusHalfPos,
			      flux_type & fluxL,
			      flux_type & fluxR)
    : m_fluxEnum(fluxEnum),
      m_gamma(gamma),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos),
      m_fluxL(fluxL),
      m_fluxR(fluxR)
  {}

  template<int _ndpc = ndpc>
  std::enable_if_t<_ndpc == 3> operator()()
  {
    namespace pda = ::pressiodemoapps;
    switch(m_fluxEnum)
    {
    case pda::InviscidFluxScheme::Rusanov:
      ee::impl::eeRusanovFluxThreeDof(m_fluxL, m_uMinusHalfNeg,
				      m_uMinusHalfPos, m_gamma);
      ee::impl::eeRusanovFluxThreeDof(m_fluxR, m_uPlusHalfNeg,
				      m_uPlusHalfPos,  m_gamma);
      break;
    }
  }

  template<int _ndpc = ndpc>
  std::enable_if_t<_ndpc == 4> operator()()
  {
    namespace pda = ::pressiodemoapps;
    switch(m_fluxEnum)
    {
    case pda::InviscidFluxScheme::Rusanov:
      ee::impl::eeRusanovFluxFourDof(m_fluxL, m_uMinusHalfNeg,
				     m_uMinusHalfPos, m_gamma);
      ee::impl::eeRusanovFluxFourDof(m_fluxR, m_uPlusHalfNeg,
				     m_uPlusHalfPos,  m_gamma);
      break;
    }
  }
};

template<int ndpc, class scalar_type, class edge_rec_type, class flux_type, class flux_jac_type>
struct FluxFunctorWithGradients
{

  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  const edge_rec_type & m_uMinusHalfNeg;
  const edge_rec_type & m_uMinusHalfPos;
  const edge_rec_type & m_uPlusHalfNeg;
  const edge_rec_type & m_uPlusHalfPos;
  flux_type & m_fluxL;
  flux_type & m_fluxR;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;

  FluxFunctorWithGradients(InviscidFluxScheme fluxEnum,
			   scalar_type gamma,
			   const edge_rec_type & uMinusHalfNeg,
			   const edge_rec_type & uMinusHalfPos,
			   const edge_rec_type & uPlusHalfNeg,
			   const edge_rec_type & uPlusHalfPos,
			   flux_type & fluxL,
			   flux_type & fluxR,
			   flux_jac_type & fluxJacLNeg,
			   flux_jac_type & fluxJacLPos,
			   flux_jac_type & fluxJacRNeg,
			   flux_jac_type & fluxJacRPos)
    : m_fluxEnum(fluxEnum),
      m_gamma(gamma),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos),
      m_fluxL(fluxL),
      m_fluxR(fluxR),
      m_fluxJacLNeg(fluxJacLNeg),
      m_fluxJacLPos(fluxJacLPos),
      m_fluxJacRNeg(fluxJacRNeg),
      m_fluxJacRPos(fluxJacRPos)
  {}

  template<int _ndpc = ndpc>
  std::enable_if_t<_ndpc == 3> operator()()
  {
    namespace pda = ::pressiodemoapps;

    switch(m_fluxEnum)
    {
    case pda::InviscidFluxScheme::Rusanov:
      ee::impl::eeRusanovFluxThreeDof(m_fluxL, m_uMinusHalfNeg,
				      m_uMinusHalfPos, m_gamma);
      ee::impl::eeRusanovFluxThreeDof(m_fluxR, m_uPlusHalfNeg,
				      m_uPlusHalfPos,  m_gamma);

      ee::impl::eeRusanovFluxJacobianThreeDof(m_fluxJacLNeg, m_fluxJacLPos,
					      m_uMinusHalfNeg, m_uMinusHalfPos,
					      m_gamma);
      ee::impl::eeRusanovFluxJacobianThreeDof(m_fluxJacRNeg, m_fluxJacRPos,
					      m_uPlusHalfNeg, m_uPlusHalfPos,
					      m_gamma);
      break;
    }
  }
};

template<int ndpc, class scalar_type, class edge_rec_type, class flux_type, class flux_jac_type>
struct FluxFunctorGradUsingDifferentRecontructions
{

  InviscidFluxScheme m_fluxEnum;
  scalar_type m_gamma;
  const edge_rec_type & m_uMinusHalfNeg;
  const edge_rec_type & m_uMinusHalfPos;
  const edge_rec_type & m_uPlusHalfNeg;
  const edge_rec_type & m_uPlusHalfPos;
  const edge_rec_type & m_uMinusHalfNegForJ;
  const edge_rec_type & m_uMinusHalfPosForJ;
  const edge_rec_type & m_uPlusHalfNegForJ;
  const edge_rec_type & m_uPlusHalfPosForJ;
  flux_type & m_fluxL;
  flux_type & m_fluxR;
  flux_jac_type & m_fluxJacLNeg;
  flux_jac_type & m_fluxJacLPos;
  flux_jac_type & m_fluxJacRNeg;
  flux_jac_type & m_fluxJacRPos;

  FluxFunctorGradUsingDifferentRecontructions(InviscidFluxScheme fluxEnum,
					      scalar_type gamma,
					      const edge_rec_type & uMinusHalfNeg,
					      const edge_rec_type & uMinusHalfPos,
					      const edge_rec_type & uPlusHalfNeg,
					      const edge_rec_type & uPlusHalfPos,
					      const edge_rec_type & uMinusHalfNegForJ,
					      const edge_rec_type & uMinusHalfPosForJ,
					      const edge_rec_type & uPlusHalfNegForJ,
					      const edge_rec_type & uPlusHalfPosForJ,
					      flux_type & fluxL,
					      flux_type & fluxR,
					      flux_jac_type & fluxJacLNeg,
					      flux_jac_type & fluxJacLPos,
					      flux_jac_type & fluxJacRNeg,
					      flux_jac_type & fluxJacRPos)
  : m_fluxEnum(fluxEnum),
    m_gamma(gamma),
    m_uMinusHalfNeg(uMinusHalfNeg),
    m_uMinusHalfPos(uMinusHalfPos),
    m_uPlusHalfNeg(uPlusHalfNeg),
    m_uPlusHalfPos(uPlusHalfPos),
    m_uMinusHalfNegForJ(uMinusHalfNegForJ),
    m_uMinusHalfPosForJ(uMinusHalfPosForJ),
    m_uPlusHalfNegForJ(uPlusHalfNegForJ),
    m_uPlusHalfPosForJ(uPlusHalfPosForJ),
    m_fluxL(fluxL),
    m_fluxR(fluxR),
    m_fluxJacLNeg(fluxJacLNeg),
    m_fluxJacLPos(fluxJacLPos),
    m_fluxJacRNeg(fluxJacRNeg),
    m_fluxJacRPos(fluxJacRPos)
  {}

  template<int _ndpc = ndpc>
  std::enable_if_t<_ndpc == 3> operator()()
  {
    namespace pda = ::pressiodemoapps;

    switch(m_fluxEnum)
    {
    case pda::InviscidFluxScheme::Rusanov:
      ee::impl::eeRusanovFluxThreeDof(m_fluxL, m_uMinusHalfNeg,
				      m_uMinusHalfPos, m_gamma);
      ee::impl::eeRusanovFluxThreeDof(m_fluxR, m_uPlusHalfNeg,
				      m_uPlusHalfPos,  m_gamma);

      ee::impl::eeRusanovFluxJacobianThreeDof(m_fluxJacLNeg, m_fluxJacLPos,
					      m_uMinusHalfNegForJ, m_uMinusHalfPosForJ,
					      m_gamma);
      ee::impl::eeRusanovFluxJacobianThreeDof(m_fluxJacRNeg, m_fluxJacRPos,
					      m_uPlusHalfNegForJ, m_uPlusHalfPosForJ,
					      m_gamma);
      break;
    }
  }
};

}}}
#endif
