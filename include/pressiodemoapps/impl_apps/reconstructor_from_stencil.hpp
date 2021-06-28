
#ifndef PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STENCIL_HPP_
#define PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STENCIL_HPP_

namespace pressiodemoapps{ namespace impl{

template<class edge_rec_t, class stencil_values>
class ReconstructorFromStencil
{

public:
  ReconstructorFromStencil() = delete;
  ReconstructorFromStencil(::pressiodemoapps::InviscidFluxReconstruction recEn,
			   const stencil_values & stencilVals,
			   edge_rec_t & uMinusHalfNeg,
			   edge_rec_t & uMinusHalfPos,
			   edge_rec_t & uPlusHalfNeg,
			   edge_rec_t & uPlusHalfPos)
    : m_recEn(recEn),
      m_stencilVals(stencilVals),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos)
  {}

  // one dof per cell
  template<int ndpc>
  typename std::enable_if<ndpc == 1>::type
  operator()()
  {
    switch(m_recEn){
    case ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder:
      {
	m_uMinusHalfNeg = m_stencilVals(0);
	m_uMinusHalfPos = m_stencilVals(1);
	m_uPlusHalfNeg  = m_stencilVals(1);
	m_uPlusHalfPos  = m_stencilVals(2);
	break;
      }

    case ::pressiodemoapps::InviscidFluxReconstruction::Weno3:
      {
	pressiodemoapps::weno3(m_uMinusHalfNeg,
			       m_uMinusHalfPos,
			       m_uPlusHalfNeg,
			       m_uPlusHalfPos,
			       m_stencilVals);
	break;
      }

    case ::pressiodemoapps::InviscidFluxReconstruction::Weno5:
      {
	pressiodemoapps::weno5(m_uMinusHalfNeg,
			       m_uMinusHalfPos,
			       m_uPlusHalfNeg,
			       m_uPlusHalfPos,
			       m_stencilVals);
	break;
      }
    }
  }

  // multiple dofs per cell
  template<int ndpc>
  typename std::enable_if< (ndpc > 1) >::type
  operator()()
  {

    switch(m_recEn)
      {
      case ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder:
	{
	  for (int i=0; i<ndpc; ++i){
	    m_uMinusHalfNeg(i) = m_stencilVals(0+i);
	    m_uMinusHalfPos(i) = m_stencilVals(ndpc+i);
	    m_uPlusHalfNeg(i)  = m_stencilVals(ndpc+i);
	    m_uPlusHalfPos(i)  = m_stencilVals(2*ndpc+i);
	  }
	  break;
	}

      case ::pressiodemoapps::InviscidFluxReconstruction::Weno3:
	{
	  using scalar_type = typename edge_rec_t::value_type;
	  std::array<scalar_type,5> mys;
	  for (int dof=0; dof<ndpc; ++dof)
	    {
	      int start = dof;
	      for (int k=0; k<5; ++k){
		mys[k] = m_stencilVals(start + ndpc*k);
	      }
	      pressiodemoapps::weno3(m_uMinusHalfNeg(dof),
				     m_uMinusHalfPos(dof),
				     m_uPlusHalfNeg(dof),
				     m_uPlusHalfPos(dof),
				     mys);
	    }
	  break;
	}

      case ::pressiodemoapps::InviscidFluxReconstruction::Weno5:
	{
	  using scalar_type = typename edge_rec_t::value_type;
	  std::array<scalar_type,7> mys;
	  for (int dof=0; dof<ndpc; ++dof)
	    {
	      int start = dof;
	      for (int k=0; k<7; ++k){
		mys[k] = m_stencilVals(start + ndpc*k);
	      }
	      pressiodemoapps::weno5(m_uMinusHalfNeg(dof),
				     m_uMinusHalfPos(dof),
				     m_uPlusHalfNeg(dof),
				     m_uPlusHalfPos(dof),
				     mys);
	    }
	  break;
	}
      }
  }

private:
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  const stencil_values & m_stencilVals;
  edge_rec_t & m_uMinusHalfNeg;
  edge_rec_t & m_uMinusHalfPos;
  edge_rec_t & m_uPlusHalfNeg;
  edge_rec_t & m_uPlusHalfPos;

};

}}
#endif
