
#ifndef PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_HPP_
#define PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_HPP_

namespace pressiodemoapps{ namespace impl{

template<int numDofPerCell, class edge_rec_t, class stencil_values>
class Reconstructor
{

public:
  Reconstructor() = delete;
  Reconstructor(::pressiodemoapps::reconstructionEnum recEn,
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

  template<int _ndpc = numDofPerCell>
  typename std::enable_if<_ndpc == 1>::type
  operator()()
  {
    switch(m_recEn){
    case pressiodemoapps::reconstructionEnum::firstOrder:
      {
	m_uMinusHalfNeg = m_stencilVals(0);
	m_uMinusHalfPos = m_stencilVals(1);
	m_uPlusHalfNeg  = m_stencilVals(1);
	m_uPlusHalfPos  = m_stencilVals(2);
	break;
      }

    case pressiodemoapps::reconstructionEnum::fifthOrderWeno:
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

  template<int _ndpc = numDofPerCell>
  typename std::enable_if< (_ndpc > 1) >::type
  operator()()
  {
    switch(m_recEn){
    case pressiodemoapps::reconstructionEnum::firstOrder:
      {
	for (int i=0; i<_ndpc; ++i){
	  m_uMinusHalfNeg(i) = m_stencilVals(0+i);
	  m_uMinusHalfPos(i) = m_stencilVals(_ndpc+i);
	  m_uPlusHalfNeg(i)  = m_stencilVals(_ndpc+i);
	  m_uPlusHalfPos(i)  = m_stencilVals(2*_ndpc+i);
	}
	break;
      }

    case pressiodemoapps::reconstructionEnum::fifthOrderWeno:
      {
	using scalar_type = typename edge_rec_t::value_type;
	std::array<scalar_type,7> mys;
	for (int dof=0; dof<_ndpc; ++dof)
	{
	  int start = dof;
	  for (int k=0; k<7; ++k){
	    mys[k] = m_stencilVals(start + _ndpc*k);
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
  ::pressiodemoapps::reconstructionEnum m_recEn;
  const stencil_values & m_stencilVals;
  edge_rec_t & m_uMinusHalfNeg;
  edge_rec_t & m_uMinusHalfPos;
  edge_rec_t & m_uPlusHalfNeg;
  edge_rec_t & m_uPlusHalfPos;

};

}}
#endif
