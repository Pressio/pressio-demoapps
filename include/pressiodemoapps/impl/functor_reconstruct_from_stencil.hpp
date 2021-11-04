
#ifndef PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STENCIL_HPP_
#define PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STENCIL_HPP_

namespace pressiodemoapps{ namespace impl{

template<class ValueType, class StencilValuesType>
class ReconstructorFromStencilOneDofPerCell
{

public:
  ReconstructorFromStencilOneDofPerCell() = delete;
  ReconstructorFromStencilOneDofPerCell(::pressiodemoapps::InviscidFluxReconstruction recEn,
					const StencilValuesType & stencilVals,
					ValueType & uMinusHalfNeg,
					ValueType & uMinusHalfPos,
					ValueType & uPlusHalfNeg,
					ValueType & uPlusHalfPos)
    : m_recEn(recEn),
      m_stencilVals(stencilVals),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos)
  {}

  void operator()()
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

private:
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  const StencilValuesType & m_stencilVals;
  ValueType & m_uMinusHalfNeg;
  ValueType & m_uMinusHalfPos;
  ValueType & m_uPlusHalfNeg;
  ValueType & m_uPlusHalfPos;

};


template<class ValueType, class StencilValuesType>
class ReconstructorFromStencilThreeDofPerCell
{

public:
  ReconstructorFromStencilThreeDofPerCell() = delete;
  ReconstructorFromStencilThreeDofPerCell(::pressiodemoapps::InviscidFluxReconstruction recEn,
					const StencilValuesType & stencilVals,
					ValueType & uMinusHalfNeg,
					ValueType & uMinusHalfPos,
					ValueType & uPlusHalfNeg,
					ValueType & uPlusHalfPos)
    : m_recEn(recEn),
      m_stencilVals(stencilVals),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos)
  {}

  void operator()()
  {
    switch(m_recEn)
      {
      case ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder:
	{

	  m_uMinusHalfNeg(0) = m_stencilVals(0);
	  m_uMinusHalfPos(0) = m_stencilVals(3);
	  m_uPlusHalfNeg(0)  = m_stencilVals(3);
	  m_uPlusHalfPos(0)  = m_stencilVals(6);

	  m_uMinusHalfNeg(1) = m_stencilVals(1);
	  m_uMinusHalfPos(1) = m_stencilVals(4);
	  m_uPlusHalfNeg(1)  = m_stencilVals(4);
	  m_uPlusHalfPos(1)  = m_stencilVals(7);

	  m_uMinusHalfNeg(2) = m_stencilVals(2);
	  m_uMinusHalfPos(2) = m_stencilVals(5);
	  m_uPlusHalfNeg(2)  = m_stencilVals(5);
	  m_uPlusHalfPos(2)  = m_stencilVals(8);

	  break;
	}

      case ::pressiodemoapps::InviscidFluxReconstruction::Weno3:
	{

	  pressiodemoapps::weno3(m_uMinusHalfNeg(0),
				 m_uMinusHalfPos(0),
				 m_uPlusHalfNeg(0),
				 m_uPlusHalfPos(0),
				 m_stencilVals(0),
				 m_stencilVals(3),
				 m_stencilVals(6),
				 m_stencilVals(9),
				 m_stencilVals(12));

	  pressiodemoapps::weno3(m_uMinusHalfNeg(1),
				 m_uMinusHalfPos(1),
				 m_uPlusHalfNeg(1),
				 m_uPlusHalfPos(1),
				 m_stencilVals(1),
				 m_stencilVals(4),
				 m_stencilVals(7),
				 m_stencilVals(10),
				 m_stencilVals(13));

	  pressiodemoapps::weno3(m_uMinusHalfNeg(2),
				 m_uMinusHalfPos(2),
				 m_uPlusHalfNeg(2),
				 m_uPlusHalfPos(2),
				 m_stencilVals(2),
				 m_stencilVals(5),
				 m_stencilVals(8),
				 m_stencilVals(11),
				 m_stencilVals(14));

	  break;
	}

      case ::pressiodemoapps::InviscidFluxReconstruction::Weno5:
	{

	  pressiodemoapps::weno5(m_uMinusHalfNeg(0),
				 m_uMinusHalfPos(0),
				 m_uPlusHalfNeg(0),
				 m_uPlusHalfPos(0),
				 m_stencilVals(0),
				 m_stencilVals(3),
				 m_stencilVals(6),
				 m_stencilVals(9),
				 m_stencilVals(12),
				 m_stencilVals(15),
				 m_stencilVals(18));

	  pressiodemoapps::weno5(m_uMinusHalfNeg(1),
				 m_uMinusHalfPos(1),
				 m_uPlusHalfNeg(1),
				 m_uPlusHalfPos(1),
				 m_stencilVals(1),
				 m_stencilVals(4),
				 m_stencilVals(7),
				 m_stencilVals(10),
				 m_stencilVals(13),
				 m_stencilVals(16),
				 m_stencilVals(19));

	  pressiodemoapps::weno5(m_uMinusHalfNeg(2),
				 m_uMinusHalfPos(2),
				 m_uPlusHalfNeg(2),
				 m_uPlusHalfPos(2),
				 m_stencilVals(2),
				 m_stencilVals(5),
				 m_stencilVals(8),
				 m_stencilVals(11),
				 m_stencilVals(14),
				 m_stencilVals(17),
				 m_stencilVals(20));

	  break;
	}
      }
  }

private:
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  const StencilValuesType & m_stencilVals;
  ValueType & m_uMinusHalfNeg;
  ValueType & m_uMinusHalfPos;
  ValueType & m_uPlusHalfNeg;
  ValueType & m_uPlusHalfPos;
};


template<class ValueType, class StencilValuesType>
class ReconstructorFromStencilFourDofPerCell
{

public:
  ReconstructorFromStencilFourDofPerCell() = delete;
  ReconstructorFromStencilFourDofPerCell(::pressiodemoapps::InviscidFluxReconstruction recEn,
					const StencilValuesType & stencilVals,
					ValueType & uMinusHalfNeg,
					ValueType & uMinusHalfPos,
					ValueType & uPlusHalfNeg,
					ValueType & uPlusHalfPos)
    : m_recEn(recEn),
      m_stencilVals(stencilVals),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos)
  {}

  void operator()()
  {

    switch(m_recEn)
      {
      case ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder:
	{

	  m_uMinusHalfNeg(0) = m_stencilVals(0);
	  m_uMinusHalfPos(0) = m_stencilVals(4);
	  m_uPlusHalfNeg(0)  = m_stencilVals(4);
	  m_uPlusHalfPos(0)  = m_stencilVals(8);

	  m_uMinusHalfNeg(1) = m_stencilVals(1);
	  m_uMinusHalfPos(1) = m_stencilVals(5);
	  m_uPlusHalfNeg(1)  = m_stencilVals(5);
	  m_uPlusHalfPos(1)  = m_stencilVals(9);

	  m_uMinusHalfNeg(2) = m_stencilVals(2);
	  m_uMinusHalfPos(2) = m_stencilVals(6);
	  m_uPlusHalfNeg(2)  = m_stencilVals(6);
	  m_uPlusHalfPos(2)  = m_stencilVals(10);

	  m_uMinusHalfNeg(3) = m_stencilVals(3);
	  m_uMinusHalfPos(3) = m_stencilVals(7);
	  m_uPlusHalfNeg(3)  = m_stencilVals(7);
	  m_uPlusHalfPos(3)  = m_stencilVals(11);

	  break;
	}

      case ::pressiodemoapps::InviscidFluxReconstruction::Weno3:
	{

	  pressiodemoapps::weno3(m_uMinusHalfNeg(0),
				 m_uMinusHalfPos(0),
				 m_uPlusHalfNeg(0),
				 m_uPlusHalfPos(0),
				 m_stencilVals(0),
				 m_stencilVals(4),
				 m_stencilVals(8),
				 m_stencilVals(12),
				 m_stencilVals(16));

	  pressiodemoapps::weno3(m_uMinusHalfNeg(1),
				 m_uMinusHalfPos(1),
				 m_uPlusHalfNeg(1),
				 m_uPlusHalfPos(1),
				 m_stencilVals(1),
				 m_stencilVals(5),
				 m_stencilVals(9),
				 m_stencilVals(13),
				 m_stencilVals(17));

	  pressiodemoapps::weno3(m_uMinusHalfNeg(2),
				 m_uMinusHalfPos(2),
				 m_uPlusHalfNeg(2),
				 m_uPlusHalfPos(2),
				 m_stencilVals(2),
				 m_stencilVals(6),
				 m_stencilVals(10),
				 m_stencilVals(14),
				 m_stencilVals(18));

	  pressiodemoapps::weno3(m_uMinusHalfNeg(3),
				 m_uMinusHalfPos(3),
				 m_uPlusHalfNeg(3),
				 m_uPlusHalfPos(3),
				 m_stencilVals(3),
				 m_stencilVals(7),
				 m_stencilVals(11),
				 m_stencilVals(15),
				 m_stencilVals(19));

	  break;
	}

      case ::pressiodemoapps::InviscidFluxReconstruction::Weno5:
	{

	  pressiodemoapps::weno5(m_uMinusHalfNeg(0),
				 m_uMinusHalfPos(0),
				 m_uPlusHalfNeg(0),
				 m_uPlusHalfPos(0),
				 m_stencilVals(0),
				 m_stencilVals(4),
				 m_stencilVals(8),
				 m_stencilVals(12),
				 m_stencilVals(16),
				 m_stencilVals(20),
				 m_stencilVals(24));

	  pressiodemoapps::weno5(m_uMinusHalfNeg(1),
				 m_uMinusHalfPos(1),
				 m_uPlusHalfNeg(1),
				 m_uPlusHalfPos(1),
				 m_stencilVals(1),
				 m_stencilVals(5),
				 m_stencilVals(9),
				 m_stencilVals(13),
				 m_stencilVals(17),
				 m_stencilVals(21),
				 m_stencilVals(25));

	  pressiodemoapps::weno5(m_uMinusHalfNeg(2),
				 m_uMinusHalfPos(2),
				 m_uPlusHalfNeg(2),
				 m_uPlusHalfPos(2),
				 m_stencilVals(2),
				 m_stencilVals(6),
				 m_stencilVals(10),
				 m_stencilVals(14),
				 m_stencilVals(18),
				 m_stencilVals(22),
				 m_stencilVals(26));

	  pressiodemoapps::weno5(m_uMinusHalfNeg(3),
				 m_uMinusHalfPos(3),
				 m_uPlusHalfNeg(3),
				 m_uPlusHalfPos(3),
				 m_stencilVals(3),
				 m_stencilVals(7),
				 m_stencilVals(11),
				 m_stencilVals(15),
				 m_stencilVals(19),
				 m_stencilVals(23),
				 m_stencilVals(27));

	  break;
	}
      }
  }

private:
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  const StencilValuesType & m_stencilVals;
  ValueType & m_uMinusHalfNeg;
  ValueType & m_uMinusHalfPos;
  ValueType & m_uPlusHalfNeg;
  ValueType & m_uPlusHalfPos;
};


template<class ValueType, class StencilValuesType>
class ReconstructorFromStencilFiveDofPerCell
{

public:
  ReconstructorFromStencilFiveDofPerCell() = delete;
  ReconstructorFromStencilFiveDofPerCell(::pressiodemoapps::InviscidFluxReconstruction recEn,
					const StencilValuesType & stencilVals,
					ValueType & uMinusHalfNeg,
					ValueType & uMinusHalfPos,
					ValueType & uPlusHalfNeg,
					ValueType & uPlusHalfPos)
    : m_recEn(recEn),
      m_stencilVals(stencilVals),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos)
  {}

  void operator()()
  {

    switch(m_recEn)
      {
      case ::pressiodemoapps::InviscidFluxReconstruction::FirstOrder:
	{

	  m_uMinusHalfNeg(0) = m_stencilVals(0);
	  m_uMinusHalfPos(0) = m_stencilVals(5);
	  m_uPlusHalfNeg(0)  = m_stencilVals(5);
	  m_uPlusHalfPos(0)  = m_stencilVals(10);

	  m_uMinusHalfNeg(1) = m_stencilVals(1);
	  m_uMinusHalfPos(1) = m_stencilVals(6);
	  m_uPlusHalfNeg(1)  = m_stencilVals(6);
	  m_uPlusHalfPos(1)  = m_stencilVals(11);

	  m_uMinusHalfNeg(2) = m_stencilVals(2);
	  m_uMinusHalfPos(2) = m_stencilVals(7);
	  m_uPlusHalfNeg(2)  = m_stencilVals(7);
	  m_uPlusHalfPos(2)  = m_stencilVals(12);

	  m_uMinusHalfNeg(3) = m_stencilVals(3);
	  m_uMinusHalfPos(3) = m_stencilVals(8);
	  m_uPlusHalfNeg(3)  = m_stencilVals(8);
	  m_uPlusHalfPos(3)  = m_stencilVals(13);

	  m_uMinusHalfNeg(4) = m_stencilVals(4);
	  m_uMinusHalfPos(4) = m_stencilVals(9);
	  m_uPlusHalfNeg(4)  = m_stencilVals(9);
	  m_uPlusHalfPos(4)  = m_stencilVals(14);

	  break;
	}

      case ::pressiodemoapps::InviscidFluxReconstruction::Weno3:
	{

	  pressiodemoapps::weno3(m_uMinusHalfNeg(0),
				 m_uMinusHalfPos(0),
				 m_uPlusHalfNeg(0),
				 m_uPlusHalfPos(0),
				 m_stencilVals(0),
				 m_stencilVals(5),
				 m_stencilVals(10),
				 m_stencilVals(15),
				 m_stencilVals(20));

	  pressiodemoapps::weno3(m_uMinusHalfNeg(1),
				 m_uMinusHalfPos(1),
				 m_uPlusHalfNeg(1),
				 m_uPlusHalfPos(1),
				 m_stencilVals(1),
				 m_stencilVals(6),
				 m_stencilVals(11),
				 m_stencilVals(16),
				 m_stencilVals(21));

	  pressiodemoapps::weno3(m_uMinusHalfNeg(2),
				 m_uMinusHalfPos(2),
				 m_uPlusHalfNeg(2),
				 m_uPlusHalfPos(2),
				 m_stencilVals(2),
				 m_stencilVals(7),
				 m_stencilVals(12),
				 m_stencilVals(17),
				 m_stencilVals(22));

	  pressiodemoapps::weno3(m_uMinusHalfNeg(3),
				 m_uMinusHalfPos(3),
				 m_uPlusHalfNeg(3),
				 m_uPlusHalfPos(3),
				 m_stencilVals(3),
				 m_stencilVals(8),
				 m_stencilVals(13),
				 m_stencilVals(18),
				 m_stencilVals(23));

	  pressiodemoapps::weno3(m_uMinusHalfNeg(4),
				 m_uMinusHalfPos(4),
				 m_uPlusHalfNeg(4),
				 m_uPlusHalfPos(4),
				 m_stencilVals(4),
				 m_stencilVals(9),
				 m_stencilVals(14),
				 m_stencilVals(19),
				 m_stencilVals(24));

	  break;
	}

      case ::pressiodemoapps::InviscidFluxReconstruction::Weno5:
	{

	  pressiodemoapps::weno5(m_uMinusHalfNeg(0),
				 m_uMinusHalfPos(0),
				 m_uPlusHalfNeg(0),
				 m_uPlusHalfPos(0),
				 m_stencilVals(0),
				 m_stencilVals(5),
				 m_stencilVals(10),
				 m_stencilVals(15),
				 m_stencilVals(20),
				 m_stencilVals(25),
				 m_stencilVals(30));

	  pressiodemoapps::weno5(m_uMinusHalfNeg(1),
				 m_uMinusHalfPos(1),
				 m_uPlusHalfNeg(1),
				 m_uPlusHalfPos(1),
				 m_stencilVals(1),
				 m_stencilVals(6),
				 m_stencilVals(11),
				 m_stencilVals(16),
				 m_stencilVals(21),
				 m_stencilVals(26),
				 m_stencilVals(31));

	  pressiodemoapps::weno5(m_uMinusHalfNeg(2),
				 m_uMinusHalfPos(2),
				 m_uPlusHalfNeg(2),
				 m_uPlusHalfPos(2),
				 m_stencilVals(2),
				 m_stencilVals(7),
				 m_stencilVals(12),
				 m_stencilVals(17),
				 m_stencilVals(22),
				 m_stencilVals(27),
				 m_stencilVals(32));

	  pressiodemoapps::weno5(m_uMinusHalfNeg(3),
				 m_uMinusHalfPos(3),
				 m_uPlusHalfNeg(3),
				 m_uPlusHalfPos(3),
				 m_stencilVals(3),
				 m_stencilVals(8),
				 m_stencilVals(13),
				 m_stencilVals(18),
				 m_stencilVals(23),
				 m_stencilVals(28),
				 m_stencilVals(33));

	  pressiodemoapps::weno5(m_uMinusHalfNeg(4),
				 m_uMinusHalfPos(4),
				 m_uPlusHalfNeg(4),
				 m_uPlusHalfPos(4),
				 m_stencilVals(4),
				 m_stencilVals(13),
				 m_stencilVals(18),
				 m_stencilVals(23),
				 m_stencilVals(28),
				 m_stencilVals(33),
				 m_stencilVals(38));

	  break;
	}
      }
  }

private:
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;
  const StencilValuesType & m_stencilVals;
  ValueType & m_uMinusHalfNeg;
  ValueType & m_uMinusHalfPos;
  ValueType & m_uPlusHalfNeg;
  ValueType & m_uPlusHalfPos;
};

}}
#endif
