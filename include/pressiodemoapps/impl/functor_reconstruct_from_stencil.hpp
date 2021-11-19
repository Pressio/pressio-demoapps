
#ifndef PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STENCIL_HPP_
#define PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STENCIL_HPP_

namespace pressiodemoapps{ namespace impl{

namespace
{
template<class ReconstructedValueType, class StencilDataType>
struct _ReconstructorFromStencilMembers
{
  _ReconstructorFromStencilMembers() = delete;
  _ReconstructorFromStencilMembers(ReconstructionScheme recEn,
				   const StencilDataType & stencilVals,
				   ReconstructedValueType & uMinusHalfNeg,
				   ReconstructedValueType & uMinusHalfPos,
				   ReconstructedValueType & uPlusHalfNeg,
				   ReconstructedValueType & uPlusHalfPos)
    : m_recEn(recEn),
      m_stencilVals(stencilVals),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos)
  {}

  ::pressiodemoapps::ReconstructionScheme m_recEn;
  const StencilDataType & m_stencilVals;
  ReconstructedValueType & m_uMinusHalfNeg;
  ReconstructedValueType & m_uMinusHalfPos;
  ReconstructedValueType & m_uPlusHalfNeg;
  ReconstructedValueType & m_uPlusHalfPos;
};
} // anonym namespace

template<int numDofPerCell, class ReconstructedValueType, class StencilDataType>
class ReconstructorFromStencil;

template<class ReconstructedValueType, class StencilDataType>
class ReconstructorFromStencil<1, ReconstructedValueType, StencilDataType>
{
  using members_t = _ReconstructorFromStencilMembers<ReconstructedValueType, StencilDataType>;
  members_t m_memb;

public:
  template<typename ...Args>
  ReconstructorFromStencil(Args && ... args)
    : m_memb(std::forward<Args>(args)...){}

  const ReconstructionScheme reconstructionScheme() const{
    return m_memb.m_recEn;
  }
  const ReconstructedValueType & reconstructionLeftNeg() const{
    return m_memb.m_uMinusHalfNeg;
  }
  const ReconstructedValueType & reconstructionLeftPos() const{
    return m_memb.m_uMinusHalfPos;
  }
  const ReconstructedValueType & reconstructionRightNeg() const{
    return m_memb.m_uPlusHalfNeg;
  }
  const ReconstructedValueType & reconstructionRightPos() const{
    return m_memb.m_uPlusHalfPos;
  }

  template<class IndexType>
  void operator()(IndexType /*unused*/)
  {
    switch(m_memb.m_recEn){
    case ::pressiodemoapps::ReconstructionScheme::FirstOrder:
      {
	m_memb.m_uMinusHalfNeg = m_memb.m_stencilVals(0);
	m_memb.m_uMinusHalfPos = m_memb.m_stencilVals(1);
	m_memb.m_uPlusHalfNeg  = m_memb.m_stencilVals(1);
	m_memb.m_uPlusHalfPos  = m_memb.m_stencilVals(2);
	break;
      }

    case ::pressiodemoapps::ReconstructionScheme::Weno3:
      {
	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg,
			       m_memb.m_uMinusHalfPos,
			       m_memb.m_uPlusHalfNeg,
			       m_memb.m_uPlusHalfPos,
			       m_memb.m_stencilVals);
	break;
      }

    case ::pressiodemoapps::ReconstructionScheme::Weno5:
      {
	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg,
			       m_memb.m_uMinusHalfPos,
			       m_memb.m_uPlusHalfNeg,
			       m_memb.m_uPlusHalfPos,
			       m_memb.m_stencilVals);
	break;
      }
    }
  }
};


template<class ReconstructedValueType, class StencilDataType>
class ReconstructorFromStencil<3, ReconstructedValueType, StencilDataType>
{
  using members_t = _ReconstructorFromStencilMembers<ReconstructedValueType, StencilDataType>;
  members_t m_memb;

public:
  template<typename ...Args>
  ReconstructorFromStencil(Args && ... args)
    : m_memb(std::forward<Args>(args)...){}

  const ReconstructionScheme reconstructionScheme() const{
    return m_memb.m_recEn;
  }

  const ReconstructedValueType & reconstructionLeftNeg() const{
    return m_memb.m_uMinusHalfNeg;
  }
  const ReconstructedValueType & reconstructionLeftPos() const{
    return m_memb.m_uMinusHalfPos;
  }
  const ReconstructedValueType & reconstructionRightNeg() const{
    return m_memb.m_uPlusHalfNeg;
  }
  const ReconstructedValueType & reconstructionRightPos() const{
    return m_memb.m_uPlusHalfPos;
  }

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    switch(m_memb.m_recEn)
      {
      case ::pressiodemoapps::ReconstructionScheme::FirstOrder:{
	m_memb.m_uMinusHalfNeg(0) = m_memb.m_stencilVals(0);
	m_memb.m_uMinusHalfPos(0) = m_memb.m_stencilVals(3);
	m_memb.m_uPlusHalfNeg(0)  = m_memb.m_stencilVals(3);
	m_memb.m_uPlusHalfPos(0)  = m_memb.m_stencilVals(6);

	m_memb.m_uMinusHalfNeg(1) = m_memb.m_stencilVals(1);
	m_memb.m_uMinusHalfPos(1) = m_memb.m_stencilVals(4);
	m_memb.m_uPlusHalfNeg(1)  = m_memb.m_stencilVals(4);
	m_memb.m_uPlusHalfPos(1)  = m_memb.m_stencilVals(7);

	m_memb.m_uMinusHalfNeg(2) = m_memb.m_stencilVals(2);
	m_memb.m_uMinusHalfPos(2) = m_memb.m_stencilVals(5);
	m_memb.m_uPlusHalfNeg(2)  = m_memb.m_stencilVals(5);
	m_memb.m_uPlusHalfPos(2)  = m_memb.m_stencilVals(8);
	break;
      }

      case ::pressiodemoapps::ReconstructionScheme::Weno3:{
	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(0),
			       m_memb.m_uMinusHalfPos(0),
			       m_memb.m_uPlusHalfNeg(0),
			       m_memb.m_uPlusHalfPos(0),
			       m_memb.m_stencilVals(0),
			       m_memb.m_stencilVals(3),
			       m_memb.m_stencilVals(6),
			       m_memb.m_stencilVals(9),
			       m_memb.m_stencilVals(12));

	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(1),
			       m_memb.m_uMinusHalfPos(1),
			       m_memb.m_uPlusHalfNeg(1),
			       m_memb.m_uPlusHalfPos(1),
			       m_memb.m_stencilVals(1),
			       m_memb.m_stencilVals(4),
			       m_memb.m_stencilVals(7),
			       m_memb.m_stencilVals(10),
			       m_memb.m_stencilVals(13));

	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(2),
			       m_memb.m_uMinusHalfPos(2),
			       m_memb.m_uPlusHalfNeg(2),
			       m_memb.m_uPlusHalfPos(2),
			       m_memb.m_stencilVals(2),
			       m_memb.m_stencilVals(5),
			       m_memb.m_stencilVals(8),
			       m_memb.m_stencilVals(11),
			       m_memb.m_stencilVals(14));
	break;
      }

      case ::pressiodemoapps::ReconstructionScheme::Weno5:{
	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(0),
			       m_memb.m_uMinusHalfPos(0),
			       m_memb.m_uPlusHalfNeg(0),
			       m_memb.m_uPlusHalfPos(0),
			       m_memb.m_stencilVals(0),
			       m_memb.m_stencilVals(3),
			       m_memb.m_stencilVals(6),
			       m_memb.m_stencilVals(9),
			       m_memb.m_stencilVals(12),
			       m_memb.m_stencilVals(15),
			       m_memb.m_stencilVals(18));

	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(1),
			       m_memb.m_uMinusHalfPos(1),
			       m_memb.m_uPlusHalfNeg(1),
			       m_memb.m_uPlusHalfPos(1),
			       m_memb.m_stencilVals(1),
			       m_memb.m_stencilVals(4),
			       m_memb.m_stencilVals(7),
			       m_memb.m_stencilVals(10),
			       m_memb.m_stencilVals(13),
			       m_memb.m_stencilVals(16),
			       m_memb.m_stencilVals(19));

	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(2),
			       m_memb.m_uMinusHalfPos(2),
			       m_memb.m_uPlusHalfNeg(2),
			       m_memb.m_uPlusHalfPos(2),
			       m_memb.m_stencilVals(2),
			       m_memb.m_stencilVals(5),
			       m_memb.m_stencilVals(8),
			       m_memb.m_stencilVals(11),
			       m_memb.m_stencilVals(14),
			       m_memb.m_stencilVals(17),
			       m_memb.m_stencilVals(20));
	break;
      }
      }
  }
};

template<class ReconstructedValueType, class StencilDataType>
class ReconstructorFromStencil<4, ReconstructedValueType, StencilDataType>
{

  using members_t = _ReconstructorFromStencilMembers<ReconstructedValueType, StencilDataType>;
  members_t m_memb;

public:
  template<typename ...Args>
  ReconstructorFromStencil(Args && ... args)
    : m_memb(std::forward<Args>(args)...){}

  const ReconstructionScheme reconstructionScheme() const{
    return m_memb.m_recEn;
  }
  const ReconstructedValueType & reconstructionLeftNeg() const{
    return m_memb.m_uMinusHalfNeg;
  }
  const ReconstructedValueType & reconstructionLeftPos() const{
    return m_memb.m_uMinusHalfPos;
  }
  const ReconstructedValueType & reconstructionRightNeg() const{
    return m_memb.m_uPlusHalfNeg;
  }
  const ReconstructedValueType & reconstructionRightPos() const{
    return m_memb.m_uPlusHalfPos;
  }

  template<class IndexType>
  void operator()(IndexType /*unused*/)
  {
    switch(m_memb.m_recEn)
      {
      case ::pressiodemoapps::ReconstructionScheme::FirstOrder:{
	m_memb.m_uMinusHalfNeg(0) = m_memb.m_stencilVals(0);
	m_memb.m_uMinusHalfPos(0) = m_memb.m_stencilVals(4);
	m_memb.m_uPlusHalfNeg(0)  = m_memb.m_stencilVals(4);
	m_memb.m_uPlusHalfPos(0)  = m_memb.m_stencilVals(8);

	m_memb.m_uMinusHalfNeg(1) = m_memb.m_stencilVals(1);
	m_memb.m_uMinusHalfPos(1) = m_memb.m_stencilVals(5);
	m_memb.m_uPlusHalfNeg(1)  = m_memb.m_stencilVals(5);
	m_memb.m_uPlusHalfPos(1)  = m_memb.m_stencilVals(9);

	m_memb.m_uMinusHalfNeg(2) = m_memb.m_stencilVals(2);
	m_memb.m_uMinusHalfPos(2) = m_memb.m_stencilVals(6);
	m_memb.m_uPlusHalfNeg(2)  = m_memb.m_stencilVals(6);
	m_memb.m_uPlusHalfPos(2)  = m_memb.m_stencilVals(10);

	m_memb.m_uMinusHalfNeg(3) = m_memb.m_stencilVals(3);
	m_memb.m_uMinusHalfPos(3) = m_memb.m_stencilVals(7);
	m_memb.m_uPlusHalfNeg(3)  = m_memb.m_stencilVals(7);
	m_memb.m_uPlusHalfPos(3)  = m_memb.m_stencilVals(11);
	break;
      }

      case ::pressiodemoapps::ReconstructionScheme::Weno3:{
	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(0),
			       m_memb.m_uMinusHalfPos(0),
			       m_memb.m_uPlusHalfNeg(0),
			       m_memb.m_uPlusHalfPos(0),
			       m_memb.m_stencilVals(0),
			       m_memb.m_stencilVals(4),
			       m_memb.m_stencilVals(8),
			       m_memb.m_stencilVals(12),
			       m_memb.m_stencilVals(16));

	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(1),
			       m_memb.m_uMinusHalfPos(1),
			       m_memb.m_uPlusHalfNeg(1),
			       m_memb.m_uPlusHalfPos(1),
			       m_memb.m_stencilVals(1),
			       m_memb.m_stencilVals(5),
			       m_memb.m_stencilVals(9),
			       m_memb.m_stencilVals(13),
			       m_memb.m_stencilVals(17));

	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(2),
			       m_memb.m_uMinusHalfPos(2),
			       m_memb.m_uPlusHalfNeg(2),
			       m_memb.m_uPlusHalfPos(2),
			       m_memb.m_stencilVals(2),
			       m_memb.m_stencilVals(6),
			       m_memb.m_stencilVals(10),
			       m_memb.m_stencilVals(14),
			       m_memb.m_stencilVals(18));

	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(3),
			       m_memb.m_uMinusHalfPos(3),
			       m_memb.m_uPlusHalfNeg(3),
			       m_memb.m_uPlusHalfPos(3),
			       m_memb.m_stencilVals(3),
			       m_memb.m_stencilVals(7),
			       m_memb.m_stencilVals(11),
			       m_memb.m_stencilVals(15),
			       m_memb.m_stencilVals(19));
	break;
      }

      case ::pressiodemoapps::ReconstructionScheme::Weno5:{
	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(0),
			       m_memb.m_uMinusHalfPos(0),
			       m_memb.m_uPlusHalfNeg(0),
			       m_memb.m_uPlusHalfPos(0),
			       m_memb.m_stencilVals(0),
			       m_memb.m_stencilVals(4),
			       m_memb.m_stencilVals(8),
			       m_memb.m_stencilVals(12),
			       m_memb.m_stencilVals(16),
			       m_memb.m_stencilVals(20),
			       m_memb.m_stencilVals(24));

	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(1),
			       m_memb.m_uMinusHalfPos(1),
			       m_memb.m_uPlusHalfNeg(1),
			       m_memb.m_uPlusHalfPos(1),
			       m_memb.m_stencilVals(1),
			       m_memb.m_stencilVals(5),
			       m_memb.m_stencilVals(9),
			       m_memb.m_stencilVals(13),
			       m_memb.m_stencilVals(17),
			       m_memb.m_stencilVals(21),
			       m_memb.m_stencilVals(25));

	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(2),
			       m_memb.m_uMinusHalfPos(2),
			       m_memb.m_uPlusHalfNeg(2),
			       m_memb.m_uPlusHalfPos(2),
			       m_memb.m_stencilVals(2),
			       m_memb.m_stencilVals(6),
			       m_memb.m_stencilVals(10),
			       m_memb.m_stencilVals(14),
			       m_memb.m_stencilVals(18),
			       m_memb.m_stencilVals(22),
			       m_memb.m_stencilVals(26));

	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(3),
			       m_memb.m_uMinusHalfPos(3),
			       m_memb.m_uPlusHalfNeg(3),
			       m_memb.m_uPlusHalfPos(3),
			       m_memb.m_stencilVals(3),
			       m_memb.m_stencilVals(7),
			       m_memb.m_stencilVals(11),
			       m_memb.m_stencilVals(15),
			       m_memb.m_stencilVals(19),
			       m_memb.m_stencilVals(23),
			       m_memb.m_stencilVals(27));
	break;
      }
      }
  }
};


// template<class ReconstructedValueType, class StencilDataType>
// class ReconstructorFromStencil<5, ReconstructedValueType, StencilDataType>
// {
//   using members_t = _ReconstructorFromStencilMembers<ReconstructedValueType, StencilDataType>;
//   members_t m_memb;

// public:
//   template<typename ...Args>
//   ReconstructorFromStencil(Args && ... args)
//     : m_memb(std::forward<Args>(args)...){}

//   const ReconstructionScheme reconstructionScheme() const{
//     return m_memb.m_recEn;
//   }

//   const ReconstructedValueType & reconstructionLeftNeg() const{
//     return m_memb.m_uMinusHalfNeg;
//   }
//   const ReconstructedValueType & reconstructionLeftPos() const{
//     return m_memb.m_uMinusHalfPos;
//   }
//   const ReconstructedValueType & reconstructionRightNeg() const{
//     return m_memb.m_uPlusHalfNeg;
//   }
//   const ReconstructedValueType & reconstructionRightPos() const{
//     return m_memb.m_uPlusHalfPos;
//   }

// template<class IndexType>
// void operator()(IndexType /*unused*/)
//   {
//     switch(m_memb.m_recEn)
//       {
//       case ::pressiodemoapps::ReconstructionScheme::FirstOrder:{
// 	m_memb.m_uMinusHalfNeg(0) = m_memb.m_stencilVals(0);
// 	m_memb.m_uMinusHalfPos(0) = m_memb.m_stencilVals(5);
// 	m_memb.m_uPlusHalfNeg(0)  = m_memb.m_stencilVals(5);
// 	m_memb.m_uPlusHalfPos(0)  = m_memb.m_stencilVals(10);

// 	m_memb.m_uMinusHalfNeg(1) = m_memb.m_stencilVals(1);
// 	m_memb.m_uMinusHalfPos(1) = m_memb.m_stencilVals(6);
// 	m_memb.m_uPlusHalfNeg(1)  = m_memb.m_stencilVals(6);
// 	m_memb.m_uPlusHalfPos(1)  = m_memb.m_stencilVals(11);

// 	m_memb.m_uMinusHalfNeg(2) = m_memb.m_stencilVals(2);
// 	m_memb.m_uMinusHalfPos(2) = m_memb.m_stencilVals(7);
// 	m_memb.m_uPlusHalfNeg(2)  = m_memb.m_stencilVals(7);
// 	m_memb.m_uPlusHalfPos(2)  = m_memb.m_stencilVals(12);

// 	m_memb.m_uMinusHalfNeg(3) = m_memb.m_stencilVals(3);
// 	m_memb.m_uMinusHalfPos(3) = m_memb.m_stencilVals(8);
// 	m_memb.m_uPlusHalfNeg(3)  = m_memb.m_stencilVals(8);
// 	m_memb.m_uPlusHalfPos(3)  = m_memb.m_stencilVals(13);

// 	m_memb.m_uMinusHalfNeg(4) = m_memb.m_stencilVals(4);
// 	m_memb.m_uMinusHalfPos(4) = m_memb.m_stencilVals(9);
// 	m_memb.m_uPlusHalfNeg(4)  = m_memb.m_stencilVals(9);
// 	m_memb.m_uPlusHalfPos(4)  = m_memb.m_stencilVals(14);
// 	break;
//       }

//       case ::pressiodemoapps::ReconstructionScheme::Weno3:{
// 	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(0),
// 			       m_memb.m_uMinusHalfPos(0),
// 			       m_memb.m_uPlusHalfNeg(0),
// 			       m_memb.m_uPlusHalfPos(0),
// 			       m_memb.m_stencilVals(0),
// 			       m_memb.m_stencilVals(5),
// 			       m_memb.m_stencilVals(10),
// 			       m_memb.m_stencilVals(15),
// 			       m_memb.m_stencilVals(20));

// 	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(1),
// 			       m_memb.m_uMinusHalfPos(1),
// 			       m_memb.m_uPlusHalfNeg(1),
// 			       m_memb.m_uPlusHalfPos(1),
// 			       m_memb.m_stencilVals(1),
// 			       m_memb.m_stencilVals(6),
// 			       m_memb.m_stencilVals(11),
// 			       m_memb.m_stencilVals(16),
// 			       m_memb.m_stencilVals(21));

// 	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(2),
// 			       m_memb.m_uMinusHalfPos(2),
// 			       m_memb.m_uPlusHalfNeg(2),
// 			       m_memb.m_uPlusHalfPos(2),
// 			       m_memb.m_stencilVals(2),
// 			       m_memb.m_stencilVals(7),
// 			       m_memb.m_stencilVals(12),
// 			       m_memb.m_stencilVals(17),
// 			       m_memb.m_stencilVals(22));

// 	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(3),
// 			       m_memb.m_uMinusHalfPos(3),
// 			       m_memb.m_uPlusHalfNeg(3),
// 			       m_memb.m_uPlusHalfPos(3),
// 			       m_memb.m_stencilVals(3),
// 			       m_memb.m_stencilVals(8),
// 			       m_memb.m_stencilVals(13),
// 			       m_memb.m_stencilVals(18),
// 			       m_memb.m_stencilVals(23));

// 	pressiodemoapps::weno3(m_memb.m_uMinusHalfNeg(4),
// 			       m_memb.m_uMinusHalfPos(4),
// 			       m_memb.m_uPlusHalfNeg(4),
// 			       m_memb.m_uPlusHalfPos(4),
// 			       m_memb.m_stencilVals(4),
// 			       m_memb.m_stencilVals(9),
// 			       m_memb.m_stencilVals(14),
// 			       m_memb.m_stencilVals(19),
// 			       m_memb.m_stencilVals(24));
// 	break;
//       }

//       case ::pressiodemoapps::ReconstructionScheme::Weno5:{
// 	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(0),
// 			       m_memb.m_uMinusHalfPos(0),
// 			       m_memb.m_uPlusHalfNeg(0),
// 			       m_memb.m_uPlusHalfPos(0),
// 			       m_memb.m_stencilVals(0),
// 			       m_memb.m_stencilVals(5),
// 			       m_memb.m_stencilVals(10),
// 			       m_memb.m_stencilVals(15),
// 			       m_memb.m_stencilVals(20),
// 			       m_memb.m_stencilVals(25),
// 			       m_memb.m_stencilVals(30));

// 	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(1),
// 			       m_memb.m_uMinusHalfPos(1),
// 			       m_memb.m_uPlusHalfNeg(1),
// 			       m_memb.m_uPlusHalfPos(1),
// 			       m_memb.m_stencilVals(1),
// 			       m_memb.m_stencilVals(6),
// 			       m_memb.m_stencilVals(11),
// 			       m_memb.m_stencilVals(16),
// 			       m_memb.m_stencilVals(21),
// 			       m_memb.m_stencilVals(26),
// 			       m_memb.m_stencilVals(31));

// 	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(2),
// 			       m_memb.m_uMinusHalfPos(2),
// 			       m_memb.m_uPlusHalfNeg(2),
// 			       m_memb.m_uPlusHalfPos(2),
// 			       m_memb.m_stencilVals(2),
// 			       m_memb.m_stencilVals(7),
// 			       m_memb.m_stencilVals(12),
// 			       m_memb.m_stencilVals(17),
// 			       m_memb.m_stencilVals(22),
// 			       m_memb.m_stencilVals(27),
// 			       m_memb.m_stencilVals(32));

// 	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(3),
// 			       m_memb.m_uMinusHalfPos(3),
// 			       m_memb.m_uPlusHalfNeg(3),
// 			       m_memb.m_uPlusHalfPos(3),
// 			       m_memb.m_stencilVals(3),
// 			       m_memb.m_stencilVals(8),
// 			       m_memb.m_stencilVals(13),
// 			       m_memb.m_stencilVals(18),
// 			       m_memb.m_stencilVals(23),
// 			       m_memb.m_stencilVals(28),
// 			       m_memb.m_stencilVals(33));

// 	pressiodemoapps::weno5(m_memb.m_uMinusHalfNeg(4),
// 			       m_memb.m_uMinusHalfPos(4),
// 			       m_memb.m_uPlusHalfNeg(4),
// 			       m_memb.m_uPlusHalfPos(4),
// 			       m_memb.m_stencilVals(4),
// 			       m_memb.m_stencilVals(13),
// 			       m_memb.m_stencilVals(18),
// 			       m_memb.m_stencilVals(23),
// 			       m_memb.m_stencilVals(28),
// 			       m_memb.m_stencilVals(33),
// 			       m_memb.m_stencilVals(38));
// 	break;
//       }
//       }
//   }
// };

}}
#endif
