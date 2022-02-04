
#ifndef PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STATE_HPP_
#define PRESSIODEMOAPPS_EDGE_RECONSTRUCTOR_FROM_STATE_HPP_

#include <iomanip>

namespace pressiodemoapps{ namespace impl{

namespace{

template<class ...Args> struct _ReconstructorMembers;

// no gradients needed, only values
template<
  class ReconstructedValueType,
  class DataType,
  class MeshType
  >
struct _ReconstructorMembers<ReconstructedValueType, DataType, MeshType>
{
  _ReconstructorMembers() = delete;

  _ReconstructorMembers(const int axis,
			ReconstructionScheme recEn,
			const DataType & functionValues,
			const MeshType & meshObj,
			ReconstructedValueType & uMinusHalfNeg,
			ReconstructedValueType & uMinusHalfPos,
			ReconstructedValueType & uPlusHalfNeg,
			ReconstructedValueType & uPlusHalfPos)
    : m_axis(axis),
      m_recEn(recEn),
      m_functionValues(functionValues),
      m_meshObj(meshObj),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos)
  {}

  _ReconstructorMembers(ReconstructionScheme recEn,
			const DataType & functionValues,
			const MeshType & meshObj,
			ReconstructedValueType & uMinusHalfNeg,
			ReconstructedValueType & uMinusHalfPos,
			ReconstructedValueType & uPlusHalfNeg,
			ReconstructedValueType & uPlusHalfPos)
    : _ReconstructorMembers(1, /* default axis is 1 */
			    recEn, functionValues, meshObj,
			    uMinusHalfNeg, uMinusHalfPos,
			    uPlusHalfNeg, uPlusHalfPos)
  {}

  ReconstructionScheme reconstructionScheme() const{
    return m_recEn;
  }

  const ReconstructedValueType & reconstructionLeftNeg() const{
    return m_uMinusHalfNeg;
  }

  const ReconstructedValueType & reconstructionLeftPos() const{
    return m_uMinusHalfPos;
  }

  const ReconstructedValueType & reconstructionRightNeg() const{
    return m_uPlusHalfNeg;
  }

  const ReconstructedValueType & reconstructionRightPos() const{
    return m_uPlusHalfPos;
  }

  const int m_axis = -1;
  const ReconstructionScheme m_recEn;
  const DataType & m_functionValues;
  const MeshType & m_meshObj;
  ReconstructedValueType & m_uMinusHalfNeg;
  ReconstructedValueType & m_uMinusHalfPos;
  ReconstructedValueType & m_uPlusHalfNeg;
  ReconstructedValueType & m_uPlusHalfPos;
};

// reconstruction gradients are needed
template<
  class ReconstructedValueType,
  class DataType,
  class MeshType,
  class ReconstructedGradType
  >
struct _ReconstructorMembers<
  ReconstructedValueType, DataType, MeshType, ReconstructedGradType
  >
{
  _ReconstructorMembers() = delete;

  _ReconstructorMembers(const int axis,
			ReconstructionScheme recEn,
			const DataType & functionValues,
			const MeshType & meshObj,
			ReconstructedValueType & uMinusHalfNeg,
			ReconstructedValueType & uMinusHalfPos,
			ReconstructedValueType & uPlusHalfNeg,
			ReconstructedValueType & uPlusHalfPos,
			ReconstructedGradType & gradLNeg,
			ReconstructedGradType & gradLPos,
			ReconstructedGradType & gradRNeg,
			ReconstructedGradType & gradRPos)
    : m_axis(axis),
      m_recEn(recEn),
      m_functionValues(functionValues),
      m_meshObj(meshObj),
      m_uMinusHalfNeg(uMinusHalfNeg),
      m_uMinusHalfPos(uMinusHalfPos),
      m_uPlusHalfNeg(uPlusHalfNeg),
      m_uPlusHalfPos(uPlusHalfPos),
      m_gradLNeg(gradLNeg),
      m_gradLPos(gradLPos),
      m_gradRNeg(gradRNeg),
      m_gradRPos(gradRPos)
  {}


  _ReconstructorMembers(ReconstructionScheme recEn,
			const DataType & functionValues,
			const MeshType & meshObj,
			ReconstructedValueType & uMinusHalfNeg,
			ReconstructedValueType & uMinusHalfPos,
			ReconstructedValueType & uPlusHalfNeg,
			ReconstructedValueType & uPlusHalfPos,
			ReconstructedGradType & gradLNeg,
			ReconstructedGradType & gradLPos,
			ReconstructedGradType & gradRNeg,
			ReconstructedGradType & gradRPos)
    : _ReconstructorMembers(1, /* default axis is 1 */
			    recEn, functionValues, meshObj,
			    uMinusHalfNeg, uMinusHalfPos,
			    uPlusHalfNeg, uPlusHalfPos,
			    gradLNeg, gradLPos,
			    gradRNeg, gradRPos)
  {}

  ReconstructionScheme reconstructionScheme() const{
    return m_recEn;
  }

  const ReconstructedValueType & reconstructionLeftNeg() const{
    return m_uMinusHalfNeg;
  }
  const ReconstructedValueType & reconstructionLeftPos() const{
    return m_uMinusHalfPos;
  }
  const ReconstructedValueType & reconstructionRightNeg() const{
    return m_uPlusHalfNeg;
  }
  const ReconstructedValueType & reconstructionRightPos() const{
    return m_uPlusHalfPos;
  }

  const ReconstructedGradType & reconstructionGradLeftNeg() const{
    return m_gradLNeg;
  }
  const ReconstructedGradType & reconstructionGradLeftPos() const{
    return m_gradLPos;
  }
  const ReconstructedGradType & reconstructionGradRightNeg() const{
    return m_gradRNeg;
  }
  const ReconstructedGradType & reconstructionGradRightPos() const{
    return m_gradRPos;
  }

  const int m_axis = -1;
  const ReconstructionScheme m_recEn;
  const DataType & m_functionValues;
  const MeshType & m_meshObj;
  ReconstructedValueType & m_uMinusHalfNeg;
  ReconstructedValueType & m_uMinusHalfPos;
  ReconstructedValueType & m_uPlusHalfNeg;
  ReconstructedValueType & m_uPlusHalfPos;
  ReconstructedGradType & m_gradLNeg ;
  ReconstructedGradType & m_gradLPos ;
  ReconstructedGradType & m_gradRNeg ;
  ReconstructedGradType & m_gradRPos ;
};
} // anonym namespace


template<int dim, int numDofPerCell, class ...Args>
class ReconstructorForDiscreteFunction;

// -----------------------------------------------------------
// 1d, multiple dofsPerCell, only values, no gradients
// -----------------------------------------------------------
template<
  int ndpc,
  class MeshType,
  class DataType,
  class ReconstructedValueType
  >
class ReconstructorForDiscreteFunction<
  1, ndpc, MeshType, DataType, ReconstructedValueType
  > : public _ReconstructorMembers<ReconstructedValueType, DataType, MeshType>
{
  using members_t = _ReconstructorMembers<ReconstructedValueType, DataType, MeshType>;

public:
  template<typename ...Args>
  ReconstructorForDiscreteFunction(Args && ... args)
    : members_t(std::forward<Args>(args)...)
  {
    // this is for 1d, axis should always be 1
    assert(members_t::m_axis == 1);
  }

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = members_t::m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0)*ndpc;

    /* For 1d, graph is:
       ptID left0 right0 [left1 right1 left2 right2]
        0     1     2    [  3     4      5     6   ]

      where [] are optionally present, depending on the stencil chosen.
    */

    const auto l0i = graph(smPt, 1)*ndpc;
    const auto r0i = graph(smPt, 2)*ndpc;

    switch(members_t::m_recEn){
    case ReconstructionScheme::FirstOrder:{
      for (int i=0; i<ndpc; ++i){
	members_t::m_uMinusHalfNeg(i) = members_t::m_functionValues(l0i+i);
	members_t::m_uMinusHalfPos(i) = members_t::m_functionValues(uIndex+i);
	members_t::m_uPlusHalfNeg(i)  = members_t::m_functionValues(uIndex+i);
	members_t::m_uPlusHalfPos(i)  = members_t::m_functionValues(r0i+i);
      }
      break;
    }

    case ReconstructionScheme::Weno3:{
      const auto l1i  = graph(smPt, 3)*ndpc;
      const auto r1i  = graph(smPt, 4)*ndpc;
      ::pressiodemoapps::weno3(members_t::m_uMinusHalfNeg,
			       members_t::m_uMinusHalfPos,
			       members_t::m_uPlusHalfNeg,
			       members_t::m_uPlusHalfPos,
			       members_t::m_functionValues,
			       l1i, l0i, r0i, r1i,
			       uIndex,
			       ndpc);
      break;
    }

    case ReconstructionScheme::Weno5:{
      const auto l1i = graph(smPt, 3)*ndpc;
      const auto r1i = graph(smPt, 4)*ndpc;
      const auto l2i = graph(smPt, 5)*ndpc;
      const auto r2i = graph(smPt, 6)*ndpc;
      ::pressiodemoapps::weno5(members_t::m_uMinusHalfNeg,
			       members_t::m_uMinusHalfPos,
			       members_t::m_uPlusHalfNeg,
			       members_t::m_uPlusHalfPos,
			       members_t::m_functionValues,
			       l2i, l1i, l0i,
			       r0i, r1i, r2i,
			       uIndex,
			       ndpc);
      break;
    }
    }//end switch
  }
};

// -----------------------------------------------------------
// 1d, multiple dofsPerCell, WITH reconstruction gradients
// -----------------------------------------------------------
template<
  int ndpc,
  class MeshType,
  class DataType,
  class ReconstructedValueType,
  class ReconstructedGradType
  >
class ReconstructorForDiscreteFunction<
  1, ndpc, MeshType, DataType, ReconstructedValueType, ReconstructedGradType
  > : public _ReconstructorMembers<ReconstructedValueType, DataType, MeshType, ReconstructedGradType>
{
  using members_t = _ReconstructorMembers<ReconstructedValueType, DataType,
					  MeshType, ReconstructedGradType>;

public:
  template<typename ...Args>
  ReconstructorForDiscreteFunction(Args && ... args)
    : members_t(std::forward<Args>(args)...)
  {
    // this is for 1d, axis should always be 1
    assert(members_t::m_axis == 1);
  }

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = members_t::m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0)*ndpc;
    const auto l0i = graph(smPt, 1)*ndpc;
    const auto r0i = graph(smPt, 2)*ndpc;

    switch(members_t::m_recEn){
    case ReconstructionScheme::FirstOrder:{
      for (int i=0; i<ndpc; ++i){
	members_t::m_uMinusHalfNeg(i) = members_t::m_functionValues(l0i+i);
	members_t::m_uMinusHalfPos(i) = members_t::m_functionValues(uIndex+i);
	members_t::m_uPlusHalfNeg(i)  = members_t::m_functionValues(uIndex+i);
	members_t::m_uPlusHalfPos(i)  = members_t::m_functionValues(r0i+i);
	members_t::m_gradLNeg(i, 0) = 1;
	members_t::m_gradLNeg(i, 1) = 0;
	members_t::m_gradLPos(i, 0) = 0;
	members_t::m_gradLPos(i, 1) = 1;
	members_t::m_gradRNeg(i, 0) = 1;
	members_t::m_gradRNeg(i, 1) = 0;
	members_t::m_gradRPos(i, 0) = 0;
	members_t::m_gradRPos(i, 1) = 1;
      }
      break;
    }

    case ReconstructionScheme::Weno3:{
      const auto l1i  = graph(smPt, 3)*ndpc;
      const auto r1i  = graph(smPt, 4)*ndpc;
      ::pressiodemoapps::weno3(members_t::m_uMinusHalfNeg,
			       members_t::m_uMinusHalfPos,
			       members_t::m_uPlusHalfNeg,
			       members_t::m_uPlusHalfPos,
			       members_t::m_gradLNeg,
			       members_t::m_gradLPos,
			       members_t::m_gradRNeg,
			       members_t::m_gradRPos,
			       members_t::m_functionValues,
			       l1i, l0i, r0i, r1i,
			       uIndex,
			       ndpc);
      break;
    }

    case ReconstructionScheme::Weno5:{
      const auto l1i = graph(smPt, 3)*ndpc;
      const auto r1i = graph(smPt, 4)*ndpc;
      const auto l2i = graph(smPt, 5)*ndpc;
      const auto r2i = graph(smPt, 6)*ndpc;
      ::pressiodemoapps::weno5(members_t::m_uMinusHalfNeg,
			       members_t::m_uMinusHalfPos,
			       members_t::m_uPlusHalfNeg,
			       members_t::m_uPlusHalfPos,
			       members_t::m_gradLNeg,
			       members_t::m_gradLPos,
			       members_t::m_gradRNeg,
			       members_t::m_gradRPos,
			       members_t::m_functionValues,
			       l2i, l1i, l0i,
			       r0i, r1i, r2i,
			       uIndex,
			       ndpc);
      break;
    }
    }// end switch
  }
};

// -----------------------------------------------------------
// 1d, 1 dof, no reconstruction gradients
// -----------------------------------------------------------
template<class MeshType, class DataType, class ReconstructedValueType>
class ReconstructorForDiscreteFunction<
  1, 1, MeshType, DataType, ReconstructedValueType
  > : public _ReconstructorMembers<ReconstructedValueType, DataType, MeshType>
{
  using members_t = _ReconstructorMembers<ReconstructedValueType, DataType, MeshType>;

public:
  template<typename ...Args>
  ReconstructorForDiscreteFunction(Args && ... args)
    : members_t(std::forward<Args>(args)...)
  {
    // this is for 1d, axis should always be 1
    assert(members_t::m_axis == 1);
  }

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = members_t::m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0);
    const auto l0i = graph(smPt, 1);
    const auto r0i = graph(smPt, 2);

    switch(members_t::m_recEn){
    case ReconstructionScheme::FirstOrder:{
      members_t::m_uMinusHalfNeg = members_t::m_functionValues(l0i);
      members_t::m_uMinusHalfPos = members_t::m_functionValues(uIndex);
      members_t::m_uPlusHalfNeg  = members_t::m_functionValues(uIndex);
      members_t::m_uPlusHalfPos  = members_t::m_functionValues(r0i);
      break;
    }

    case ReconstructionScheme::Weno3:{
      const auto l1i = graph(smPt, 3);
      const auto r1i = graph(smPt, 4);
      ::pressiodemoapps::weno3(members_t::m_uMinusHalfNeg,
			       members_t::m_uMinusHalfPos,
			       members_t::m_uPlusHalfNeg,
			       members_t::m_uPlusHalfPos,
			       members_t::m_functionValues(l1i),
			       members_t::m_functionValues(l0i),
			       members_t::m_functionValues(uIndex),
			       members_t::m_functionValues(r0i),
			       members_t::m_functionValues(r1i));
      break;
    }

    case ReconstructionScheme::Weno5:{
      const auto l1i = graph(smPt, 3);
      const auto r1i = graph(smPt, 4);
      const auto l2i = graph(smPt, 5);
      const auto r2i = graph(smPt, 6);
      ::pressiodemoapps::weno5(members_t::m_uMinusHalfNeg,
			       members_t::m_uMinusHalfPos,
			       members_t::m_uPlusHalfNeg,
			       members_t::m_uPlusHalfPos,
			       members_t::m_functionValues(l2i),
			       members_t::m_functionValues(l1i),
			       members_t::m_functionValues(l0i),
			       members_t::m_functionValues(uIndex),
			       members_t::m_functionValues(r0i),
			       members_t::m_functionValues(r1i),
			       members_t::m_functionValues(r2i));
      break;
    }
    }
  }
};

// -----------------------------------------------------------
// 1d, 1 dofsPerCell, WITH reconstruction gradients
// -----------------------------------------------------------
template<
  class MeshType,
  class DataType,
  class ReconstructedValueType,
  class ReconstructedGradType
  >
class ReconstructorForDiscreteFunction<
  1, 1, MeshType, DataType, ReconstructedValueType, ReconstructedGradType
  > : public _ReconstructorMembers<ReconstructedValueType, DataType, MeshType, ReconstructedGradType>
{
  using members_t = _ReconstructorMembers<ReconstructedValueType, DataType,
					  MeshType, ReconstructedGradType>;

public:
  template<typename ...Args>
  ReconstructorForDiscreteFunction(Args && ... args)
    : members_t(std::forward<Args>(args)...)
  {
    // this is for 1d, axis should always be 1
    assert(members_t::m_axis == 1);
  }

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = members_t::m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0);
    const auto l0i  = graph(smPt, 1);
    const auto r0i  = graph(smPt, 2);

    switch(members_t::m_recEn){
    case ReconstructionScheme::FirstOrder:
      {
	members_t::m_uMinusHalfNeg = members_t::m_functionValues(l0i);
	members_t::m_uMinusHalfPos = members_t::m_functionValues(uIndex);
	members_t::m_uPlusHalfNeg  = members_t::m_functionValues(uIndex);
	members_t::m_uPlusHalfPos  = members_t::m_functionValues(r0i);
	members_t::m_gradLNeg[0] = 1;
	members_t::m_gradLNeg[1] = 0;
	members_t::m_gradLPos[0] = 0;
	members_t::m_gradLPos[1] = 1;
	members_t::m_gradRNeg[0] = 1;
	members_t::m_gradRNeg[1] = 0;
	members_t::m_gradRPos[0] = 0;
	members_t::m_gradRPos[1] = 1;
	break;
      }

    case ReconstructionScheme::Weno3:
      {
	const auto l1i = graph(smPt, 3);
	const auto r1i = graph(smPt, 4);
	::pressiodemoapps::weno3(members_t::m_uMinusHalfNeg,
				 members_t::m_uMinusHalfPos,
				 members_t::m_uPlusHalfNeg,
				 members_t::m_uPlusHalfPos,
				 members_t::m_gradLNeg,
				 members_t::m_gradLPos,
				 members_t::m_gradRNeg,
				 members_t::m_gradRPos,
				 members_t::m_functionValues(l1i),
				 members_t::m_functionValues(l0i),
				 members_t::m_functionValues(uIndex),
				 members_t::m_functionValues(r0i),
				 members_t::m_functionValues(r1i));
	break;
      }

    case ReconstructionScheme::Weno5:
      {
	const auto l1i  = graph(smPt, 3);
	const auto r1i  = graph(smPt, 4);
	const auto l2i  = graph(smPt, 5);
	const auto r2i  = graph(smPt, 6);
	::pressiodemoapps::weno5(members_t::m_uMinusHalfNeg,
				 members_t::m_uMinusHalfPos,
				 members_t::m_uPlusHalfNeg,
				 members_t::m_uPlusHalfPos,
				 members_t::m_gradLNeg,
				 members_t::m_gradLPos,
				 members_t::m_gradRNeg,
				 members_t::m_gradRPos,
				 members_t::m_functionValues(l2i),
				 members_t::m_functionValues(l1i),
				 members_t::m_functionValues(l0i),
				 members_t::m_functionValues(uIndex),
				 members_t::m_functionValues(r0i),
				 members_t::m_functionValues(r1i),
				 members_t::m_functionValues(r2i));
	break;
      }
    }// end switch
  }
};


// -----------------------------------------------------------
// 2d, multiple dofsPerCell, no reconstruction gradients
// -----------------------------------------------------------
template<
  int ndpc,
  class MeshType,
  class DataType,
  class ReconstructedValueType
  >
class ReconstructorForDiscreteFunction<
  2, ndpc, MeshType, DataType, ReconstructedValueType
  > : public _ReconstructorMembers<ReconstructedValueType, DataType, MeshType>
{
  using members_t = _ReconstructorMembers<ReconstructedValueType, DataType, MeshType>;

public:
  template<typename ...Args>
  ReconstructorForDiscreteFunction(Args && ... args)
    : members_t(std::forward<Args>(args)...){}

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = members_t::m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0)*ndpc;
    const auto axis = members_t::m_axis;

    /*
      For 2d, graph is:
      ptID left0 front0 right0 back0 [left1 front1 right1 back1 left2 front2 right2 back2]
       0     1     2      3     4    [  5     6      7      8    9      10     11    12  ]

      where [] are optionally present, depending on the stencil chosen.
    */

    switch(members_t::m_recEn){
    case ReconstructionScheme::FirstOrder:{
      const auto l0i = (axis == 1) ? graph(smPt, 1)*ndpc : graph(smPt, 4)*ndpc;
      const auto r0i = (axis == 1) ? graph(smPt, 3)*ndpc : graph(smPt, 2)*ndpc;
      for (int i=0; i<ndpc; ++i){
	members_t::m_uMinusHalfNeg(i) = members_t::m_functionValues(l0i+i);
	members_t::m_uMinusHalfPos(i) = members_t::m_functionValues(uIndex+i);
	members_t::m_uPlusHalfNeg(i)  = members_t::m_functionValues(uIndex+i);
	members_t::m_uPlusHalfPos(i)  = members_t::m_functionValues(r0i+i);
      }
      break;
    }

    case ReconstructionScheme::Weno3:{
      const auto l0i = (axis == 1) ? graph(smPt, 1)*ndpc  : graph(smPt, 4)*ndpc;
      const auto r0i = (axis == 1) ? graph(smPt, 3)*ndpc  : graph(smPt, 2)*ndpc;
      const auto l1i = (axis == 1) ? graph(smPt, 5)*ndpc  : graph(smPt, 8)*ndpc;
      const auto r1i = (axis == 1) ? graph(smPt, 7)*ndpc  : graph(smPt, 6)*ndpc;
      ::pressiodemoapps::weno3(members_t::m_uMinusHalfNeg,
			       members_t::m_uMinusHalfPos,
			       members_t::m_uPlusHalfNeg,
			       members_t::m_uPlusHalfPos,
			       members_t::m_functionValues,
			       l1i, l0i, r0i, r1i,
			       uIndex,
			       ndpc);
      break;
    }

    case ReconstructionScheme::Weno5:{
      const auto l0i = (axis == 1) ? graph(smPt, 1)*ndpc  : graph(smPt, 4)*ndpc;
      const auto r0i = (axis == 1) ? graph(smPt, 3)*ndpc  : graph(smPt, 2)*ndpc;
      const auto l1i = (axis == 1) ? graph(smPt, 5)*ndpc  : graph(smPt, 8)*ndpc;
      const auto r1i = (axis == 1) ? graph(smPt, 7)*ndpc  : graph(smPt, 6)*ndpc;
      const auto l2i = (axis == 1) ? graph(smPt, 9)*ndpc  : graph(smPt, 12)*ndpc;
      const auto r2i = (axis == 1) ? graph(smPt, 11)*ndpc : graph(smPt, 10)*ndpc;
      ::pressiodemoapps::weno5(members_t::m_uMinusHalfNeg,
			       members_t::m_uMinusHalfPos,
			       members_t::m_uPlusHalfNeg,
			       members_t::m_uPlusHalfPos,
			       members_t::m_functionValues,
			       l2i, l1i, l0i,
			       r0i, r1i, r2i,
			       uIndex,
			       ndpc);
      break;
    }
    }// end switch
  }
};

// -----------------------------------------------------------
// 2d, multiple dofsPerCell, WITH reconstruction gradients
// -----------------------------------------------------------
template<
  int ndpc,
  class MeshType,
  class DataType,
  class ReconstructedValueType,
  class ReconstructedGradType
  >
class ReconstructorForDiscreteFunction<
  2, ndpc, MeshType, DataType, ReconstructedValueType, ReconstructedGradType
  > : public _ReconstructorMembers<ReconstructedValueType, DataType, MeshType, ReconstructedGradType>
{
  using members_t = _ReconstructorMembers<ReconstructedValueType, DataType,
					  MeshType, ReconstructedGradType>;

public:
  template<typename ...Args>
  ReconstructorForDiscreteFunction(Args && ... args)
    : members_t(std::forward<Args>(args)...){}

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = members_t::m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0)*ndpc;
    const auto axis = members_t::m_axis;

    switch(members_t::m_recEn){
    case ReconstructionScheme::FirstOrder:{
      const auto l0i = (axis == 1) ? graph(smPt, 1)*ndpc : graph(smPt, 4)*ndpc;
      const auto r0i = (axis == 1) ? graph(smPt, 3)*ndpc : graph(smPt, 2)*ndpc;
      for (int i=0; i<ndpc; ++i){
	members_t::m_uMinusHalfNeg(i) = members_t::m_functionValues(l0i+i);
	members_t::m_uMinusHalfPos(i) = members_t::m_functionValues(uIndex+i);
	members_t::m_uPlusHalfNeg(i)  = members_t::m_functionValues(uIndex+i);
	members_t::m_uPlusHalfPos(i)  = members_t::m_functionValues(r0i+i);
	members_t::m_gradLNeg(i, 0) = 1;
	members_t::m_gradLNeg(i, 1) = 0;
	members_t::m_gradLPos(i, 0) = 0;
	members_t::m_gradLPos(i, 1) = 1;
	members_t::m_gradRNeg(i, 0) = 1;
	members_t::m_gradRNeg(i, 1) = 0;
	members_t::m_gradRPos(i, 0) = 0;
	members_t::m_gradRPos(i, 1) = 1;
      }

      break;
    }

    case ReconstructionScheme::Weno3:{
      const auto l0i = (axis == 1) ? graph(smPt, 1)*ndpc  : graph(smPt, 4)*ndpc;
      const auto r0i = (axis == 1) ? graph(smPt, 3)*ndpc  : graph(smPt, 2)*ndpc;
      const auto l1i = (axis == 1) ? graph(smPt, 5)*ndpc  : graph(smPt, 8)*ndpc;
      const auto r1i = (axis == 1) ? graph(smPt, 7)*ndpc  : graph(smPt, 6)*ndpc;
      ::pressiodemoapps::weno3(members_t::m_uMinusHalfNeg,
			       members_t::m_uMinusHalfPos,
			       members_t::m_uPlusHalfNeg,
			       members_t::m_uPlusHalfPos,
			       members_t::m_gradLNeg,
			       members_t::m_gradLPos,
			       members_t::m_gradRNeg,
			       members_t::m_gradRPos,
			       members_t::m_functionValues,
			       l1i, l0i, r0i, r1i,
			       uIndex,
			       ndpc);
      break;
    }

    case ReconstructionScheme::Weno5:{
      const auto l0i = (axis == 1) ? graph(smPt, 1)*ndpc  : graph(smPt, 4)*ndpc;
      const auto r0i = (axis == 1) ? graph(smPt, 3)*ndpc  : graph(smPt, 2)*ndpc;
      const auto l1i = (axis == 1) ? graph(smPt, 5)*ndpc  : graph(smPt, 8)*ndpc;
      const auto r1i = (axis == 1) ? graph(smPt, 7)*ndpc  : graph(smPt, 6)*ndpc;
      const auto l2i = (axis == 1) ? graph(smPt, 9)*ndpc  : graph(smPt, 12)*ndpc;
      const auto r2i = (axis == 1) ? graph(smPt, 11)*ndpc : graph(smPt, 10)*ndpc;
      ::pressiodemoapps::weno5(members_t::m_uMinusHalfNeg,
			       members_t::m_uMinusHalfPos,
			       members_t::m_uPlusHalfNeg,
			       members_t::m_uPlusHalfPos,
			       members_t::m_gradLNeg,
			       members_t::m_gradLPos,
			       members_t::m_gradRNeg,
			       members_t::m_gradRPos,
			       members_t::m_functionValues,
			       l2i, l1i, l0i,
			       r0i, r1i, r2i,
			       uIndex,
			       ndpc);
      break;
    }
    }//end switch
  }
};

// -----------------------------------------------------------
// 3d, multiple dofsPerCell, no reconstruction gradients
// -----------------------------------------------------------
template<
  int ndpc,
  class MeshType,
  class DataType,
  class ReconstructedValueType
  >
class ReconstructorForDiscreteFunction<
  3, ndpc, MeshType, DataType, ReconstructedValueType
  > : public _ReconstructorMembers<ReconstructedValueType, DataType, MeshType>
{
  using members_t = _ReconstructorMembers<ReconstructedValueType, DataType, MeshType>;

public:
  template<typename ...Args>
  ReconstructorForDiscreteFunction(Args && ... args)
    : members_t(std::forward<Args>(args)...){}

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = members_t::m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0)*ndpc;
    const auto axis = members_t::m_axis;

    /*
      For 3d, graph is:
      ptID l0 f0 r0 ba0 bot0 top0 [l1 f1 r1 ba1 bot1 top1 l2 f2 r2 ba2 bot2 top2]
       0    1  2  3  4    5   6   [ 7  8  9  10  11   12  13 14 15  16  17  18  ]

      where [] are optionally present, depending on the stencil chosen.
    */

    switch(members_t::m_recEn)
      {
      case ReconstructionScheme::FirstOrder:{
	auto l0i = (axis == 1) ? graph(smPt, 1)*ndpc : (axis==2) ? graph(smPt, 4)*ndpc : graph(smPt, 5)*ndpc;
	auto r0i = (axis == 1) ? graph(smPt, 3)*ndpc : (axis==2) ? graph(smPt, 2)*ndpc : graph(smPt, 6)*ndpc;
	for (int i=0; i<ndpc; ++i){
	  members_t::m_uMinusHalfNeg(i) = members_t::m_functionValues(l0i+i);
	  members_t::m_uMinusHalfPos(i) = members_t::m_functionValues(uIndex+i);
	  members_t::m_uPlusHalfNeg(i)  = members_t::m_functionValues(uIndex+i);
	  members_t::m_uPlusHalfPos(i)  = members_t::m_functionValues(r0i+i);
	}
	break;
      }

      case ReconstructionScheme::Weno3:{
	auto l0i = (axis == 1) ? graph(smPt, 1)*ndpc : (axis==2) ? graph(smPt, 4)*ndpc : graph(smPt, 5)*ndpc;
	auto r0i = (axis == 1) ? graph(smPt, 3)*ndpc : (axis==2) ? graph(smPt, 2)*ndpc : graph(smPt, 6)*ndpc;
	auto l1i = (axis == 1) ? graph(smPt, 7)*ndpc : (axis==2) ? graph(smPt, 10)*ndpc : graph(smPt, 11)*ndpc;
	auto r1i = (axis == 1) ? graph(smPt, 9)*ndpc : (axis==2) ? graph(smPt, 8)*ndpc : graph(smPt, 12)*ndpc;
	::pressiodemoapps::weno3(members_t::m_uMinusHalfNeg,
				 members_t::m_uMinusHalfPos,
				 members_t::m_uPlusHalfNeg,
				 members_t::m_uPlusHalfPos,
				 members_t::m_functionValues,
				 l1i, l0i, r0i, r1i,
				 uIndex,
				 ndpc);
	break;
      }

      case ReconstructionScheme::Weno5:
	{
	  throw std::runtime_error("missing impl");
	  break;
	}
      }
  }
};

// -----------------------------------------------------------
// 3d, multiple dofsPerCell, WITH reconstruction gradients
// -----------------------------------------------------------
template<
  int ndpc,
  class MeshType,
  class DataType,
  class ReconstructedValueType,
  class ReconstructedGradType
  >
class ReconstructorForDiscreteFunction<
  3, ndpc, MeshType, DataType, ReconstructedValueType, ReconstructedGradType
  > : public _ReconstructorMembers<ReconstructedValueType, DataType, MeshType, ReconstructedGradType>
{
  using members_t = _ReconstructorMembers<
    ReconstructedValueType, DataType, MeshType, ReconstructedGradType>;

public:
  template<typename ...Args>
  ReconstructorForDiscreteFunction(Args && ... args)
    : members_t(std::forward<Args>(args)...){}

  template<class IndexType>
  void operator()(IndexType smPt)
  {
    const auto & graph  = members_t::m_meshObj.graph();
    const auto & uIndex = graph(smPt, 0)*ndpc;
    const auto axis = members_t::m_axis;

    switch(members_t::m_recEn)
      {
      case ReconstructionScheme::FirstOrder:{
	auto l0i = (axis == 1) ? graph(smPt, 1)*ndpc : (axis==2) ? graph(smPt, 4)*ndpc : graph(smPt, 5)*ndpc;
	auto r0i = (axis == 1) ? graph(smPt, 3)*ndpc : (axis==2) ? graph(smPt, 2)*ndpc : graph(smPt, 6)*ndpc;
	for (int i=0; i<ndpc; ++i){
	  members_t::m_uMinusHalfNeg(i) = members_t::m_functionValues(l0i+i);
	  members_t::m_uMinusHalfPos(i) = members_t::m_functionValues(uIndex+i);
	  members_t::m_uPlusHalfNeg(i)  = members_t::m_functionValues(uIndex+i);
	  members_t::m_uPlusHalfPos(i)  = members_t::m_functionValues(r0i+i);
	}
	break;
      }

      case ReconstructionScheme::Weno3:{
	auto l0i = (axis == 1) ? graph(smPt, 1)*ndpc : (axis==2) ? graph(smPt, 4)*ndpc : graph(smPt, 5)*ndpc;
	auto r0i = (axis == 1) ? graph(smPt, 3)*ndpc : (axis==2) ? graph(smPt, 2)*ndpc : graph(smPt, 6)*ndpc;
	auto l1i = (axis == 1) ? graph(smPt, 7)*ndpc : (axis==2) ? graph(smPt, 10)*ndpc : graph(smPt, 11)*ndpc;
	auto r1i = (axis == 1) ? graph(smPt, 9)*ndpc : (axis==2) ? graph(smPt, 8)*ndpc : graph(smPt, 12)*ndpc;
	::pressiodemoapps::weno3(members_t::m_uMinusHalfNeg,
				 members_t::m_uMinusHalfPos,
				 members_t::m_uPlusHalfNeg,
				 members_t::m_uPlusHalfPos,
				 members_t::m_gradLNeg,
				 members_t::m_gradLPos,
				 members_t::m_gradRNeg,
				 members_t::m_gradRPos,
				 members_t::m_functionValues,
				 l1i, l0i, r0i, r1i,
				 uIndex,
				 ndpc);
	break;
      }

      case ReconstructionScheme::Weno5:
	{
	  throw std::runtime_error("missing impl");
	  break;
	}
      }
  }
};


}} // end namespaces
#endif
