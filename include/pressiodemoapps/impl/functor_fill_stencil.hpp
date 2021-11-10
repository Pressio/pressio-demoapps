
#ifndef PRESSIODEMOAPPS_STENCIL_FILLER_HPP_
#define PRESSIODEMOAPPS_STENCIL_FILLER_HPP_

namespace pressiodemoapps{ namespace impl{

template<
  int dim,
  int numDofPerCell,
  class stencil_values_t,
  class state_t,
  class mesh_t,
  class ghost_t
  >
class StencilFiller;

//------------------------------
// dim=1, 1 dofs/cell
//------------------------------
template<class stencil_values_t, class state_t, class mesh_t, class ghost_t>
class StencilFiller<
  1, 1, stencil_values_t, state_t, mesh_t, ghost_t
  >
{

public:
  StencilFiller() = delete;
  StencilFiller(const int stencilSize,
		const state_t & stateIn,
		const mesh_t & meshIn,
		const ghost_t & ghostLeft,
		const ghost_t & ghostRight,
		stencil_values_t & stencilVals)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostRight(ghostRight),
      m_stencilVals(stencilVals)
  {}

  template<class index_t>
  void operator()(const index_t smPt, int ghostRow)
  {
    constexpr int numDofPerCell = 1;

    const auto & graph = m_meshObj.graph();
    const auto uIndex  = graph(smPt, 0)*numDofPerCell;
    const auto w0      = graph(smPt, 1);
    const auto e0      = graph(smPt, 2);
    const auto w0i     = w0*numDofPerCell;
    const auto e0i     = e0*numDofPerCell;

    switch(m_stencilSize)
      {
      case 3:
	{
	  m_stencilVals(0) = (w0==-1) ? m_ghostLeft(ghostRow, 0)  : m_state(w0i);
	  m_stencilVals(1) = m_state(uIndex);
	  m_stencilVals(2) = (e0==-1) ? m_ghostRight(ghostRow, 0) : m_state(e0i);
	  break;
	}
      }
  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  const ghost_t & m_ghostLeft;
  const ghost_t & m_ghostRight;
  stencil_values_t & m_stencilVals;
};


//------------------------------
// dim=1, 3 dofs/cell
//------------------------------
template<class stencil_values_t, class state_t, class mesh_t, class ghost_t>
class StencilFiller<
  1, 3, stencil_values_t, state_t, mesh_t, ghost_t
  >
{

public:
  StencilFiller() = delete;
  StencilFiller(const int stencilSize,
		const state_t & stateIn,
		const mesh_t & meshIn,
		const ghost_t & ghostLeft,
		const ghost_t & ghostRight,
		stencil_values_t & stencilVals)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostRight(ghostRight),
      m_stencilVals(stencilVals)
  {}

  template<class index_t>
  void operator()(const index_t smPt, int ghostRow)
  {
    constexpr int numDofPerCell = 3;

    const auto & graph = m_meshObj.graph();
    const auto uIndex  = graph(smPt, 0)*numDofPerCell;
    const auto w0      = graph(smPt, 1);
    const auto e0      = graph(smPt, 2);
    const auto w0i     = w0*numDofPerCell;
    const auto e0i     = e0*numDofPerCell;

    switch(m_stencilSize)
      {
      case 3:
	{
	  m_stencilVals(0) = (w0==-1) ? m_ghostLeft(ghostRow,0) : m_state(w0i);
	  m_stencilVals(1) = (w0==-1) ? m_ghostLeft(ghostRow,1) : m_state(w0i+1);
	  m_stencilVals(2) = (w0==-1) ? m_ghostLeft(ghostRow,2) : m_state(w0i+2);
	  m_stencilVals(3) = m_state(uIndex);
	  m_stencilVals(4) = m_state(uIndex+1);
	  m_stencilVals(5) = m_state(uIndex+2);
	  m_stencilVals(6) = (e0==-1) ? m_ghostRight(ghostRow,0) : m_state(e0i);
	  m_stencilVals(7) = (e0==-1) ? m_ghostRight(ghostRow,1) : m_state(e0i+1);
	  m_stencilVals(8) = (e0==-1) ? m_ghostRight(ghostRow,2) : m_state(e0i+2);
	  break;
	}

      case 5:
	{
	  const auto w1 = graph(smPt, 3);
	  const auto e1 = graph(smPt, 4);
	  const auto w1i = w1*numDofPerCell;
	  const auto e1i = e1*numDofPerCell;

	  m_stencilVals(0)  = (w1==-1) ? m_ghostLeft(ghostRow,3) : m_state(w1i);
	  m_stencilVals(1)  = (w1==-1) ? m_ghostLeft(ghostRow,4) : m_state(w1i+1);
	  m_stencilVals(2)  = (w1==-1) ? m_ghostLeft(ghostRow,5) : m_state(w1i+2);
	  m_stencilVals(3)  = (w0==-1) ? m_ghostLeft(ghostRow,0) : m_state(w0i);
	  m_stencilVals(4)  = (w0==-1) ? m_ghostLeft(ghostRow,1) : m_state(w0i+1);
	  m_stencilVals(5)  = (w0==-1) ? m_ghostLeft(ghostRow,2) : m_state(w0i+2);
	  m_stencilVals(6)  = m_state(uIndex);
	  m_stencilVals(7)  = m_state(uIndex+1);
	  m_stencilVals(8)  = m_state(uIndex+2);
	  m_stencilVals(9)  = (e0==-1) ? m_ghostRight(ghostRow,0) : m_state(e0i);
	  m_stencilVals(10) = (e0==-1) ? m_ghostRight(ghostRow,1) : m_state(e0i+1);
	  m_stencilVals(11) = (e0==-1) ? m_ghostRight(ghostRow,2) : m_state(e0i+2);
	  m_stencilVals(12) = (e1==-1) ? m_ghostRight(ghostRow,3) : m_state(e1i);
	  m_stencilVals(13) = (e1==-1) ? m_ghostRight(ghostRow,4) : m_state(e1i+1);
	  m_stencilVals(14) = (e1==-1) ? m_ghostRight(ghostRow,5) : m_state(e1i+2);
	  break;
	}

      case 7:
	{
	  const auto w1 = graph(smPt, 3);
	  const auto e1 = graph(smPt, 4);
	  const auto w2 = graph(smPt, 5);
	  const auto e2 = graph(smPt, 6);
	  const auto w1i = w1*numDofPerCell;
	  const auto e1i = e1*numDofPerCell;
	  const auto w2i = w2*numDofPerCell;
	  const auto e2i = e2*numDofPerCell;

	  m_stencilVals(0)  = (w2==-1) ? m_ghostLeft(ghostRow,6) : m_state(w2i);
	  m_stencilVals(1)  = (w2==-1) ? m_ghostLeft(ghostRow,7) : m_state(w2i+1);
	  m_stencilVals(2)  = (w2==-1) ? m_ghostLeft(ghostRow,8) : m_state(w2i+2);
	  m_stencilVals(3)  = (w1==-1) ? m_ghostLeft(ghostRow,3) : m_state(w1i);
	  m_stencilVals(4)  = (w1==-1) ? m_ghostLeft(ghostRow,4) : m_state(w1i+1);
	  m_stencilVals(5)  = (w1==-1) ? m_ghostLeft(ghostRow,5) : m_state(w1i+2);
	  m_stencilVals(6)  = (w0==-1) ? m_ghostLeft(ghostRow,0) : m_state(w0i);
	  m_stencilVals(7)  = (w0==-1) ? m_ghostLeft(ghostRow,1) : m_state(w0i+1);
	  m_stencilVals(8)  = (w0==-1) ? m_ghostLeft(ghostRow,2) : m_state(w0i+2);
	  m_stencilVals(9)  = m_state(uIndex);
	  m_stencilVals(10) = m_state(uIndex+1);
	  m_stencilVals(11) = m_state(uIndex+2);
	  m_stencilVals(12) = (e0==-1) ? m_ghostRight(ghostRow,0) : m_state(e0i);
	  m_stencilVals(13) = (e0==-1) ? m_ghostRight(ghostRow,1) : m_state(e0i+1);
	  m_stencilVals(14) = (e0==-1) ? m_ghostRight(ghostRow,2) : m_state(e0i+2);
	  m_stencilVals(15) = (e1==-1) ? m_ghostRight(ghostRow,3) : m_state(e1i);
	  m_stencilVals(16) = (e1==-1) ? m_ghostRight(ghostRow,4) : m_state(e1i+1);
	  m_stencilVals(17) = (e1==-1) ? m_ghostRight(ghostRow,5) : m_state(e1i+2);
	  m_stencilVals(18) = (e2==-1) ? m_ghostRight(ghostRow,6) : m_state(e2i);
	  m_stencilVals(19) = (e2==-1) ? m_ghostRight(ghostRow,7) : m_state(e2i+1);
	  m_stencilVals(20) = (e2==-1) ? m_ghostRight(ghostRow,8) : m_state(e2i+2);
	  break;
	}

      }
  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  const ghost_t & m_ghostLeft;
  const ghost_t & m_ghostRight;
  stencil_values_t & m_stencilVals;
};


//------------------------------
// dim=2, 1 dofs/cell
//------------------------------
template<class stencil_values_t, class state_t, class mesh_t, class ghost_t>
class StencilFiller<
  2, 1, stencil_values_t, state_t, mesh_t, ghost_t
  >
{

public:
  StencilFiller() = delete;
  StencilFiller(const int stencilSize,
		const state_t & stateIn,
		const mesh_t & meshIn,
		const ghost_t & ghostLeft,
		const ghost_t & ghostRight,
		stencil_values_t & stencilVals,
		const int axis)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostRight(ghostRight),
      m_stencilVals(stencilVals),
      m_axis(axis)
  {}

  template<class index_t>
  void operator()(const index_t smPt, int ghostRow)
  {
    constexpr auto numDofPerCell = 1;
    const auto & graph = m_meshObj.graph();
    const auto uIndex = graph(smPt, 0)*numDofPerCell;

    switch(m_stencilSize)
    {
    case 3:
      {
	const auto l0 = (m_axis == 1) ? graph(smPt, 1) : graph(smPt, 4);
	const auto r0 = (m_axis == 1) ? graph(smPt, 3) : graph(smPt, 2);

	if (l0 == -1){
	  m_stencilVals(0) = m_ghostLeft(ghostRow, 0);
	}else{
	  const auto index = l0*numDofPerCell;
	  m_stencilVals(0) = m_state(index);
	}

	m_stencilVals(1)  = m_state(uIndex);

	if (r0 == -1){
	  m_stencilVals(2)  = m_ghostRight(ghostRow, 0);
	}else{
	  const auto index = r0*numDofPerCell;
	  m_stencilVals(2)  = m_state(index);
	}
	break;
      } // case 3
    }
  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  const ghost_t & m_ghostLeft;
  const ghost_t & m_ghostRight;
  stencil_values_t & m_stencilVals;
  int m_axis;
};


//------------------------------
// dim=2, 3 dofs/cell
//------------------------------
template<class stencil_values_t, class state_t, class mesh_t, class ghost_t>
class StencilFiller<
  2, 3, stencil_values_t, state_t, mesh_t, ghost_t
  >
{

public:
  StencilFiller() = delete;
  StencilFiller(const int stencilSize,
		const state_t & stateIn,
		const mesh_t & meshIn,
		const ghost_t & ghostLeft,
		const ghost_t & ghostRight,
		stencil_values_t & stencilVals,
		const int axis)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostRight(ghostRight),
      m_stencilVals(stencilVals),
      m_axis(axis)
  {}

  template<class index_t>
  void operator()(const index_t smPt, int ghostRow)
  {
    constexpr auto numDofPerCell = 3;
    const auto & graph = m_meshObj.graph();
    const auto uIndex = graph(smPt, 0)*numDofPerCell;

    switch(m_stencilSize)
    {
    case 3:
      {
	const auto l0 = (m_axis == 1) ? graph(smPt, 1) : graph(smPt, 4);
	const auto r0 = (m_axis == 1) ? graph(smPt, 3) : graph(smPt, 2);

	if (l0 == -1){
	  m_stencilVals(0) = m_ghostLeft(ghostRow, 0);
	  m_stencilVals(1) = m_ghostLeft(ghostRow, 1);
	  m_stencilVals(2) = m_ghostLeft(ghostRow, 2);
	}else{
	  const auto index = l0*3;
	  m_stencilVals(0) = m_state(index);
	  m_stencilVals(1) = m_state(index+1);
	  m_stencilVals(2) = m_state(index+2);
	}

	m_stencilVals(3)  = m_state(uIndex);
	m_stencilVals(4)  = m_state(uIndex+1);
	m_stencilVals(5)  = m_state(uIndex+2);

	if (r0 == -1){
	  m_stencilVals(6)  = m_ghostRight(ghostRow, 0);
	  m_stencilVals(7)  = m_ghostRight(ghostRow, 1);
	  m_stencilVals(8) = m_ghostRight(ghostRow, 2);
	}else{
	  const auto index = r0*3;
	  m_stencilVals(6)  = m_state(index);
	  m_stencilVals(7)  = m_state(index+1);
	  m_stencilVals(8) = m_state(index+2);
	}
	break;
      } // case 3

    case 5:
      {
	const auto l0 = (m_axis == 1) ? graph(smPt, 1)  : graph(smPt, 4);
	const auto r0 = (m_axis == 1) ? graph(smPt, 3)  : graph(smPt, 2);
	const auto l1 = (m_axis == 1) ? graph(smPt, 5)  : graph(smPt, 8);
	const auto r1 = (m_axis == 1) ? graph(smPt, 7)  : graph(smPt, 6);
	const auto uIndex = graph(smPt, 0)*numDofPerCell;

	if (l1 == -1){
	  m_stencilVals(0) = m_ghostLeft(ghostRow, 3);
	  m_stencilVals(1) = m_ghostLeft(ghostRow, 4);
	  m_stencilVals(2) = m_ghostLeft(ghostRow, 5);
	}else{
	  const auto index = l1*numDofPerCell;
	  m_stencilVals(0) = m_state(index);
	  m_stencilVals(1) = m_state(index+1);
	  m_stencilVals(2) = m_state(index+2);
	}

	if (l0 == -1){
	  m_stencilVals(3) = m_ghostLeft(ghostRow, 0);
	  m_stencilVals(4) = m_ghostLeft(ghostRow, 1);
	  m_stencilVals(5) = m_ghostLeft(ghostRow, 2);
	}else{
	  const auto index = l0*numDofPerCell;
	  m_stencilVals(3) = m_state(index);
	  m_stencilVals(4) = m_state(index+1);
	  m_stencilVals(5) = m_state(index+2);
	}

	m_stencilVals(6)  = m_state(uIndex);
	m_stencilVals(7)  = m_state(uIndex+1);
	m_stencilVals(8) = m_state(uIndex+2);

	if (r0 == -1){
	  m_stencilVals(9) = m_ghostRight(ghostRow, 0);
	  m_stencilVals(10) = m_ghostRight(ghostRow, 1);
	  m_stencilVals(11) = m_ghostRight(ghostRow, 2);
	}else{
	  const auto index = r0*numDofPerCell;
	  m_stencilVals(9) = m_state(index);
	  m_stencilVals(10) = m_state(index+1);
	  m_stencilVals(11) = m_state(index+2);
	}

	if (r1 == -1){
	  m_stencilVals(12) = m_ghostRight(ghostRow, 3);
	  m_stencilVals(13) = m_ghostRight(ghostRow, 4);
	  m_stencilVals(14) = m_ghostRight(ghostRow, 5);
	}else{
	  const auto index = r1*numDofPerCell;
	  m_stencilVals(12) = m_state(index);
	  m_stencilVals(13) = m_state(index+1);
	  m_stencilVals(14) = m_state(index+2);
	}

	break;
      }//case 5

    case 7:
      {
	const auto l0 = (m_axis == 1) ? graph(smPt, 1)  : graph(smPt, 4);
	const auto r0 = (m_axis == 1) ? graph(smPt, 3)  : graph(smPt, 2);
	const auto l1 = (m_axis == 1) ? graph(smPt, 5)  : graph(smPt, 8);
	const auto r1 = (m_axis == 1) ? graph(smPt, 7)  : graph(smPt, 6);
	const auto l2 = (m_axis == 1) ? graph(smPt, 9)  : graph(smPt, 12);
	const auto r2 = (m_axis == 1) ? graph(smPt, 11) : graph(smPt, 10);
	const auto uIndex = graph(smPt, 0)*numDofPerCell;

	if (l2 == -1){
	  m_stencilVals(0) = m_ghostLeft(ghostRow, 6);
	  m_stencilVals(1) = m_ghostLeft(ghostRow, 7);
	  m_stencilVals(2) = m_ghostLeft(ghostRow, 8);
	}else{
	  const auto index = l2*numDofPerCell;
	  m_stencilVals(0) = m_state(index);
	  m_stencilVals(1) = m_state(index+1);
	  m_stencilVals(2) = m_state(index+2);
	  m_stencilVals(3) = m_state(index+3);
	}

	if (l1 == -1){
	  m_stencilVals(3) = m_ghostLeft(ghostRow, 3);
	  m_stencilVals(4) = m_ghostLeft(ghostRow, 4);
	  m_stencilVals(5) = m_ghostLeft(ghostRow, 5);
	}else{
	  const auto index = l1*numDofPerCell;
	  m_stencilVals(3) = m_state(index);
	  m_stencilVals(4) = m_state(index+1);
	  m_stencilVals(5) = m_state(index+2);
	}

	if (l0 == -1){
	  m_stencilVals(6)  = m_ghostLeft(ghostRow, 0);
	  m_stencilVals(7)  = m_ghostLeft(ghostRow, 1);
	  m_stencilVals(8) = m_ghostLeft(ghostRow, 2);
	}else{
	  const auto index = l0*numDofPerCell;
	  m_stencilVals(6)  = m_state(index);
	  m_stencilVals(7)  = m_state(index+1);
	  m_stencilVals(8) = m_state(index+2);
	}

	m_stencilVals(9) = m_state(uIndex);
	m_stencilVals(10) = m_state(uIndex+1);
	m_stencilVals(11) = m_state(uIndex+2);

	if (r0 == -1){
	  m_stencilVals(12) = m_ghostRight(ghostRow, 0);
	  m_stencilVals(13) = m_ghostRight(ghostRow, 1);
	  m_stencilVals(14) = m_ghostRight(ghostRow, 2);
	}else{
	  const auto index = r0*numDofPerCell;
	  m_stencilVals(12) = m_state(index);
	  m_stencilVals(13) = m_state(index+1);
	  m_stencilVals(14) = m_state(index+2);
	}

	if (r1 == -1){
	  m_stencilVals(15) = m_ghostRight(ghostRow, 3);
	  m_stencilVals(16) = m_ghostRight(ghostRow, 4);
	  m_stencilVals(17) = m_ghostRight(ghostRow, 5);
	}else{
	  const auto index = r1*numDofPerCell;
	  m_stencilVals(15) = m_state(index);
	  m_stencilVals(16) = m_state(index+1);
	  m_stencilVals(17) = m_state(index+2);
	}

	if (r2 == -1){
	  m_stencilVals(18) = m_ghostRight(ghostRow, 6);
	  m_stencilVals(19) = m_ghostRight(ghostRow, 7);
	  m_stencilVals(20) = m_ghostRight(ghostRow, 8);
	}else{
	  const auto index = r2*numDofPerCell;
	  m_stencilVals(18) = m_state(index);
	  m_stencilVals(19) = m_state(index+1);
	  m_stencilVals(20) = m_state(index+2);
	}
	break;
      }//case 7
    }
  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  const ghost_t & m_ghostLeft;
  const ghost_t & m_ghostRight;
  stencil_values_t & m_stencilVals;
  int m_axis;
};


//------------------------------
// dim=2, 4 dofs/cell
//------------------------------
template<class stencil_values_t, class state_t, class mesh_t, class ghost_t>
class StencilFiller<
  2, 4, stencil_values_t, state_t, mesh_t, ghost_t
  >
{

public:
  StencilFiller() = delete;
  StencilFiller(const int stencilSize,
		const state_t & stateIn,
		const mesh_t & meshIn,
		const ghost_t & ghostLeft,
		const ghost_t & ghostRight,
		stencil_values_t & stencilVals,
		const int axis)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostRight(ghostRight),
      m_stencilVals(stencilVals),
      m_axis(axis)
  {}

  template<class index_t>
  void operator()(const index_t smPt, int ghostRow)
  {
    constexpr auto numDofPerCell = 4;
    const auto & graph = m_meshObj.graph();
    const auto uIndex = graph(smPt, 0)*numDofPerCell;

    switch(m_stencilSize)
    {
    case 3:
      {
	const auto l0 = (m_axis == 1) ? graph(smPt, 1) : graph(smPt, 4);
	const auto r0 = (m_axis == 1) ? graph(smPt, 3) : graph(smPt, 2);

	if (l0 == -1){
	  m_stencilVals(0) = m_ghostLeft(ghostRow, 0);
	  m_stencilVals(1) = m_ghostLeft(ghostRow, 1);
	  m_stencilVals(2) = m_ghostLeft(ghostRow, 2);
	  m_stencilVals(3) = m_ghostLeft(ghostRow, 3);
	}else{
	  const auto index = l0*4;
	  m_stencilVals(0) = m_state(index);
	  m_stencilVals(1) = m_state(index+1);
	  m_stencilVals(2) = m_state(index+2);
	  m_stencilVals(3) = m_state(index+3);
	}

	m_stencilVals(4)  = m_state(uIndex);
	m_stencilVals(5)  = m_state(uIndex+1);
	m_stencilVals(6)  = m_state(uIndex+2);
	m_stencilVals(7)  = m_state(uIndex+3);

	if (r0 == -1){
	  m_stencilVals(8)  = m_ghostRight(ghostRow, 0);
	  m_stencilVals(9)  = m_ghostRight(ghostRow, 1);
	  m_stencilVals(10) = m_ghostRight(ghostRow, 2);
	  m_stencilVals(11) = m_ghostRight(ghostRow, 3);
	}else{
	  const auto index = r0*4;
	  m_stencilVals(8)  = m_state(index);
	  m_stencilVals(9)  = m_state(index+1);
	  m_stencilVals(10) = m_state(index+2);
	  m_stencilVals(11) = m_state(index+3);
	}
	break;
      } // case 3

    case 5:
      {
	const auto l0 = (m_axis == 1) ? graph(smPt, 1)  : graph(smPt, 4);
	const auto r0 = (m_axis == 1) ? graph(smPt, 3)  : graph(smPt, 2);
	const auto l1 = (m_axis == 1) ? graph(smPt, 5)  : graph(smPt, 8);
	const auto r1 = (m_axis == 1) ? graph(smPt, 7)  : graph(smPt, 6);
	const auto uIndex = graph(smPt, 0)*numDofPerCell;

	if (l1 == -1){
	  m_stencilVals(0) = m_ghostLeft(ghostRow, 4);
	  m_stencilVals(1) = m_ghostLeft(ghostRow, 5);
	  m_stencilVals(2) = m_ghostLeft(ghostRow, 6);
	  m_stencilVals(3) = m_ghostLeft(ghostRow, 7);
	}else{
	  const auto index = l1*numDofPerCell;
	  m_stencilVals(0) = m_state(index);
	  m_stencilVals(1) = m_state(index+1);
	  m_stencilVals(2) = m_state(index+2);
	  m_stencilVals(3) = m_state(index+3);
	}

	if (l0 == -1){
	  m_stencilVals(4) = m_ghostLeft(ghostRow, 0);
	  m_stencilVals(5) = m_ghostLeft(ghostRow, 1);
	  m_stencilVals(6) = m_ghostLeft(ghostRow, 2);
	  m_stencilVals(7) = m_ghostLeft(ghostRow, 3);
	}else{
	  const auto index = l0*numDofPerCell;
	  m_stencilVals(4) = m_state(index);
	  m_stencilVals(5) = m_state(index+1);
	  m_stencilVals(6) = m_state(index+2);
	  m_stencilVals(7) = m_state(index+3);
	}

	m_stencilVals(8)  = m_state(uIndex);
	m_stencilVals(9)  = m_state(uIndex+1);
	m_stencilVals(10) = m_state(uIndex+2);
	m_stencilVals(11) = m_state(uIndex+3);

	if (r0 == -1){
	  m_stencilVals(12) = m_ghostRight(ghostRow, 0);
	  m_stencilVals(13) = m_ghostRight(ghostRow, 1);
	  m_stencilVals(14) = m_ghostRight(ghostRow, 2);
	  m_stencilVals(15) = m_ghostRight(ghostRow, 3);
	}else{
	  const auto index = r0*numDofPerCell;
	  m_stencilVals(12) = m_state(index);
	  m_stencilVals(13) = m_state(index+1);
	  m_stencilVals(14) = m_state(index+2);
	  m_stencilVals(15) = m_state(index+3);
	}

	if (r1 == -1){
	  m_stencilVals(16) = m_ghostRight(ghostRow, 4);
	  m_stencilVals(17) = m_ghostRight(ghostRow, 5);
	  m_stencilVals(18) = m_ghostRight(ghostRow, 6);
	  m_stencilVals(19) = m_ghostRight(ghostRow, 7);
	}else{
	  const auto index = r1*numDofPerCell;
	  m_stencilVals(16) = m_state(index);
	  m_stencilVals(17) = m_state(index+1);
	  m_stencilVals(18) = m_state(index+2);
	  m_stencilVals(19) = m_state(index+3);
	}

	break;
      }//case 5

    case 7:
      {
	const auto l0 = (m_axis == 1) ? graph(smPt, 1)  : graph(smPt, 4);
	const auto r0 = (m_axis == 1) ? graph(smPt, 3)  : graph(smPt, 2);
	const auto l1 = (m_axis == 1) ? graph(smPt, 5)  : graph(smPt, 8);
	const auto r1 = (m_axis == 1) ? graph(smPt, 7)  : graph(smPt, 6);
	const auto l2 = (m_axis == 1) ? graph(smPt, 9)  : graph(smPt, 12);
	const auto r2 = (m_axis == 1) ? graph(smPt, 11) : graph(smPt, 10);
	const auto uIndex = graph(smPt, 0)*numDofPerCell;

	if (l2 == -1){
	  m_stencilVals(0) = m_ghostLeft(ghostRow, 8);
	  m_stencilVals(1) = m_ghostLeft(ghostRow, 9);
	  m_stencilVals(2) = m_ghostLeft(ghostRow, 10);
	  m_stencilVals(3) = m_ghostLeft(ghostRow, 11);
	}else{
	  const auto index = l2*numDofPerCell;
	  m_stencilVals(0) = m_state(index);
	  m_stencilVals(1) = m_state(index+1);
	  m_stencilVals(2) = m_state(index+2);
	  m_stencilVals(3) = m_state(index+3);
	}

	if (l1 == -1){
	  m_stencilVals(4) = m_ghostLeft(ghostRow, 4);
	  m_stencilVals(5) = m_ghostLeft(ghostRow, 5);
	  m_stencilVals(6) = m_ghostLeft(ghostRow, 6);
	  m_stencilVals(7) = m_ghostLeft(ghostRow, 7);
	}else{
	  const auto index = l1*numDofPerCell;
	  m_stencilVals(4) = m_state(index);
	  m_stencilVals(5) = m_state(index+1);
	  m_stencilVals(6) = m_state(index+2);
	  m_stencilVals(7) = m_state(index+3);
	}

	if (l0 == -1){
	  m_stencilVals(8)  = m_ghostLeft(ghostRow, 0);
	  m_stencilVals(9)  = m_ghostLeft(ghostRow, 1);
	  m_stencilVals(10) = m_ghostLeft(ghostRow, 2);
	  m_stencilVals(11) = m_ghostLeft(ghostRow, 3);
	}else{
	  const auto index = l0*numDofPerCell;
	  m_stencilVals(8)  = m_state(index);
	  m_stencilVals(9)  = m_state(index+1);
	  m_stencilVals(10) = m_state(index+2);
	  m_stencilVals(11) = m_state(index+3);
	}

	m_stencilVals(12) = m_state(uIndex);
	m_stencilVals(13) = m_state(uIndex+1);
	m_stencilVals(14) = m_state(uIndex+2);
	m_stencilVals(15) = m_state(uIndex+3);

	if (r0 == -1){
	  m_stencilVals(16) = m_ghostRight(ghostRow, 0);
	  m_stencilVals(17) = m_ghostRight(ghostRow, 1);
	  m_stencilVals(18) = m_ghostRight(ghostRow, 2);
	  m_stencilVals(19) = m_ghostRight(ghostRow, 3);
	}else{
	  const auto index = r0*numDofPerCell;
	  m_stencilVals(16) = m_state(index);
	  m_stencilVals(17) = m_state(index+1);
	  m_stencilVals(18) = m_state(index+2);
	  m_stencilVals(19) = m_state(index+3);
	}

	if (r1 == -1){
	  m_stencilVals(20) = m_ghostRight(ghostRow, 4);
	  m_stencilVals(21) = m_ghostRight(ghostRow, 5);
	  m_stencilVals(22) = m_ghostRight(ghostRow, 6);
	  m_stencilVals(23) = m_ghostRight(ghostRow, 7);
	}else{
	  const auto index = r1*numDofPerCell;
	  m_stencilVals(20) = m_state(index);
	  m_stencilVals(21) = m_state(index+1);
	  m_stencilVals(22) = m_state(index+2);
	  m_stencilVals(23) = m_state(index+3);
	}

	if (r2 == -1){
	  m_stencilVals(24) = m_ghostRight(ghostRow, 8);
	  m_stencilVals(25) = m_ghostRight(ghostRow, 9);
	  m_stencilVals(26) = m_ghostRight(ghostRow, 10);
	  m_stencilVals(27) = m_ghostRight(ghostRow, 11);
	}else{
	  const auto index = r2*numDofPerCell;
	  m_stencilVals(24) = m_state(index);
	  m_stencilVals(25) = m_state(index+1);
	  m_stencilVals(26) = m_state(index+2);
	  m_stencilVals(27) = m_state(index+3);
	}
	break;
      }//case 7
    }
  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  const ghost_t & m_ghostLeft;
  const ghost_t & m_ghostRight;
  stencil_values_t & m_stencilVals;
  int m_axis;
};


//------------------------------
// dim=3, 5 dofs/cell
//------------------------------
template<class stencil_values_t, class state_t, class mesh_t, class ghost_t>
class StencilFiller<
  3, 5, stencil_values_t, state_t, mesh_t, ghost_t
  >
{

public:
  StencilFiller() = delete;
  StencilFiller(const int stencilSize,
		const state_t & stateIn,
		const mesh_t & meshIn,
		const ghost_t & ghostLeft,
		const ghost_t & ghostRight,
		stencil_values_t & stencilVals,
		const int axis)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostRight(ghostRight),
      m_stencilVals(stencilVals),
      m_axis(axis)
  {}

  template<class index_t>
  void operator()(const index_t smPt, int ghostRow)
  {
    constexpr auto numDofPerCell = 5;
    const auto & graph = m_meshObj.graph();
    const auto uIndex = graph(smPt, 0)*numDofPerCell;

    switch(m_stencilSize)
    {

    case 3:
      {
	auto l0 = (m_axis == 1) ? graph(smPt, 1) : (m_axis==2) ? graph(smPt, 4) : graph(smPt, 5);
	auto r0 = (m_axis == 1) ? graph(smPt, 3) : (m_axis==2) ? graph(smPt, 2) : graph(smPt, 6);

	if (l0 == -1){
	  m_stencilVals(0) = m_ghostLeft(ghostRow, 0);
	  m_stencilVals(1) = m_ghostLeft(ghostRow, 1);
	  m_stencilVals(2) = m_ghostLeft(ghostRow, 2);
	  m_stencilVals(3) = m_ghostLeft(ghostRow, 3);
	  m_stencilVals(4) = m_ghostLeft(ghostRow, 4);
	}else{
	  const auto index = l0*numDofPerCell;
	  m_stencilVals(0) = m_state(index);
	  m_stencilVals(1) = m_state(index+1);
	  m_stencilVals(2) = m_state(index+2);
	  m_stencilVals(3) = m_state(index+3);
	  m_stencilVals(4) = m_state(index+4);
	}

	m_stencilVals(5)  = m_state(uIndex);
	m_stencilVals(6)  = m_state(uIndex+1);
	m_stencilVals(7)  = m_state(uIndex+2);
	m_stencilVals(8)  = m_state(uIndex+3);
	m_stencilVals(9)  = m_state(uIndex+4);

	if (r0 == -1){
	  m_stencilVals(10) = m_ghostRight(ghostRow, 0);
	  m_stencilVals(11) = m_ghostRight(ghostRow, 1);
	  m_stencilVals(12) = m_ghostRight(ghostRow, 2);
	  m_stencilVals(13) = m_ghostRight(ghostRow, 3);
	  m_stencilVals(14) = m_ghostRight(ghostRow, 4);
	}else{
	  const auto index = r0*numDofPerCell;
	  m_stencilVals(10) = m_state(index);
	  m_stencilVals(11) = m_state(index+1);
	  m_stencilVals(12) = m_state(index+2);
	  m_stencilVals(13) = m_state(index+3);
	  m_stencilVals(14) = m_state(index+4);
	}
	break;
      } // case 3

    case 5:
      {
	auto l0 = (m_axis == 1) ? graph(smPt, 1) : (m_axis==2) ? graph(smPt, 4)  : graph(smPt, 5);
	auto r0 = (m_axis == 1) ? graph(smPt, 3) : (m_axis==2) ? graph(smPt, 2)  : graph(smPt, 6);
	auto l1 = (m_axis == 1) ? graph(smPt, 7) : (m_axis==2) ? graph(smPt, 10) : graph(smPt, 11);
	auto r1 = (m_axis == 1) ? graph(smPt, 9) : (m_axis==2) ? graph(smPt, 8)  : graph(smPt, 12);

	if (l1 == -1){
	  m_stencilVals(0) = m_ghostLeft(ghostRow, 5);
	  m_stencilVals(1) = m_ghostLeft(ghostRow, 6);
	  m_stencilVals(2) = m_ghostLeft(ghostRow, 7);
	  m_stencilVals(3) = m_ghostLeft(ghostRow, 8);
	  m_stencilVals(4) = m_ghostLeft(ghostRow, 9);
	}else{
	  const auto index = l1*numDofPerCell;
	  m_stencilVals(0) = m_state(index);
	  m_stencilVals(1) = m_state(index+1);
	  m_stencilVals(2) = m_state(index+2);
	  m_stencilVals(3) = m_state(index+3);
	  m_stencilVals(4) = m_state(index+4);
	}

	if (l0 == -1){
	  m_stencilVals(5) = m_ghostLeft(ghostRow, 0);
	  m_stencilVals(6) = m_ghostLeft(ghostRow, 1);
	  m_stencilVals(7) = m_ghostLeft(ghostRow, 2);
	  m_stencilVals(8) = m_ghostLeft(ghostRow, 3);
	  m_stencilVals(9) = m_ghostLeft(ghostRow, 4);
	}else{
	  const auto index = l0*numDofPerCell;
	  m_stencilVals(5) = m_state(index);
	  m_stencilVals(6) = m_state(index+1);
	  m_stencilVals(7) = m_state(index+2);
	  m_stencilVals(8) = m_state(index+3);
	  m_stencilVals(9) = m_state(index+4);
	}

	m_stencilVals(10)  = m_state(uIndex);
	m_stencilVals(11)  = m_state(uIndex+1);
	m_stencilVals(12)  = m_state(uIndex+2);
	m_stencilVals(13)  = m_state(uIndex+3);
	m_stencilVals(14)  = m_state(uIndex+4);

	if (r0 == -1){
	  m_stencilVals(15) = m_ghostRight(ghostRow, 0);
	  m_stencilVals(16) = m_ghostRight(ghostRow, 1);
	  m_stencilVals(17) = m_ghostRight(ghostRow, 2);
	  m_stencilVals(18) = m_ghostRight(ghostRow, 3);
	  m_stencilVals(19) = m_ghostRight(ghostRow, 4);
	}else{
	  const auto index = r0*numDofPerCell;
	  m_stencilVals(15) = m_state(index);
	  m_stencilVals(16) = m_state(index+1);
	  m_stencilVals(17) = m_state(index+2);
	  m_stencilVals(18) = m_state(index+3);
	  m_stencilVals(19) = m_state(index+4);
	}

	if (r1 == -1){
	  m_stencilVals(20) = m_ghostRight(ghostRow, 5);
	  m_stencilVals(21) = m_ghostRight(ghostRow, 6);
	  m_stencilVals(22) = m_ghostRight(ghostRow, 7);
	  m_stencilVals(23) = m_ghostRight(ghostRow, 8);
	  m_stencilVals(24) = m_ghostRight(ghostRow, 9);
	}else{
	  const auto index = r1*numDofPerCell;
	  m_stencilVals(20) = m_state(index);
	  m_stencilVals(21) = m_state(index+1);
	  m_stencilVals(22) = m_state(index+2);
	  m_stencilVals(23) = m_state(index+3);
	  m_stencilVals(24) = m_state(index+4);
	}

	break;
      } // case 5

    case 7:
      {
	throw std::runtime_error("For 3d, only 1st order is currently working");
	break;
      } // case 7

    }
  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  const ghost_t & m_ghostLeft;
  const ghost_t & m_ghostRight;
  stencil_values_t & m_stencilVals;
  int m_axis;
};

}}
#endif
