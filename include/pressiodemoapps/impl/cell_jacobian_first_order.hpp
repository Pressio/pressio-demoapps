
#ifndef PRESSIODEMOAPPS_EDGE_CELL_JACOBIAN_FIRST_ORDER_HPP_
#define PRESSIODEMOAPPS_EDGE_CELL_JACOBIAN_FIRST_ORDER_HPP_

namespace pressiodemoapps{ namespace impl{

template<int dim, class IndexT, class GraphType>
void _get_left_right_cells_indices(IndexT & l0, IndexT & r0,
				   const IndexT & smPt,
				   const GraphType & graph,
				   int axis)
{
  if (dim == 1){
    l0 = graph(smPt, 1);
    r0 = graph(smPt, 2);
  }
  else if(dim==2){
    l0 = (axis == 1) ? graph(smPt, 1) : graph(smPt, 4);
    r0 = (axis == 1) ? graph(smPt, 3) : graph(smPt, 2);
  }
  else if(dim==3){
    l0 = (axis == 1) ? graph(smPt, 1) : (axis==2) ? graph(smPt, 4) : graph(smPt, 5);
    r0 = (axis == 1) ? graph(smPt, 3) : (axis==2) ? graph(smPt, 2) : graph(smPt, 6);
  }
}

template<
  class ScalarType,
  class SparseMatrix,
  class CellJacobianType,
  class MeshType
  >
struct CellJacobianMembers
{
  // used to figure out which direction we are dealing with
  // if ==1, then we are doing "x"
  // if ==2, then we are doing "y"
  // if ==3, then we are doing "z"
  int m_axis = 1;

  SparseMatrix & m_J;
  const MeshType & m_meshObj;
  const CellJacobianType & m_JLneg;
  const CellJacobianType & m_JLpos;
  const CellJacobianType & m_JRneg;
  const CellJacobianType & m_JRpos;

  // m_hInv = 1/h where h is the cell width
  // depends on which axis we are doing
  ScalarType m_hInv = {};

  CellJacobianMembers(SparseMatrix & J,
		      const MeshType & meshObj,
		      const CellJacobianType & JLneg,
		      const CellJacobianType & JLpos,
		      const CellJacobianType & JRneg,
		      const CellJacobianType & JRpos,
		      int axis = 1)
    : m_axis(axis), m_J(J), m_meshObj(meshObj),
      m_JLneg(JLneg), m_JLpos(JLpos), m_JRneg(JRneg), m_JRpos(JRpos)
  {
    m_hInv = (m_axis == 1) ? m_meshObj.dxInv() : (m_axis==2) ? m_meshObj.dyInv() : m_meshObj.dzInv();
  }
};


template<
  int dim, int numDofPerCell,
  class ScalarType,
  class SparseMatrix,
  class CellJacobianType,
  class MeshType
  >
class FirstOrderInnerCellJacobianFunctor
  : private CellJacobianMembers<ScalarType, SparseMatrix, CellJacobianType, MeshType>
{
  using members = CellJacobianMembers<ScalarType, SparseMatrix, CellJacobianType, MeshType>;
  using members::m_axis;
  using members::m_J;
  using members::m_meshObj;
  using members::m_JLneg;
  using members::m_JLpos;
  using members::m_JRneg;
  using members::m_JRpos;
  using members::m_hInv;

public:
  FirstOrderInnerCellJacobianFunctor(SparseMatrix & J,
				     const MeshType & meshObj,
				     const CellJacobianType & JLneg,
				     const CellJacobianType & JLpos,
				     const CellJacobianType & JRneg,
				     const CellJacobianType & JRpos,
				     int axis = 1)
    : members(J, meshObj, JLneg, JLpos, JRneg, JRpos){}

  template<class index_t>
  void operator()(const index_t smPt)
  {
    const auto & graph = m_meshObj.graph();
    index_t l0 = {};
    index_t r0 = {};
    _get_left_right_cells_indices<dim>(l0, r0, smPt, graph, m_axis);

    index_t rowIndex  = smPt*numDofPerCell;
    index_t col_im1   = l0*numDofPerCell;
    index_t col_i     = graph(smPt, 0)*numDofPerCell;
    index_t col_ip1   = r0*numDofPerCell;

    for (int k=0; k<numDofPerCell; ++k){
      for (int j=0; j<numDofPerCell; ++j){
	m_J.coeffRef(rowIndex+k, col_im1+j) +=  m_JLneg(k,j)*m_hInv;
	m_J.coeffRef(rowIndex+k, col_i+j)   += (m_JLpos(k,j)-m_JRneg(k,j))*m_hInv;
	m_J.coeffRef(rowIndex+k, col_ip1+j) += -m_JRpos(k,j)*m_hInv;
      }
    }
  }
};

template<
  int dim,
  class ScalarType,
  class SparseMatrix,
  class CellJacobianType,
  class MeshType
  >
class FirstOrderInnerCellJacobianFunctor<
  dim, 1, ScalarType, SparseMatrix, CellJacobianType, MeshType
  >
  : private CellJacobianMembers<ScalarType, SparseMatrix, CellJacobianType, MeshType>
{
  using members = CellJacobianMembers<ScalarType, SparseMatrix, CellJacobianType, MeshType>;
  using members::m_axis;
  using members::m_J;
  using members::m_meshObj;
  using members::m_JLneg;
  using members::m_JLpos;
  using members::m_JRneg;
  using members::m_JRpos;
  using members::m_hInv;

public:
  FirstOrderInnerCellJacobianFunctor(SparseMatrix & J,
				     const MeshType & meshObj,
				     const CellJacobianType & JLneg,
				     const CellJacobianType & JLpos,
				     const CellJacobianType & JRneg,
				     const CellJacobianType & JRpos,
				     int axis = 1)
    : members(J, meshObj, JLneg, JLpos, JRneg, JRpos){}

  template<class index_t>
  void operator()(const index_t smPt)
  {
    const auto & graph = m_meshObj.graph();
    index_t l0 = {};
    index_t r0 = {};
    _get_left_right_cells_indices<dim>(l0, r0, smPt, graph, m_axis);

    constexpr int numDofPerCell = 1;
    index_t rowIndex  = smPt*numDofPerCell;
    index_t col_im1   = l0*numDofPerCell;
    index_t col_i     = graph(smPt, 0)*numDofPerCell;
    index_t col_ip1   = r0*numDofPerCell;

    m_J.coeffRef(rowIndex, col_im1) +=  m_JLneg*m_hInv;
    m_J.coeffRef(rowIndex, col_i)   += (m_JLpos-m_JRneg)*m_hInv;
    m_J.coeffRef(rowIndex, col_ip1) += -m_JRpos*m_hInv;
  }
};


template<
  int dim, int numDofPerCell,
  class ScalarType,
  class SparseMatrix,
  class CellJacobianType,
  class MeshType
  >
class FirstOrderBdCellJacobianFunctor
  : private CellJacobianMembers<ScalarType, SparseMatrix, CellJacobianType, MeshType>
{
  using members = CellJacobianMembers<ScalarType, SparseMatrix, CellJacobianType, MeshType>;
  using members::m_axis;
  using members::m_J;
  using members::m_meshObj;
  using members::m_JLneg;
  using members::m_JLpos;
  using members::m_JRneg;
  using members::m_JRpos;
  using members::m_hInv;

public:
  FirstOrderBdCellJacobianFunctor(SparseMatrix & J,
				const MeshType & meshObj,
				const CellJacobianType & JLneg, const CellJacobianType & JLpos,
				const CellJacobianType & JRneg, const CellJacobianType & JRpos,
				int axis = 1)
    : members(J, meshObj, JLneg, JLpos, JRneg, JRpos){}

  template<class index_t>
  void operator()(const index_t smPt,
		  const std::array<ScalarType, numDofPerCell> & factors,
		  int bc_type = 0)
  {
    const auto & graph = m_meshObj.graph();
    index_t l0 = {};
    index_t r0 = {};
    _get_left_right_cells_indices<dim>(l0, r0, smPt, graph, m_axis);

    auto rowIndex  = smPt*numDofPerCell;
    auto col_i     = graph(smPt, 0)*numDofPerCell;

    for (int k=0; k<numDofPerCell; ++k){
      for (int j=0; j<numDofPerCell; ++j){
	m_J.coeffRef(rowIndex+k, col_i+j) += (m_JLpos(k,j) - m_JRneg(k,j))*m_hInv;
      }
    }

    if (l0 != -1){
      auto col_im1 = l0*numDofPerCell;
      for (int k=0; k<numDofPerCell; ++k){
	for (int j=0; j<numDofPerCell; ++j){
	  m_J.coeffRef(rowIndex+k, col_im1+j) += m_JLneg(k,j)*m_hInv;
	}
      }
    }

    if (r0 != -1){
      auto col_ip1 = r0*numDofPerCell;
      for (int k=0; k<numDofPerCell; ++k){
	for (int j=0; j<numDofPerCell; ++j){
	  m_J.coeffRef(rowIndex+k, col_ip1+j) += -m_JRpos(k,j)*m_hInv;
	}
      }
    }

    if (bc_type != 2)
    {
      if (l0 == -1){
	for (int k=0; k<numDofPerCell; ++k){
	  for (int j=0; j<numDofPerCell; ++j){
	    m_J.coeffRef(rowIndex+k, col_i+j) += (factors[j]*m_JLneg(k,j))*m_hInv;
	  }
	}
      }

      if (r0 == -1){
	for (int k=0; k<numDofPerCell; ++k){
	  for (int j=0; j<numDofPerCell; ++j){
	    m_J.coeffRef(rowIndex+k, col_i+j) += (factors[j]*-m_JRpos(k,j))*m_hInv;
	  }
	}
      }
    }
  }
};

}}
#endif
