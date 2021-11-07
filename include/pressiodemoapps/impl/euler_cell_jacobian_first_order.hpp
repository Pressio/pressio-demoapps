
#ifndef PRESSIODEMOAPPS_EDGE_CELL_JACOBIAN_FIRST_ORDER_HPP_
#define PRESSIODEMOAPPS_EDGE_CELL_JACOBIAN_FIRST_ORDER_HPP_

namespace pressiodemoapps{ namespace ee{ namespace impl{

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
  int dim, int numDofPerCell,
  class ScalarType,
  class SparseMatrix,
  class CellJacobianType,
  class MeshType
  >
class FirstOrderInnerCellJacobianFunctor
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

public:
  FirstOrderInnerCellJacobianFunctor(SparseMatrix & J,
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
  int dim, int numDofPerCell,
  class ScalarType,
  class SparseMatrix,
  class CellJacobianType,
  class MeshType
  >
class FirstOrderCellJacobianFunctor
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

  std::array<ScalarType, numDofPerCell> m_factorsDefault;
  std::array<ScalarType, numDofPerCell> m_factorsReflective;

public:
  FirstOrderCellJacobianFunctor(SparseMatrix & J,
				const MeshType & meshObj,
				const CellJacobianType & JLneg, const CellJacobianType & JLpos,
				const CellJacobianType & JRneg, const CellJacobianType & JRpos,
				int axis = 1)
    : m_axis(axis), m_J(J), m_meshObj(meshObj),
      m_JLneg(JLneg), m_JLpos(JLpos), m_JRneg(JRneg), m_JRpos(JRpos)
  {
    m_hInv = (m_axis == 1) ? m_meshObj.dxInv() : (m_axis==2) ? m_meshObj.dyInv() : m_meshObj.dzInv();
    m_factorsDefault.fill(static_cast<ScalarType>(1));

    // prepare factors that we use below based on BC being reflective or not
    m_factorsReflective.fill(static_cast<ScalarType>(1));
    if (numDofPerCell == 3){
      m_factorsReflective[1] = static_cast<ScalarType>(-1);
    }
    else if (numDofPerCell == 4){
      if (m_axis==1){
	m_factorsReflective[1] = static_cast<ScalarType>(-1);
      }
      else{
	m_factorsReflective[2] = static_cast<ScalarType>(-1);
      }
    }
    else if (numDofPerCell == 5){
      if (m_axis==1){
	m_factorsReflective[1] = static_cast<ScalarType>(-1);
      }
      else if (m_axis==2){
	m_factorsReflective[2] = static_cast<ScalarType>(-1);
      }
      else if (m_axis==3){
	m_factorsReflective[3] = static_cast<ScalarType>(-1);
      }
    }
  }

  template<class index_t>
  void operator()(const index_t smPt, int bc_type = 0)
  {
    const auto & graph = m_meshObj.graph();
    index_t l0 = {};
    index_t r0 = {};
    _get_left_right_cells_indices<dim>(l0, r0, smPt, graph, m_axis);

    const auto & factors = (bc_type == 1) ? m_factorsReflective : m_factorsDefault;

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

}}}
#endif


// ScalarType coeff = 1.;
// m_J.coeffRef(rowIndex,  col_i)     = (      m_JLneg(0,0) + m_JLpos(0,0) - m_JRneg(0,0))*m_hInv;
// m_J.coeffRef(rowIndex,  col_i+1)   = (coeff*m_JLneg(0,1) + m_JLpos(0,1) - m_JRneg(0,1))*m_hInv;
// m_J.coeffRef(rowIndex,  col_i+2)   = (      m_JLneg(0,2) + m_JLpos(0,2) - m_JRneg(0,2))*m_hInv;
// m_J.coeffRef(rowIndex,  col_ip1)   = (-m_JRpos(0,0))*m_hInv;
// m_J.coeffRef(rowIndex,  col_ip1+1) = (-m_JRpos(0,1))*m_hInv;
// m_J.coeffRef(rowIndex,  col_ip1+2) = (-m_JRpos(0,2))*m_hInv;

// m_J.coeffRef(rowIndex+1, col_i)     = (      m_JLneg(1,0) + m_JLpos(1,0)-m_JRneg(1,0))*m_hInv;
// m_J.coeffRef(rowIndex+1, col_i+1)   = (coeff*m_JLneg(1,1) + m_JLpos(1,1)-m_JRneg(1,1))*m_hInv;
// m_J.coeffRef(rowIndex+1, col_i+2)   = (      m_JLneg(1,2) + m_JLpos(1,2)-m_JRneg(1,2))*m_hInv;
// m_J.coeffRef(rowIndex+1, col_ip1)   = -m_JRpos(1,0)*m_hInv;
// m_J.coeffRef(rowIndex+1, col_ip1+1) = -m_JRpos(1,1)*m_hInv;
// m_J.coeffRef(rowIndex+1, col_ip1+2) = -m_JRpos(1,2)*m_hInv;

// m_J.coeffRef(rowIndex+2, col_i)     = (      m_JLneg(2,0) + m_JLpos(2,0)-m_JRneg(2,0))*m_hInv;
// m_J.coeffRef(rowIndex+2, col_i+1)   = (coeff*m_JLneg(2,1) + m_JLpos(2,1)-m_JRneg(2,1))*m_hInv;
// m_J.coeffRef(rowIndex+2, col_i+2)   = (      m_JLneg(2,2) + m_JLpos(2,2)-m_JRneg(2,2))*m_hInv;
// m_J.coeffRef(rowIndex+2, col_ip1)   = -m_JRpos(2,0)*m_hInv;
// m_J.coeffRef(rowIndex+2, col_ip1+1) = -m_JRpos(2,1)*m_hInv;
// m_J.coeffRef(rowIndex+2, col_ip1+2) = -m_JRpos(2,2)*m_hInv;
