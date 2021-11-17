
#ifndef PRESSIODEMOAPPS_EDGE_CELL_JACOBIAN_HPP_
#define PRESSIODEMOAPPS_EDGE_CELL_JACOBIAN_HPP_

namespace pressiodemoapps{ namespace impl{

namespace{
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

template<int dim, class IndexT, class GraphType>
void _get_left_right_cells_indices(IndexT & l0, IndexT & r0,
				   IndexT & l1, IndexT & r1,
				   const IndexT & smPt,
				   const GraphType & graph,
				   int axis)
{
  _get_left_right_cells_indices<dim>(l0, r0, smPt, graph, axis);

  if (dim == 1){
    l1 = graph(smPt, 3);
    r1 = graph(smPt, 4);
  }
  else if(dim==2){
    throw std::runtime_error("get_left_right_for_2dim_celljacweno MISSING");
  }
  else if(dim==3){
    throw std::runtime_error("get_left_right_for_3dim_celljacweno MISSING");
  }
}

template<int dim, class IndexT, class GraphType>
void _get_left_right_cells_indices(IndexT & l0, IndexT & r0,
				   IndexT & l1, IndexT & r1,
				   IndexT & l2, IndexT & r2,
				   const IndexT & smPt,
				   const GraphType & graph,
				   int axis)
{
  _get_left_right_cells_indices<dim>(l0, r0, l1, r1, smPt, graph, axis);

  if (dim == 1){
    l2 = graph(smPt, 5);
    r2 = graph(smPt, 6);
  }
  else if(dim==2){
    throw std::runtime_error("get_left_right_for_2dim_celljacweno MISSING");
  }
  else if(dim==3){
    throw std::runtime_error("get_left_right_for_3dim_celljacweno MISSING");
  }
}
} // end anonymous namespace


template<
  class ScalarType,
  class SparseMatrix,
  class FluxJacobianType,
  class MeshType,
  class GradType
  >
struct CellJacobianMembers
{
  // used to figure out which direction we are dealing with
  // if ==1, then we are doing "x"
  // if ==2, then we are doing "y"
  // if ==3, then we are doing "z"
  int m_axis = 1;

  InviscidFluxReconstruction m_recEn;
  SparseMatrix & m_J;
  const MeshType & m_meshObj;
  const FluxJacobianType & m_JLNeg;
  const FluxJacobianType & m_JLPos;
  const FluxJacobianType & m_JRNeg;
  const FluxJacobianType & m_JRPos;
  const GradType & m_gradLNeg;
  const GradType & m_gradLPos;
  const GradType & m_gradRNeg;
  const GradType & m_gradRPos;

  // m_hInv = 1/h where h is the cell width
  // depends on which axis we are doing
  ScalarType m_hInv = {};

  CellJacobianMembers(InviscidFluxReconstruction recEn,
			  SparseMatrix & J,
			  const MeshType & meshObj,
			  const FluxJacobianType & JLNeg,
			  const FluxJacobianType & JLPos,
			  const FluxJacobianType & JRNeg,
			  const FluxJacobianType & JRPos,
			  const GradType & gradLNeg,
			  const GradType & gradLPos,
			  const GradType & gradRNeg,
			  const GradType & gradRPos,
			  int axis = 1)
    : m_axis(axis), m_recEn(recEn),
      m_J(J), m_meshObj(meshObj),
      m_JLNeg(JLNeg), m_JLPos(JLPos),
      m_JRNeg(JRNeg), m_JRPos(JRPos),
      m_gradLNeg(gradLNeg), m_gradLPos(gradLPos),
      m_gradRNeg(gradRNeg), m_gradRPos(gradRPos)
  {
    m_hInv = (m_axis == 1) ?
      m_meshObj.dxInv() : (m_axis==2) ?
      m_meshObj.dyInv() : m_meshObj.dzInv();
  }
};

template<
  int dim, int numDofPerCell,
  class ScalarType,
  class SparseMatrix,
  class FluxJacobianType,
  class MeshType,
  class GradType
  >
class InnerCellJacobianFunctor
  : private CellJacobianMembers<
  ScalarType, SparseMatrix, FluxJacobianType, MeshType, GradType
  >
{
  using members = CellJacobianMembers<ScalarType, SparseMatrix,
				      FluxJacobianType, MeshType, GradType>;
  using members::m_axis;
  using members::m_recEn;
  using members::m_J;
  using members::m_meshObj;
  using members::m_JLNeg;
  using members::m_JLPos;
  using members::m_JRNeg;
  using members::m_JRPos;
  using members::m_gradLNeg;
  using members::m_gradLPos;
  using members::m_gradRNeg;
  using members::m_gradRPos;
  using members::m_hInv;

public:
  InnerCellJacobianFunctor(InviscidFluxReconstruction recEn,
			   SparseMatrix & J,
			   const MeshType & meshObj,
			   const FluxJacobianType & JLNeg,
			   const FluxJacobianType & JLPos,
			   const FluxJacobianType & JRNeg,
			   const FluxJacobianType & JRPos,
			   const GradType & gradLNeg,
			   const GradType & gradLPos,
			   const GradType & gradRNeg,
			   const GradType & gradRPos,
			   int axis = 1)
    : members(recEn, J, meshObj,
	      JLNeg, JLPos, JRNeg, JRPos,
	      gradLNeg, gradLPos, gradRNeg, gradRPos, axis)
  {}

  template<class index_t>
  void operator()(const index_t smPt)
  {
    // this class is specialized for numDofPerCell = 1 so we omit it below
    const auto & graph = m_meshObj.graph();
    index_t rowIndex  = smPt*numDofPerCell;

    if (m_recEn == InviscidFluxReconstruction::FirstOrder)
      {
	index_t l0 = {};
	index_t r0 = {};
	_get_left_right_cells_indices<dim>(l0, r0, smPt, graph, m_axis);

	index_t rowIndex  = smPt*numDofPerCell;
	index_t col_im1   = l0*numDofPerCell;
	index_t col_i     = graph(smPt, 0)*numDofPerCell;
	index_t col_ip1   = r0*numDofPerCell;

	for (int k=0; k<numDofPerCell; ++k){
	  for (int j=0; j<numDofPerCell; ++j){
	    m_J.coeffRef(rowIndex+k, col_im1+j) +=  m_JLNeg(k,j)*m_hInv;
	    m_J.coeffRef(rowIndex+k, col_i+j)   += (m_JLPos(k,j)-m_JRNeg(k,j))*m_hInv;
	    m_J.coeffRef(rowIndex+k, col_ip1+j) += -m_JRPos(k,j)*m_hInv;
	  }
	}
      }

    else if (m_recEn == InviscidFluxReconstruction::Weno3)
    {
      index_t l0 = {};
      index_t r0 = {};
      index_t l1 = {};
      index_t r1 = {};
      _get_left_right_cells_indices<dim>(l0, r0, l1, r1,
					 smPt, graph, m_axis);

      index_t col_im2   = l1*numDofPerCell;
      index_t col_im1   = l0*numDofPerCell;
      index_t col_i     = graph(smPt, 0)*numDofPerCell;
      index_t col_ip1   = r0*numDofPerCell;
      index_t col_ip2   = r1*numDofPerCell;

      for (int k=0; k<numDofPerCell; ++k){
	for (int j=0; j<numDofPerCell; ++j){
	  // sensitivy of flux at i-1/2
	  m_J.coeffRef(rowIndex+k, col_im2+j) += m_JLNeg(k,j)*m_gradLNeg(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im2+j) += m_JLPos(k,j)*m_gradLPos(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) += m_JLNeg(k,j)*m_gradLNeg(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) += m_JLPos(k,j)*m_gradLPos(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   += m_JLNeg(k,j)*m_gradLNeg(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   += m_JLPos(k,j)*m_gradLPos(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) += m_JLNeg(k,j)*m_gradLNeg(j,3)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) += m_JLPos(k,j)*m_gradLPos(j,3)*m_hInv;

	  // sensitivy of flux at i+1/2
	  m_J.coeffRef(rowIndex+k, col_im1+j) -= m_JRNeg(k,j)*m_gradRNeg(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) -= m_JRPos(k,j)*m_gradRPos(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   -= m_JRNeg(k,j)*m_gradRNeg(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   -= m_JRPos(k,j)*m_gradRPos(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) -= m_JRNeg(k,j)*m_gradRNeg(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) -= m_JRPos(k,j)*m_gradRPos(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip2+j) -= m_JRNeg(k,j)*m_gradRNeg(j,3)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip2+j) -= m_JRPos(k,j)*m_gradRPos(j,3)*m_hInv;
	}
      }
    }

    else if (m_recEn == InviscidFluxReconstruction::Weno5)
    {
      index_t l0 = {};
      index_t r0 = {};
      index_t l1 = {};
      index_t r1 = {};
      index_t l2 = {};
      index_t r2 = {};
      _get_left_right_cells_indices<dim>(l0, r0, l1, r1, l2, r2,
					 smPt, graph, m_axis);

      index_t col_im3   = l2*numDofPerCell;
      index_t col_im2   = l1*numDofPerCell;
      index_t col_im1   = l0*numDofPerCell;
      index_t col_i     = graph(smPt, 0)*numDofPerCell;
      index_t col_ip1   = r0*numDofPerCell;
      index_t col_ip2   = r1*numDofPerCell;
      index_t col_ip3   = r2*numDofPerCell;

      for (int k=0; k<numDofPerCell; ++k){
	for (int j=0; j<numDofPerCell; ++j){
	  // sensitivy of flux at i-1/2
	  m_J.coeffRef(rowIndex+k, col_im3+j) += m_JLNeg(k,j)*m_gradLNeg(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im3+j) += m_JLPos(k,j)*m_gradLPos(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im2+j) += m_JLNeg(k,j)*m_gradLNeg(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im2+j) += m_JLPos(k,j)*m_gradLPos(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) += m_JLNeg(k,j)*m_gradLNeg(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) += m_JLPos(k,j)*m_gradLPos(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   += m_JLNeg(k,j)*m_gradLNeg(j,3)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   += m_JLPos(k,j)*m_gradLPos(j,3)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) += m_JLNeg(k,j)*m_gradLNeg(j,4)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) += m_JLPos(k,j)*m_gradLPos(j,4)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip2+j) += m_JLNeg(k,j)*m_gradLNeg(j,5)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip2+j) += m_JLPos(k,j)*m_gradLPos(j,5)*m_hInv;

	  // sensitivy of flux at i+1/2
	  m_J.coeffRef(rowIndex+k, col_im2+j) -= m_JRNeg(k,j)*m_gradRNeg(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im2+j) -= m_JRPos(k,j)*m_gradRPos(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) -= m_JRNeg(k,j)*m_gradRNeg(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) -= m_JRPos(k,j)*m_gradRPos(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   -= m_JRNeg(k,j)*m_gradRNeg(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   -= m_JRPos(k,j)*m_gradRPos(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) -= m_JRNeg(k,j)*m_gradRNeg(j,3)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) -= m_JRPos(k,j)*m_gradRPos(j,3)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip2+j) -= m_JRNeg(k,j)*m_gradRNeg(j,4)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip2+j) -= m_JRPos(k,j)*m_gradRPos(j,4)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip3+j) -= m_JRNeg(k,j)*m_gradRNeg(j,5)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip3+j) -= m_JRPos(k,j)*m_gradRPos(j,5)*m_hInv;
	}
      }
    }
  }
};

// ---------------------------------
// partial specialize for 1 dof/cell
// ---------------------------------
template<
  int dim,
  class ScalarType,
  class SparseMatrix,
  class FluxJacobianType,
  class MeshType,
  class GradType
  >
class InnerCellJacobianFunctor<
  dim, 1, ScalarType, SparseMatrix, FluxJacobianType, MeshType, GradType
  >
  : private CellJacobianMembers<
  ScalarType, SparseMatrix, FluxJacobianType, MeshType, GradType
  >
{
  using members = CellJacobianMembers<ScalarType, SparseMatrix,
					  FluxJacobianType, MeshType, GradType>;
  using members::m_axis;
  using members::m_recEn;
  using members::m_J;
  using members::m_meshObj;
  using members::m_JLNeg;
  using members::m_JLPos;
  using members::m_JRNeg;
  using members::m_JRPos;
  using members::m_gradLNeg;
  using members::m_gradLPos;
  using members::m_gradRNeg;
  using members::m_gradRPos;
  using members::m_hInv;

public:
  InnerCellJacobianFunctor(InviscidFluxReconstruction recEn,
			   SparseMatrix & J,
			   const MeshType & meshObj,
			   const FluxJacobianType & JLNeg,
			   const FluxJacobianType & JLPos,
			   const FluxJacobianType & JRNeg,
			   const FluxJacobianType & JRPos,
			   const GradType & gradLNeg,
			   const GradType & gradLPos,
			   const GradType & gradRNeg,
			   const GradType & gradRPos,
			   int axis = 1)
    : members(recEn, J, meshObj, JLNeg, JLPos, JRNeg, JRPos,
	      gradLNeg, gradLPos, gradRNeg, gradRPos, axis)
  {}

  template<class index_t>
  void operator()(const index_t smPt)
  {
    // this class is specialized for numDofPerCell = 1 so we omit it below
    const auto & graph = m_meshObj.graph();
    const index_t rowIndex  = smPt;

    if (m_recEn == InviscidFluxReconstruction::FirstOrder)
      {
	index_t l0 = {};
	index_t r0 = {};
	_get_left_right_cells_indices<dim>(l0, r0, smPt, graph, m_axis);

	index_t rowIndex  = smPt;
	index_t col_im1   = l0;
	index_t col_i     = graph(smPt, 0);
	index_t col_ip1   = r0;

	m_J.coeffRef(rowIndex, col_im1) +=  m_JLNeg*m_hInv;
	m_J.coeffRef(rowIndex, col_i)   += (m_JLPos-m_JRNeg)*m_hInv;
	m_J.coeffRef(rowIndex, col_ip1) += -m_JRPos*m_hInv;
      }

    else if (m_recEn == InviscidFluxReconstruction::Weno3)
    {
      index_t l0 = {};
      index_t r0 = {};
      index_t l1 = {};
      index_t r1 = {};
      _get_left_right_cells_indices<dim>(l0, r0, l1, r1,
					 smPt, graph, m_axis);

      index_t col_im2   = l1;
      index_t col_im1   = l0;
      index_t col_i     = graph(smPt, 0);
      index_t col_ip1   = r0;
      index_t col_ip2   = r1;

      // sensitivy of flux at i-1/2
      m_J.coeffRef(rowIndex, col_im2) += m_JLNeg*m_gradLNeg[0]*m_hInv;
      m_J.coeffRef(rowIndex, col_im2) += m_JLPos*m_gradLPos[0]*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) += m_JLNeg*m_gradLNeg[1]*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) += m_JLPos*m_gradLPos[1]*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   += m_JLNeg*m_gradLNeg[2]*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   += m_JLPos*m_gradLPos[2]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) += m_JLNeg*m_gradLNeg[3]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) += m_JLPos*m_gradLPos[3]*m_hInv;

      // sensitivy of flux at i+1/2
      m_J.coeffRef(rowIndex, col_im1) -= m_JRNeg*m_gradRNeg[0]*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) -= m_JRPos*m_gradRPos[0]*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   -= m_JRNeg*m_gradRNeg[1]*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   -= m_JRPos*m_gradRPos[1]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) -= m_JRNeg*m_gradRNeg[2]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) -= m_JRPos*m_gradRPos[2]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip2) -= m_JRNeg*m_gradRNeg[3]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip2) -= m_JRPos*m_gradRPos[3]*m_hInv;
    }

    else if (m_recEn == InviscidFluxReconstruction::Weno5)
    {
      index_t l0 = {};
      index_t r0 = {};
      index_t l1 = {};
      index_t r1 = {};
      index_t l2 = {};
      index_t r2 = {};
      _get_left_right_cells_indices<dim>(l0, r0, l1, r1, l2, r2,
					 smPt, graph, m_axis);

      index_t col_im3   = l2;
      index_t col_im2   = l1;
      index_t col_im1   = l0;
      index_t col_i     = graph(smPt, 0);
      index_t col_ip1   = r0;
      index_t col_ip2   = r1;
      index_t col_ip3   = r2;

      // sensitivy of flux at i-1/2
      m_J.coeffRef(rowIndex, col_im3) += m_JLNeg*m_gradLNeg[0]*m_hInv;
      m_J.coeffRef(rowIndex, col_im3) += m_JLPos*m_gradLPos[0]*m_hInv;
      m_J.coeffRef(rowIndex, col_im2) += m_JLNeg*m_gradLNeg[1]*m_hInv;
      m_J.coeffRef(rowIndex, col_im2) += m_JLPos*m_gradLPos[1]*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) += m_JLNeg*m_gradLNeg[2]*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) += m_JLPos*m_gradLPos[2]*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   += m_JLNeg*m_gradLNeg[3]*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   += m_JLPos*m_gradLPos[3]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) += m_JLNeg*m_gradLNeg[4]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) += m_JLPos*m_gradLPos[4]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip2) += m_JLNeg*m_gradLNeg[5]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip2) += m_JLPos*m_gradLPos[5]*m_hInv;

      // sensitivy of flux at i+1/2
      m_J.coeffRef(rowIndex, col_im2) -= m_JRNeg*m_gradRNeg[0]*m_hInv;
      m_J.coeffRef(rowIndex, col_im2) -= m_JRPos*m_gradRPos[0]*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) -= m_JRNeg*m_gradRNeg[1]*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) -= m_JRPos*m_gradRPos[1]*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   -= m_JRNeg*m_gradRNeg[2]*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   -= m_JRPos*m_gradRPos[2]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) -= m_JRNeg*m_gradRNeg[3]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) -= m_JRPos*m_gradRPos[3]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip2) -= m_JRNeg*m_gradRNeg[4]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip2) -= m_JRPos*m_gradRPos[4]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip3) -= m_JRNeg*m_gradRNeg[5]*m_hInv;
      m_J.coeffRef(rowIndex, col_ip3) -= m_JRPos*m_gradRPos[5]*m_hInv;
    }
  }
};


template<
  class ScalarType,
  class SparseMatrix,
  class FluxJacobianType,
  class MeshType
  >
struct FirstOrderCellJacobianMembers
{
  // used to figure out which direction we are dealing with
  // if ==1, then we are doing "x"
  // if ==2, then we are doing "y"
  // if ==3, then we are doing "z"
  int m_axis = 1;

  SparseMatrix & m_J;
  const MeshType & m_meshObj;
  const FluxJacobianType & m_JLneg;
  const FluxJacobianType & m_JLpos;
  const FluxJacobianType & m_JRneg;
  const FluxJacobianType & m_JRpos;

  // m_hInv = 1/h where h is the cell width
  // depends on which axis we are doing
  ScalarType m_hInv = {};

  FirstOrderCellJacobianMembers(SparseMatrix & J,
				const MeshType & meshObj,
				const FluxJacobianType & JLneg,
				const FluxJacobianType & JLpos,
				const FluxJacobianType & JRneg,
				const FluxJacobianType & JRpos,
				int axis = 1)
    : m_axis(axis), m_J(J), m_meshObj(meshObj),
      m_JLneg(JLneg), m_JLpos(JLpos), m_JRneg(JRneg), m_JRpos(JRpos)
  {
    m_hInv = (m_axis == 1) ?
      m_meshObj.dxInv() : (m_axis==2) ?
        m_meshObj.dyInv() : m_meshObj.dzInv();
  }
};

template<
  int dim, int numDofPerCell,
  class ScalarType,
  class SparseMatrix,
  class FluxJacobianType,
  class MeshType
  >
class FirstOrderBdCellJacobianFunctor
  : private FirstOrderCellJacobianMembers<ScalarType, SparseMatrix, FluxJacobianType, MeshType>
{
  using members = FirstOrderCellJacobianMembers<ScalarType, SparseMatrix, FluxJacobianType, MeshType>;
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
				  const FluxJacobianType & JLneg,
				  const FluxJacobianType & JLpos,
				  const FluxJacobianType & JRneg,
				  const FluxJacobianType & JRpos,
				  int axis = 1)
    : members(J, meshObj, JLneg, JLpos, JRneg, JRpos, axis){}

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
