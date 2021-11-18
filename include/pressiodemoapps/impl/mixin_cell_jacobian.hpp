
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

template<class Parent, int dim, int ndpc, class MeshType, class JacobianType>
struct InnerCellJacobian : Parent
{
private:
  int m_axis = {};
  typename MeshType::scalar_t m_hInv = {};
  const MeshType & m_meshObj;
  JacobianType & m_J;

public:
  template<class ...Args>
  InnerCellJacobian(JacobianType & J,
		    int axis,
		    const MeshType & meshObj,
		    Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_axis(axis), m_meshObj(meshObj), m_J(J)
  {
    m_hInv = (axis == 1) ? meshObj.dxInv()
      : (axis==2) ? meshObj.dyInv() : meshObj.dzInv();
  }

  template<class index_t>
  void operator()(const index_t smPt)
  {
    Parent::operator()(smPt);

    const auto & graph = m_meshObj.graph();
    const index_t rowIndex  = smPt*ndpc;

    const auto & JLNeg = Parent::fluxJacLNeg();
    const auto & JLPos = Parent::fluxJacLPos();
    const auto & JRNeg = Parent::fluxJacRNeg();
    const auto & JRPos = Parent::fluxJacRPos();

    if (Parent::reconstructionScheme() == ReconstructionScheme::FirstOrder)
      {
	index_t l0 = {};
	index_t r0 = {};
	_get_left_right_cells_indices<dim>(l0, r0, smPt, graph, m_axis);
	index_t rowIndex  = smPt*ndpc;
	index_t col_im1   = l0*ndpc;
	index_t col_i     = graph(smPt, 0)*ndpc;
	index_t col_ip1   = r0*ndpc;
	for (int k=0; k<ndpc; ++k){
	  for (int j=0; j<ndpc; ++j){
	    m_J.coeffRef(rowIndex+k, col_im1+j) +=  JLNeg(k,j)*m_hInv;
	    m_J.coeffRef(rowIndex+k, col_i+j)   += (JLPos(k,j)-JRNeg(k,j))*m_hInv;
	    m_J.coeffRef(rowIndex+k, col_ip1+j) += -JRPos(k,j)*m_hInv;
	  }
	}
      }

    else if (Parent::reconstructionScheme() == ReconstructionScheme::Weno3)
    {
      const auto & JLNeg = Parent::fluxJacLNeg();
      const auto & JLPos = Parent::fluxJacLPos();
      const auto & JRNeg = Parent::fluxJacRNeg();
      const auto & JRPos = Parent::fluxJacRPos();
      const auto & gradLNeg = Parent::reconstructionGradLeftNeg();
      const auto & gradLPos = Parent::reconstructionGradLeftPos();
      const auto & gradRNeg = Parent::reconstructionGradRightNeg();
      const auto & gradRPos = Parent::reconstructionGradRightPos();

      index_t l0 = {};
      index_t r0 = {};
      index_t l1 = {};
      index_t r1 = {};
      _get_left_right_cells_indices<dim>(l0, r0, l1, r1,
					 smPt, graph, m_axis);

      index_t col_im2   = l1*ndpc;
      index_t col_im1   = l0*ndpc;
      index_t col_i     = graph(smPt, 0)*ndpc;
      index_t col_ip1   = r0*ndpc;
      index_t col_ip2   = r1*ndpc;

      for (int k=0; k<ndpc; ++k){
	for (int j=0; j<ndpc; ++j){
	  // sensitivy of flux at i-1/2
	  m_J.coeffRef(rowIndex+k, col_im2+j) += JLNeg(k,j)*gradLNeg(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im2+j) += JLPos(k,j)*gradLPos(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) += JLNeg(k,j)*gradLNeg(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) += JLPos(k,j)*gradLPos(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   += JLNeg(k,j)*gradLNeg(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   += JLPos(k,j)*gradLPos(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) += JLNeg(k,j)*gradLNeg(j,3)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) += JLPos(k,j)*gradLPos(j,3)*m_hInv;

	  // sensitivy of flux at i+1/2
	  m_J.coeffRef(rowIndex+k, col_im1+j) -= JRNeg(k,j)*gradRNeg(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) -= JRPos(k,j)*gradRPos(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   -= JRNeg(k,j)*gradRNeg(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   -= JRPos(k,j)*gradRPos(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) -= JRNeg(k,j)*gradRNeg(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) -= JRPos(k,j)*gradRPos(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip2+j) -= JRNeg(k,j)*gradRNeg(j,3)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip2+j) -= JRPos(k,j)*gradRPos(j,3)*m_hInv;
	}
      }
    }

    else if (Parent::reconstructionScheme() == ReconstructionScheme::Weno5)
    {

      const auto & JLNeg = Parent::fluxJacLNeg();
      const auto & JLPos = Parent::fluxJacLPos();
      const auto & JRNeg = Parent::fluxJacRNeg();
      const auto & JRPos = Parent::fluxJacRPos();
      const auto & gradLNeg = Parent::reconstructionGradLeftNeg();
      const auto & gradLPos = Parent::reconstructionGradLeftPos();
      const auto & gradRNeg = Parent::reconstructionGradRightNeg();
      const auto & gradRPos = Parent::reconstructionGradRightPos();

      index_t l0 = {};
      index_t r0 = {};
      index_t l1 = {};
      index_t r1 = {};
      index_t l2 = {};
      index_t r2 = {};
      _get_left_right_cells_indices<dim>(l0, r0, l1, r1, l2, r2,
					 smPt, graph, m_axis);

      index_t col_im3   = l2*ndpc;
      index_t col_im2   = l1*ndpc;
      index_t col_im1   = l0*ndpc;
      index_t col_i     = graph(smPt, 0)*ndpc;
      index_t col_ip1   = r0*ndpc;
      index_t col_ip2   = r1*ndpc;
      index_t col_ip3   = r2*ndpc;

      for (int k=0; k<ndpc; ++k){
	for (int j=0; j<ndpc; ++j){
	  // sensitivy of flux at i-1/2
	  m_J.coeffRef(rowIndex+k, col_im3+j) += JLNeg(k,j)*gradLNeg(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im3+j) += JLPos(k,j)*gradLPos(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im2+j) += JLNeg(k,j)*gradLNeg(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im2+j) += JLPos(k,j)*gradLPos(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) += JLNeg(k,j)*gradLNeg(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) += JLPos(k,j)*gradLPos(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   += JLNeg(k,j)*gradLNeg(j,3)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   += JLPos(k,j)*gradLPos(j,3)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) += JLNeg(k,j)*gradLNeg(j,4)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) += JLPos(k,j)*gradLPos(j,4)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip2+j) += JLNeg(k,j)*gradLNeg(j,5)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip2+j) += JLPos(k,j)*gradLPos(j,5)*m_hInv;

	  // sensitivy of flux at i+1/2
	  m_J.coeffRef(rowIndex+k, col_im2+j) -= JRNeg(k,j)*gradRNeg(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im2+j) -= JRPos(k,j)*gradRPos(j,0)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) -= JRNeg(k,j)*gradRNeg(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_im1+j) -= JRPos(k,j)*gradRPos(j,1)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   -= JRNeg(k,j)*gradRNeg(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_i+j)   -= JRPos(k,j)*gradRPos(j,2)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) -= JRNeg(k,j)*gradRNeg(j,3)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip1+j) -= JRPos(k,j)*gradRPos(j,3)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip2+j) -= JRNeg(k,j)*gradRNeg(j,4)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip2+j) -= JRPos(k,j)*gradRPos(j,4)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip3+j) -= JRNeg(k,j)*gradRNeg(j,5)*m_hInv;
	  m_J.coeffRef(rowIndex+k, col_ip3+j) -= JRPos(k,j)*gradRPos(j,5)*m_hInv;
	}
      }
    }

  }// end operator()
};


template<class Parent, int dim, class MeshType, class JacobianType>
struct InnerCellJacobian<Parent, dim, 1, MeshType, JacobianType> : Parent
{
private:
  int m_axis = {};
  typename MeshType::scalar_t m_hInv = {};
  const MeshType & m_meshObj;
  JacobianType & m_J;

public:
  template<class ...Args>
  InnerCellJacobian(JacobianType & J,
		    int axis,
		    const MeshType & meshObj,
		    Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_axis(axis), m_meshObj(meshObj), m_J(J)
  {
    m_hInv = (axis == 1) ? meshObj.dxInv()
      : (axis==2) ? meshObj.dyInv() : meshObj.dzInv();
  }

  template<class index_t>
  void operator()(const index_t smPt)
  {
    Parent::operator()(smPt);

    const auto & graph = m_meshObj.graph();
    const index_t rowIndex  = smPt;

    const auto & JLNeg = Parent::fluxJacLNeg();
    const auto & JLPos = Parent::fluxJacLPos();
    const auto & JRNeg = Parent::fluxJacRNeg();
    const auto & JRPos = Parent::fluxJacRPos();

    if (Parent::reconstructionScheme() == ReconstructionScheme::FirstOrder)
      {
	index_t l0 = {};
	index_t r0 = {};
	_get_left_right_cells_indices<dim>(l0, r0, smPt, graph, m_axis);
	index_t rowIndex  = smPt;
	index_t col_im1   = l0;
	index_t col_i     = graph(smPt, 0);
	index_t col_ip1   = r0;

	m_J.coeffRef(rowIndex, col_im1) +=  JLNeg*m_hInv;
	m_J.coeffRef(rowIndex, col_i)   += (JLPos-JRNeg)*m_hInv;
	m_J.coeffRef(rowIndex, col_ip1) += -JRPos*m_hInv;
      }

    else if (Parent::reconstructionScheme() == ReconstructionScheme::Weno3)
    {
      const auto & JLNeg = Parent::fluxJacLNeg();
      const auto & JLPos = Parent::fluxJacLPos();
      const auto & JRNeg = Parent::fluxJacRNeg();
      const auto & JRPos = Parent::fluxJacRPos();
      const auto & gradLNeg = Parent::reconstructionGradLeftNeg();
      const auto & gradLPos = Parent::reconstructionGradLeftPos();
      const auto & gradRNeg = Parent::reconstructionGradRightNeg();
      const auto & gradRPos = Parent::reconstructionGradRightPos();

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
      m_J.coeffRef(rowIndex, col_im2) += JLNeg*gradLNeg(0)*m_hInv;
      m_J.coeffRef(rowIndex, col_im2) += JLPos*gradLPos(0)*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) += JLNeg*gradLNeg(1)*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) += JLPos*gradLPos(1)*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   += JLNeg*gradLNeg(2)*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   += JLPos*gradLPos(2)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) += JLNeg*gradLNeg(3)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) += JLPos*gradLPos(3)*m_hInv;

      // sensitivy of flux at i+1/2
      m_J.coeffRef(rowIndex, col_im1) -= JRNeg*gradRNeg(0)*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) -= JRPos*gradRPos(0)*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   -= JRNeg*gradRNeg(1)*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   -= JRPos*gradRPos(1)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) -= JRNeg*gradRNeg(2)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) -= JRPos*gradRPos(2)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip2) -= JRNeg*gradRNeg(3)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip2) -= JRPos*gradRPos(3)*m_hInv;
    }

    else if (Parent::reconstructionScheme() == ReconstructionScheme::Weno5)
    {
      const auto & JLNeg = Parent::fluxJacLNeg();
      const auto & JLPos = Parent::fluxJacLPos();
      const auto & JRNeg = Parent::fluxJacRNeg();
      const auto & JRPos = Parent::fluxJacRPos();
      const auto & gradLNeg = Parent::reconstructionGradLeftNeg();
      const auto & gradLPos = Parent::reconstructionGradLeftPos();
      const auto & gradRNeg = Parent::reconstructionGradRightNeg();
      const auto & gradRPos = Parent::reconstructionGradRightPos();

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
      m_J.coeffRef(rowIndex, col_im3) += JLNeg*gradLNeg(0)*m_hInv;
      m_J.coeffRef(rowIndex, col_im3) += JLPos*gradLPos(0)*m_hInv;
      m_J.coeffRef(rowIndex, col_im2) += JLNeg*gradLNeg(1)*m_hInv;
      m_J.coeffRef(rowIndex, col_im2) += JLPos*gradLPos(1)*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) += JLNeg*gradLNeg(2)*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) += JLPos*gradLPos(2)*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   += JLNeg*gradLNeg(3)*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   += JLPos*gradLPos(3)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) += JLNeg*gradLNeg(4)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) += JLPos*gradLPos(4)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip2) += JLNeg*gradLNeg(5)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip2) += JLPos*gradLPos(5)*m_hInv;

      // sensitivy of flux at i+1/2
      m_J.coeffRef(rowIndex, col_im2) -= JRNeg*gradRNeg(0)*m_hInv;
      m_J.coeffRef(rowIndex, col_im2) -= JRPos*gradRPos(0)*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) -= JRNeg*gradRNeg(1)*m_hInv;
      m_J.coeffRef(rowIndex, col_im1) -= JRPos*gradRPos(1)*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   -= JRNeg*gradRNeg(2)*m_hInv;
      m_J.coeffRef(rowIndex, col_i)   -= JRPos*gradRPos(2)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) -= JRNeg*gradRNeg(3)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip1) -= JRPos*gradRPos(3)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip2) -= JRNeg*gradRNeg(4)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip2) -= JRPos*gradRPos(4)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip3) -= JRNeg*gradRNeg(5)*m_hInv;
      m_J.coeffRef(rowIndex, col_ip3) -= JRPos*gradRPos(5)*m_hInv;
    }

  }// end operator()
};



template<class Parent, int dim, int ndpc, class MeshType, class JacobianType>
struct FirstOrderBdCellJacobian : Parent
{
private:
  int m_axis = {};
  typename MeshType::scalar_t m_hInv = {};
  const MeshType & m_meshObj;
  JacobianType & m_J;

  using scalar_type = typename MeshType::scalar_t;

public:
  template<class ...Args>
  FirstOrderBdCellJacobian(JacobianType & J,
			   int axis,
			   const MeshType & meshObj,
			   Args && ...args)
    : Parent(std::forward<Args>(args)...),
      m_axis(axis), m_meshObj(meshObj), m_J(J)
  {
    m_hInv = (axis == 1) ? meshObj.dxInv()
      : (axis==2) ? meshObj.dyInv() : meshObj.dzInv();
  }

  template<class index_t>
  void operator()(const index_t smPt,
		  const std::array<scalar_type, ndpc> & factors,
		  int bc_type)
  {
    Parent::operator()(smPt);

    const auto & JLNeg = Parent::fluxJacLNeg();
    const auto & JLPos = Parent::fluxJacLPos();
    const auto & JRNeg = Parent::fluxJacRNeg();
    const auto & JRPos = Parent::fluxJacRPos();

    const auto & graph = m_meshObj.graph();
    index_t l0 = {};
    index_t r0 = {};
    _get_left_right_cells_indices<dim>(l0, r0, smPt, graph, m_axis);

    auto rowIndex  = smPt*ndpc;
    auto col_i     = graph(smPt, 0)*ndpc;

    for (int k=0; k<ndpc; ++k){
      for (int j=0; j<ndpc; ++j){
	m_J.coeffRef(rowIndex+k, col_i+j) += (JLPos(k,j) - JRNeg(k,j))*m_hInv;
      }
    }

    if (l0 != -1){
      auto col_im1 = l0*ndpc;
      for (int k=0; k<ndpc; ++k){
	for (int j=0; j<ndpc; ++j){
	  m_J.coeffRef(rowIndex+k, col_im1+j) += JLNeg(k,j)*m_hInv;
	}
      }
    }

    if (r0 != -1){
      auto col_ip1 = r0*ndpc;
      for (int k=0; k<ndpc; ++k){
	for (int j=0; j<ndpc; ++j){
	  m_J.coeffRef(rowIndex+k, col_ip1+j) += -JRPos(k,j)*m_hInv;
	}
      }
    }

    if (bc_type != 2)
    {
      if (l0 == -1){
	for (int k=0; k<ndpc; ++k){
	  for (int j=0; j<ndpc; ++j){
	    m_J.coeffRef(rowIndex+k, col_i+j) += (factors[j]*JLNeg(k,j))*m_hInv;
	  }
	}
      }

      if (r0 == -1){
	for (int k=0; k<ndpc; ++k){
	  for (int j=0; j<ndpc; ++j){
	    m_J.coeffRef(rowIndex+k, col_i+j) += (factors[j]*-JRPos(k,j))*m_hInv;
	  }
	}
      }
    }

  }// end operator()
};

}}
#endif


// // ---------------------------------
// // partial specialize for 1 dof/cell
// // ---------------------------------
