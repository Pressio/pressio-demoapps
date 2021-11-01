
#ifndef PRESSIODEMOAPPS_EDGE_CELL_JACOBIAN_HPP_
#define PRESSIODEMOAPPS_EDGE_CELL_JACOBIAN_HPP_

namespace pressiodemoapps{ namespace impl{

template<int dim, class ScalarType, class MatType, class T, class MeshType>
class FirstOrderNearBDCellJacobianThreeDofFunctor
{
  int m_axis = 1;
  SparseMatrix m_J;
  T & m_JLneg;
  T & m_JLpos;
  T & m_JRneg;
  T & m_JRpos;

public:
  FirstOrderNearBDCellJacobianThreeDofFunctor(int axis,
					      SparseMatrix J,
					      T & JLneg, T & JLpos, T & JRneg, T & JRpos,
					      bool is_reflective)
  {}

  FirstOrderNearBDCellJacobianThreeDofFunctor(SparseMatrix J,
					      T & JLneg, T & JLpos, T & JRneg, T & JRpos,
					      bool is_reflective)
  {}

  template<class index_t>
  void operator()(const index_t smPt)
  {
    auto cellGIDStencilMesh = graph(smPt, 0);
    auto l0 = {};
    auto r0 = {};
    if (dim == 1){
      l0 = graph(smPt, 1);
      r0 = graph(smPt, 2);
    }
    else if(dim==2){
      l0  = (axis == 1) ? graph(smPt, 1) : graph(smPt, 4);
      r0  = (axis == 1) ? graph(smPt, 3) : graph(smPt, 2);
    }
    else if(dim==3){
      l0 = (axis == 1) ? graph(smPt, 1) : (axis==2) ? graph(smPt, 4) : graph(smPt, 5);
      r0 = (axis == 1) ? graph(smPt, 3) : (axis==2) ? graph(smPt, 2) : graph(smPt, 6);
    }

    if (thisCellHasBoundaryOnLeft)
    {
      const auto coeff = is_reflective ? -1. : 1.;

      rowIndex = smPt*numDofPerCell;
      auto c_im1 = l0*numDofPerCell;
      auto c_ip1 = r0*numDofPerCell;

      auto row = m_jacobian.row(rowIndex);
      // fill rows
      //,  c_i,     JLneg(0,0) + JLpos(0,0) - JRneg(0,0));

      m_jacobian(rowIndex,  c_i+1,   coeff*JLneg(0,1) + JLpos(0,1) - JRneg(0,0));
      m_jacobian(rowIndex,  c_i+2,   JLneg(0,2) + JLpos(0,2) - JRneg(0,0));
      m_jacobian(rowIndex,  c_ip1,   -JRpos(0,0));
      m_jacobian(rowIndex,  c_ip1+1, -JRpos(0,1));
      m_jacobian(rowIndex,  c_ip1+2, -JRpos(0,2));

      m_jacobian(rowIndex+1, c_i,	    JLneg(1,0) + JLpos(1,0)-JRneg(1,0));
      m_jacobian(rowIndex+1, c_i+1,   coeff*JLneg(1,1) + JLpos(1,1)-JRneg(1,0));
      m_jacobian(rowIndex+1, c_i+2,   JLneg(1,2) + JLpos(1,2)-JRneg(1,0));
      m_jacobian(rowIndex+1, c_ip1,   -JRpos(1,0));
      m_jacobian(rowIndex+1, c_ip1+1, -JRpos(1,1));
      m_jacobian(rowIndex+1, c_ip1+2, -JRpos(1,2));

      m_jacobian(rowIndex+2, c_i,     JLneg(2,0) + JLpos(2,0)-JRneg(2,0));
      m_jacobian(rowIndex+2, c_i+1,   coeff*JLneg(2,1) + JLpos(2,1)-JRneg(2,0));
      m_jacobian(rowIndex+2, c_i+2,   JLneg(2,2) + JLpos(2,2)-JRneg(2,0));
      m_jacobian(rowIndex+2, c_ip1,   -JRpos(2,0));
      m_jacobian(rowIndex+2, c_ip1+1, -JRpos(2,1));
      m_jacobian(rowIndex+2, c_ip1+2, -JRpos(2,2));
    }

    if (thisCellHasCoundaryOnRight)
      {
	m_jacobian(rowIndex,  c_im1,   JLneg(0,0));
	m_jacobian(rowIndex,  c_im1+1, JLneg(0,1));
	m_jacobian(rowIndex,  c_im1+2, JLneg(0,2));
	m_jacobian(rowIndex,  c_i,     -JRpos(0,0) + JLpos(0,0) - JRneg(0,0));
	m_jacobian(rowIndex,  c_i+1,   coeff*-JRpos(0,1) + JLpos(0,1) - JRneg(0,0));
	m_jacobian(rowIndex,  c_i+2,   -JRpos(0,2) + JLpos(0,2) - JRneg(0,0));

	m_jacobian(rowIndex+1, c_im1,   JLneg(1,0));
	m_jacobian(rowIndex+1, c_im1+1, JLneg(1,1));
	m_jacobian(rowIndex+1, c_im1+2, JLneg(1,2));
	m_jacobian(rowIndex+1, c_i,     -JRpos(1,0) + JLpos(1,0)-JRneg(1,0));
	m_jacobian(rowIndex+1, c_i+1,   coeff*-JRpos(1,1) + JLpos(1,1)-JRneg(1,0));
	m_jacobian(rowIndex+1, c_i+2,   -JRpos(1,2) + JLpos(1,2)-JRneg(1,0));

	m_jacobian(rowIndex+2, c_im1,   JLneg(2,0));
	m_jacobian(rowIndex+2, c_im1+1, JLneg(2,1));
	m_jacobian(rowIndex+2, c_im1+2, JLneg(2,2));
	m_jacobian(rowIndex+2, c_i,     -JRpos(2,0) + JLpos(2,0)-JRneg(2,0));
	m_jacobian(rowIndex+2, c_i+1,   coeff*-JRpos(2,1) + JLpos(2,1)-JRneg(2,0));
	m_jacobian(rowIndex+2, c_i+2,   -JRpos(2,2) + JLpos(2,2)-JRneg(2,0));
    }
  }
};

template<class ScalarType, class T>
class CellJacobianThreeDofPerCellFirstOrderFunctor
{
  SparseMatrix m_J;
  T & m_JLneg;
  T & m_JLpos;
  T & m_JRneg;
  T & m_JRpos;

public:
  CellJacobianThreeDofPerCellFirstOrderFunctor() = delete;
  CellJacobianThreeDofPerCellFirstOrderFunctor(::pressiodemoapps::InviscidFluxReconstruction recEn,
				      SparseMatrix J,
				      T & JLneg, T & JLpos, T & JRneg, T & JRpos)
  {}

  template<class index_t>
  void operator()(const index_t smPt)
  {
    const auto cellGIDStencilMesh      = graph(smPt, 0);
    // const auto leftCellGIDStencilMesh  = graph(smPt, 1);
    // const auto rightCellGIDStencilMesh = graph(smPt, 2);

    rowIndex = smPt*numDofPerCell;
    auto c_im1 = leftCellGIDStencilMesh*numDofPerCell;
    auto c_ip1 = rightCellGIDStencilMesh*numDofPerCell;

    m_jac.insert(rowIndex,  c_im1,   JLneg(0,0));
    m_jac.insert(rowIndex,  c_im1+1, JLneg(0,1));
    m_jac.insert(rowIndex,  c_im1+2, JLneg(0,2));
    m_jac.insert(rowIndex,  c_i,     JLpos(0,0)-JRneg(0,0));
    m_jac.insert(rowIndex,  c_i+1,   JLpos(0,1)-JRneg(0,0));
    m_jac.insert(rowIndex,  c_i+2,   JLpos(0,2)-JRneg(0,0));
    m_jac.insert(rowIndex,  c_ip1,   -JRpos(0,0));
    m_jac.insert(rowIndex,  c_ip1+1, -JRpos(0,1));
    m_jac.insert(rowIndex,  c_ip1+2, -JRpos(0,2));

    m_jac.insert(rowIndex+1, c_im1,   JLneg(1,0));
    m_jac.insert(rowIndex+1, c_im1+1, JLneg(1,1));
    m_jac.insert(rowIndex+1, c_im1+2, JLneg(1,2));
    m_jac.insert(rowIndex+1, c_i,   JLpos(1,0)-JRneg(1,0));
    m_jac.insert(rowIndex+1, c_i+1, JLpos(1,1)-JRneg(1,0));
    m_jac.insert(rowIndex+1, c_i+2, JLpos(1,2)-JRneg(1,0));
    m_jac.insert(rowIndex+1, c_ip1,   -JRpos(1,0));
    m_jac.insert(rowIndex+1, c_ip1+1, -JRpos(1,1));
    m_jac.insert(rowIndex+1, c_ip1+2, -JRpos(1,2));

    m_jac.insert(rowIndex+2, c_im1,   JLneg(2,0));
    m_jac.insert(rowIndex+2, c_im1+1, JLneg(2,1));
    m_jac.insert(rowIndex+2, c_im1+2, JLneg(2,2));
    m_jac.insert(rowIndex+2, c_i,   JLpos(2,0)-JRneg(2,0));
    m_jac.insert(rowIndex+2, c_i+1, JLpos(2,1)-JRneg(2,0));
    m_jac.insert(rowIndex+2, c_i+2, JLpos(2,2)-JRneg(2,0));
    m_jac.insert(rowIndex+2, c_ip1,   -JRpos(2,0));
    m_jac.insert(rowIndex+2, c_ip1+1, -JRpos(2,1));
    m_jac.insert(rowIndex+2, c_ip1+2, -JRpos(2,2));
  }
};



template<class ScalarType, class T>
class EigenJacobianThreeDofPerCellFunctorInnerCell
{
  T & m_JLneg;
  T & m_JLpos;
  T & m_JRneg;
  T & m_JRpos;

public:
  EigenJacobianThreeDofPerCellFunctorInnerCell() = delete;
  EigenJacobianThreeDofPerCellFunctorInnerCell(::pressiodemoapps::InviscidFluxReconstruction recEn,
				      T & JLneg, T & JLpos, T & JRneg, T & JRpos)
  {}

  template<class index_t>
  void operator()(const index_t smPt)
  {
    const auto cellGIDStencilMesh      = graph(smPt, 0);
    const auto leftCellGIDStencilMesh  = graph(smPt, 1);
    const auto rightCellGIDStencilMesh = graph(smPt, 2);

    Matrix<double,-1,-1> q(3, u_stencilVals.size());
    // q(:, 0) is d_uLneg/d_ui-2
    // q(:, 1) is d_uLneg/d_ui-1
    // q(:, 2) is d_uLneg/d_ui
    // q(:, 3) is d_uLpos/d_ui-1
    // q(:, 4) is d_uLpos/d_ui
    // q(:, 5) is d_uLpos/d_ui+1

    // q(:, 6)  is d_uRneg/d_ui-1
    // q(:, 7)  is d_uRneg/d_ui
    // q(:, 8)  is d_uRneg/d_ui+1
    // q(:, 9)  is d_uRpos/d_ui
    // q(:, 10) is d_uRpos/d_ui+1
    // q(:, 11) is d_uRpos/d_ui+2

    const auto w0 = graph(smPt, 1);
    const auto e0 = graph(smPt, 2);
    const auto w1 = graph(smPt, 3);
    const auto e1 = graph(smPt, 4);

    const auto num_pertur = u_stencilVals.size()/numDofPerCell;
    for (int j=0; j<num_perturb; ++j)
      {
	svals2 = m_stencilVals;
	// if j=0, perturb first numDofPerCell of stencilVals2
	// if j=1, perturb second numDofPerCell of stencilVals2
	// ...

	// assuming that we have svals2 that contains the perturb values
	// we need to run the reconstruction
	Reconstructor2.template operator()<numDofPerCell>();
	// here, it means that, these are overwritten:
	// uMinusHalfNegFD
	// uMinusHalfPosFD
	// uPlusHalfNegFD
	// uPlusHalfPosFD

	if (j==0){
	  q(0,0) = (uMinusHalfNegFD(0)-uMinusHalfNeg(0))/eps;
	  q(1,0) = (uMinusHalfNegFD(1)-uMinusHalfNeg(1))/eps;
	  q(2,0) = (uMinusHalfNegFD(2)-uMinusHalfNeg(2))/eps;

	  auto M = JLneg * q.col(0).asDiagonalMatrix();
	  m_tripletList.push_back( Tr(vIndexDensityDof,   w1*numDofPerCell,   M(0,0)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof,   w1*numDofPerCell+1, M(0,1)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof,   w1*numDofPerCell+2, M(0,2)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof+1, w1*numDofPerCell,   M(1,0)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof+1, w1*numDofPerCell+1, M(1,1)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof+1, w1*numDofPerCell+2, M(1,2)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof+2, w1*numDofPerCell,   M(2,0)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof+2, w1*numDofPerCell+1, M(2,1)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof+2, w1*numDofPerCell+2, M(2,2)) );
	}

	if (j==1){
	  q(0,1) = (uMinusHalfNegFD(0)-uMinusHalfNeg(0))/eps;
	  q(1,1) = (uMinusHalfNegFD(1)-uMinusHalfNeg(1))/eps;
	  q(2,1) = (uMinusHalfNegFD(2)-uMinusHalfNeg(2))/eps;
	  auto M1 = JLneg * q.col(1).asDiagonalMatrix();

	  q(0,3) = (uMinusHalfPosFD(0)-uMinusHalfPos(0))/eps;
	  q(1,3) = (uMinusHalfPosFD(1)-uMinusHalfPos(1))/eps;
	  q(2,3) = (uMinusHalfPosFD(2)-uMinusHalfPos(2))/eps;
	  auto M2 = JLpos * q.col(3).asDiagonalMatrix();

	  q(0,6) = (uPlusHalfNegFD(0)-uPlusHalfNeg(0))/eps;
	  q(1,6) = (uPlusHalfNegFD(1)-uPlusHalfNeg(1))/eps;
	  q(2,6) = (uPlusHalfNegFD(2)-uPlusHalfNeg(2))/eps;
	  auto M3 = JRneg * q.col(6).asDiagonalMatrix();

	  auto M = M1+M2+M3;
	  m_tripletList.push_back( Tr(vIndexDensityDof,   w0*numDofPerCell,   M(0,0)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof,   w0*numDofPerCell+1, M(0,1)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof,   w0*numDofPerCell+2, M(0,2)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof+1, w0*numDofPerCell,   M(1,0)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof+1, w0*numDofPerCell+1, M(1,1)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof+1, w0*numDofPerCell+2, M(1,2)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof+2, w0*numDofPerCell,   M(2,0)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof+2, w0*numDofPerCell+1, M(2,1)) );
	  m_tripletList.push_back( Tr(vIndexDensityDof+2, w0*numDofPerCell+2, M(2,2)) );
	}

	if (j==2){
	  q(0,2) = (uMinusHalfNegFD(0)-uMinusHalfNeg(0))/eps;
	  q(1,2) = (uMinusHalfNegFD(1)-uMinusHalfNeg(1))/eps;
	  q(2,2) = (uMinusHalfNegFD(2)-uMinusHalfNeg(2))/eps;

	  q(0,4) = (uMinusHalfPosFD(0)-uMinusHalfPos(0))/eps;
	  q(1,4) = (uMinusHalfPosFD(1)-uMinusHalfPos(1))/eps;
	  q(2,4) = (uMinusHalfPosFD(2)-uMinusHalfPos(2))/eps;

	  q(0,7) = (uPlusHalfNegFD(0)-uPlusHalfNeg(0))/eps;
	  q(1,7) = (uPlusHalfNegFD(1)-uPlusHalfNeg(1))/eps;
	  q(2,7) = (uPlusHalfNegFD(2)-uPlusHalfNeg(2))/eps;

	  q(0,9) = (uPlusHalfPosFD(0)-uPlusHalfPos(0))/eps;
	  q(1,9) = (uPlusHalfPosFD(1)-uPlusHalfPos(1))/eps;
	  q(2,9) = (uPlusHalfPosFD(2)-uPlusHalfPos(2))/eps;
	}

	if (j==3){
	  q(0,5) = (uMinusHalfPosFD(0)-uMinusHalfPos(0))/eps;
	  q(1,5) = (uMinusHalfPosFD(1)-uMinusHalfPos(1))/eps;
	  q(2,5) = (uMinusHalfPosFD(2)-uMinusHalfPos(2))/eps;

	  q(0,8) = (uPlusHalfNegFD(0)-uPlusHalfNeg(0))/eps;
	  q(1,8) = (uPlusHalfNegFD(1)-uPlusHalfNeg(1))/eps;
	  q(2,8) = (uPlusHalfNegFD(2)-uPlusHalfNeg(2))/eps;

	  q(0,10) = (uPlusHalfPosFD(0)-uPlusHalfPos(0))/eps;
	  q(1,10) = (uPlusHalfPosFD(1)-uPlusHalfPos(1))/eps;
	  q(2,10) = (uPlusHalfPosFD(2)-uPlusHalfPos(2))/eps;
	}

	if (j==4){
	  q(0,11) = (uPlusHalfPosFD(0)-uPlusHalfPos(0))/eps;
	  q(1,11) = (uPlusHalfPosFD(1)-uPlusHalfPos(1))/eps;
	  q(2,11) = (uPlusHalfPosFD(2)-uPlusHalfPos(2))/eps;
	}
      }
  }
};

}}
#endif
