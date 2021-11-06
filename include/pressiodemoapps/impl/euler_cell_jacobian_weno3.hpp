
#ifndef PRESSIODEMOAPPS_EDGE_CELL_JACOBIAN_WENO3_HPP_
#define PRESSIODEMOAPPS_EDGE_CELL_JACOBIAN_WENO3_HPP_

namespace pressiodemoapps{ namespace ee{ namespace impl{

// template<class ScalarType, class T>
// class EigenJacobianThreeDofPerCellFunctorInnerCell
// {
//   T & m_JLneg;
//   T & m_JLpos;
//   T & m_JRneg;
//   T & m_JRpos;

// public:
//   EigenJacobianThreeDofPerCellFunctorInnerCell() = delete;
//   EigenJacobianThreeDofPerCellFunctorInnerCell(::pressiodemoapps::InviscidFluxReconstruction recEn,
// 				      T & JLneg, T & JLpos, T & JRneg, T & JRpos)
//   {}

//   template<class index_t>
//   void operator()(const index_t smPt)
//   {
//     const auto cellGIDStencilMesh      = graph(smPt, 0);
//     const auto leftCellGIDStencilMesh  = graph(smPt, 1);
//     const auto rightCellGIDStencilMesh = graph(smPt, 2);

//     Matrix<double,-1,-1> q(3, u_stencilVals.size());
//     // q(:, 0) is d_uLneg/d_ui-2
//     // q(:, 1) is d_uLneg/d_ui-1
//     // q(:, 2) is d_uLneg/d_ui
//     // q(:, 3) is d_uLpos/d_ui-1
//     // q(:, 4) is d_uLpos/d_ui
//     // q(:, 5) is d_uLpos/d_ui+1

//     // q(:, 6)  is d_uRneg/d_ui-1
//     // q(:, 7)  is d_uRneg/d_ui
//     // q(:, 8)  is d_uRneg/d_ui+1
//     // q(:, 9)  is d_uRpos/d_ui
//     // q(:, 10) is d_uRpos/d_ui+1
//     // q(:, 11) is d_uRpos/d_ui+2

//     const auto w0 = graph(smPt, 1);
//     const auto e0 = graph(smPt, 2);
//     const auto w1 = graph(smPt, 3);
//     const auto e1 = graph(smPt, 4);

//     const auto num_pertur = u_stencilVals.size()/numDofPerCell;
//     for (int j=0; j<num_perturb; ++j)
//       {
// 	svals2 = m_stencilVals;
// 	// if j=0, perturb first numDofPerCell of stencilVals2
// 	// if j=1, perturb second numDofPerCell of stencilVals2
// 	// ...

// 	// assuming that we have svals2 that contains the perturb values
// 	// we need to run the reconstruction
// 	Reconstructor2.template operator()<numDofPerCell>();
// 	// here, it means that, these are overwritten:
// 	// uMinusHalfNegFD
// 	// uMinusHalfPosFD
// 	// uPlusHalfNegFD
// 	// uPlusHalfPosFD

// 	if (j==0){
// 	  q(0,0) = (uMinusHalfNegFD(0)-uMinusHalfNeg(0))/eps;
// 	  q(1,0) = (uMinusHalfNegFD(1)-uMinusHalfNeg(1))/eps;
// 	  q(2,0) = (uMinusHalfNegFD(2)-uMinusHalfNeg(2))/eps;

// 	  auto M = JLneg * q.col(0).asDiagonalMatrix();
// 	  m_tripletList.push_back( Tr(vIndexDensityDof,   w1*numDofPerCell,   M(0,0)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof,   w1*numDofPerCell+1, M(0,1)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof,   w1*numDofPerCell+2, M(0,2)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof+1, w1*numDofPerCell,   M(1,0)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof+1, w1*numDofPerCell+1, M(1,1)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof+1, w1*numDofPerCell+2, M(1,2)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof+2, w1*numDofPerCell,   M(2,0)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof+2, w1*numDofPerCell+1, M(2,1)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof+2, w1*numDofPerCell+2, M(2,2)) );
// 	}

// 	if (j==1){
// 	  q(0,1) = (uMinusHalfNegFD(0)-uMinusHalfNeg(0))/eps;
// 	  q(1,1) = (uMinusHalfNegFD(1)-uMinusHalfNeg(1))/eps;
// 	  q(2,1) = (uMinusHalfNegFD(2)-uMinusHalfNeg(2))/eps;
// 	  auto M1 = JLneg * q.col(1).asDiagonalMatrix();

// 	  q(0,3) = (uMinusHalfPosFD(0)-uMinusHalfPos(0))/eps;
// 	  q(1,3) = (uMinusHalfPosFD(1)-uMinusHalfPos(1))/eps;
// 	  q(2,3) = (uMinusHalfPosFD(2)-uMinusHalfPos(2))/eps;
// 	  auto M2 = JLpos * q.col(3).asDiagonalMatrix();

// 	  q(0,6) = (uPlusHalfNegFD(0)-uPlusHalfNeg(0))/eps;
// 	  q(1,6) = (uPlusHalfNegFD(1)-uPlusHalfNeg(1))/eps;
// 	  q(2,6) = (uPlusHalfNegFD(2)-uPlusHalfNeg(2))/eps;
// 	  auto M3 = JRneg * q.col(6).asDiagonalMatrix();

// 	  auto M = M1+M2+M3;
// 	  m_tripletList.push_back( Tr(vIndexDensityDof,   w0*numDofPerCell,   M(0,0)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof,   w0*numDofPerCell+1, M(0,1)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof,   w0*numDofPerCell+2, M(0,2)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof+1, w0*numDofPerCell,   M(1,0)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof+1, w0*numDofPerCell+1, M(1,1)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof+1, w0*numDofPerCell+2, M(1,2)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof+2, w0*numDofPerCell,   M(2,0)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof+2, w0*numDofPerCell+1, M(2,1)) );
// 	  m_tripletList.push_back( Tr(vIndexDensityDof+2, w0*numDofPerCell+2, M(2,2)) );
// 	}

// 	if (j==2){
// 	  q(0,2) = (uMinusHalfNegFD(0)-uMinusHalfNeg(0))/eps;
// 	  q(1,2) = (uMinusHalfNegFD(1)-uMinusHalfNeg(1))/eps;
// 	  q(2,2) = (uMinusHalfNegFD(2)-uMinusHalfNeg(2))/eps;

// 	  q(0,4) = (uMinusHalfPosFD(0)-uMinusHalfPos(0))/eps;
// 	  q(1,4) = (uMinusHalfPosFD(1)-uMinusHalfPos(1))/eps;
// 	  q(2,4) = (uMinusHalfPosFD(2)-uMinusHalfPos(2))/eps;

// 	  q(0,7) = (uPlusHalfNegFD(0)-uPlusHalfNeg(0))/eps;
// 	  q(1,7) = (uPlusHalfNegFD(1)-uPlusHalfNeg(1))/eps;
// 	  q(2,7) = (uPlusHalfNegFD(2)-uPlusHalfNeg(2))/eps;

// 	  q(0,9) = (uPlusHalfPosFD(0)-uPlusHalfPos(0))/eps;
// 	  q(1,9) = (uPlusHalfPosFD(1)-uPlusHalfPos(1))/eps;
// 	  q(2,9) = (uPlusHalfPosFD(2)-uPlusHalfPos(2))/eps;
// 	}

// 	if (j==3){
// 	  q(0,5) = (uMinusHalfPosFD(0)-uMinusHalfPos(0))/eps;
// 	  q(1,5) = (uMinusHalfPosFD(1)-uMinusHalfPos(1))/eps;
// 	  q(2,5) = (uMinusHalfPosFD(2)-uMinusHalfPos(2))/eps;

// 	  q(0,8) = (uPlusHalfNegFD(0)-uPlusHalfNeg(0))/eps;
// 	  q(1,8) = (uPlusHalfNegFD(1)-uPlusHalfNeg(1))/eps;
// 	  q(2,8) = (uPlusHalfNegFD(2)-uPlusHalfNeg(2))/eps;

// 	  q(0,10) = (uPlusHalfPosFD(0)-uPlusHalfPos(0))/eps;
// 	  q(1,10) = (uPlusHalfPosFD(1)-uPlusHalfPos(1))/eps;
// 	  q(2,10) = (uPlusHalfPosFD(2)-uPlusHalfPos(2))/eps;
// 	}

// 	if (j==4){
// 	  q(0,11) = (uPlusHalfPosFD(0)-uPlusHalfPos(0))/eps;
// 	  q(1,11) = (uPlusHalfPosFD(1)-uPlusHalfPos(1))/eps;
// 	  q(2,11) = (uPlusHalfPosFD(2)-uPlusHalfPos(2))/eps;
// 	}
//       }
//   }
// };


}}}
#endif
