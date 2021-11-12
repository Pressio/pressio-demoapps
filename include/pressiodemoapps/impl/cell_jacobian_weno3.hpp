
#ifndef PRESSIODEMOAPPS_CELL_JACOBIAN_WENO3_HPP_
#define PRESSIODEMOAPPS_CELL_JACOBIAN_WENO3_HPP_

namespace pressiodemoapps{ namespace ee{ namespace impl{


 /*
   Draft function to create a global Jacobian for a 1D, 3rd Order Weno discretization 

   Nomenclature diagram:

   cell indices:
   |   l1: u_{i-2}  |  l0: u_{i-1}  |    u_i    |   r0: u_{i+1}   |   r1: u_{i+2}   |
   reconstruction:              Lneg|Lpos   Rneg|Rpos


 */
// template<class ScalarType, class T>
// class EigenJacobianThreeDofPerCellFunctor
// {
//   SparseMatrix & J;
//   MeshType & meshObj;
//   CellJacobianType & m_JLneg;
//   CellJacobianType & m_JLpos;
//   CellJacobianType & m_JRneg;
//   CellJacobianType & m_JRpos;
//   RecontructionJacobianData m_recJacData;
//   // The recontruction jacobian data object requres the meshObj and cell_type as an input:
//   // cell types: interior, near boundary, on boundary
//   // 
//   // The object should have the following methods:
//   // setCellIndex(index_t smPt)
//   // getPertubations: get list of length 5 vectors with perturbations for states in all 5 stencil cells (these will depend on the boundary condition)
//   //     For the interior cells, the perturbation vectors are:
//   //     [eps, 0, 0, 0, 0], [ 0, eps, 0, 0, 0], [0, 0, eps, 0, 0], [0, 0, 0, eps, 0], [0, 0, 0, 0, eps]
//   //     For a neumann BC, the perturbation vectors are:
//   //     [eps, 0, 0, 0, eps], [0 , eps, 0 , eps, 0], [0, 0, eps, 0, 0]
//   //     For a left reflecting BC, the perturbation vectors are:
//   //     [-eps, 0, 0, 0, eps], [0 , -eps, 0 , eps, 0], [0, 0, eps, 0, 0]
//   //     For a left dirichlet BC:
//   //     [0, 0, eps, 0, 0], [0, 0, 0, eps, 0], [0, 0, 0, 0, eps] 
//   // getTrueIndex(1): get true cell index for l0. This will simply be cell i-1 when smPt corresponds to an inner cell. W0 will be cell i+1 when smPt is a cell on a Neumann boundary
//   // getTrueIndex(2): get true cell index for r0...
//   // getTrueIndex(3): get true cell index for l1...
//   // getTrueIndex(4): get true cell index for r1...
//
// public:
//   EigenJacobianThreeDofPerCellFunctorInnerCell() = delete;
//   EigenJacobianThreeDofPerCellFunctorInnerCell(SparseMatrix & J, MeshType & meshObj,
// 				      T & JLneg, T & JLpos, T & JRneg, T & JRpos, Enum cell_type):
//     m_J(J),
//     m_meshObj(meshObj),
//     m_JLneg(JLneg),
//     m_JLpos(JLpos),
//     m_JRneg(JRneg),
//     m_JRpos(JRpos),
//     m_recJacData(meshObj,cell_type)
//   {
//   }
//
//   template<class index_t>
//   void operator()(const index_t smPt)
//   {  
//     const auto & graph = m_meshObj.graph();
//
//     // Create matrix to store reconstruction Jacobian diagonals in
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
//
//      recJacData.setCellIndex(smPt)
//     const auto perturbations = recJacData.getPerturbations()
//     const auto l0 = recJacData.getTrueIndex(1) 
//     const auto r0 = recJacData.getTrueIndex(2) 
//     const auto l1 = recJacData.getTrueIndex(3) 
//     const auto r1 = recJacData.getTrueIndex(4) 
//     // For interior cells, this should be the same as:
//     // const auto l0 = graph(smPt, 1);
//     // const auto r0 = graph(smPt, 2);
//     // const auto l1 = graph(smPt, 3);
//     // const auto r1 = graph(smPt, 4);
//
//     index_t rowIndex  = smPt*numDofPerCell;
//     index_t col_l2    = l1*numDofPerCell;
//     index_t col_l1    = l0*numDofPerCell;
//     index_t colIndex  = graph(smPt, 0)*numDofPerCell;
//     index_t col_r0    = r0*numDofPerCell;
//     index_t col_r1    = r1*numDofPerCell;
//
//
//     // TODO warning, this is pseudo code...
//     const auto num_pertur = len(perturbations);
//     for (int j=0; j<num_perturb; ++j){
// 	   // Copy stencil cell states from unperturbed states
//       svals2 = m_stencilVals;
//
//       // TODO Perturb cell states in svals2
//       pertubation = perturbations[j]
//       for (int k=0; k<len(perturbation); k++){
//         for (int l=0; l<numDofPerCell; l++){
//         svals2(k,l) += perturbation[k] // perturb lth Dof of kth cell in stencil 
//         }
//       }
//
// 	   // assuming that we have svals2 that contains the perturb values
// 	   // we need to run the reconstruction
// 	   Reconstructor2.template operator()<numDofPerCell>();
// 	   // here, it means that, these are overwritten:
// 	   // uMinusHalfNegFD
// 	   // uMinusHalfPosFD
// 	   // uPlusHalfNegFD
// 	   // uPlusHalfPosFD
//
//       // Update global Jacobian based on cell peturbation vector
//       // Entries for Cell l1
//       if (perturbation[0]!=0){
// 	     q(0,0) = (uMinusHalfNegFD(0)-uMinusHalfNeg(0))/eps;
// 	     q(1,0) = (uMinusHalfNegFD(1)-uMinusHalfNeg(1))/eps;
// 	     q(2,0) = (uMinusHalfNegFD(2)-uMinusHalfNeg(2))/eps;
//         
//         auto M = JLneg * q.col(1).asDiagonalMatrix();
//
//         for (int k=0; k<numDofPerCell; ++k){
//           for (int l=0; l<numDofPerCell; ++l){
//	         m_J.coeffRef(rowIndex+k, col_l1+l) += M(k,l)*m_hInv;
//           }
//         }
// 	   }
//
//       // Entries for Cell l0
//       if (perturbation[1]!=0){
// 	     q(0,1) = (uMinusHalfNegFD(0)-uMinusHalfNeg(0))/eps;
// 	     q(1,1) = (uMinusHalfNegFD(1)-uMinusHalfNeg(1))/eps;
// 	     q(2,1) = (uMinusHalfNegFD(2)-uMinusHalfNeg(2))/eps;
// 	     auto M1 = JLneg * q.col(1).asDiagonalMatrix();
//
// 	     q(0,3) = (uMinusHalfPosFD(0)-uMinusHalfPos(0))/eps;
// 	     q(1,3) = (uMinusHalfPosFD(1)-uMinusHalfPos(1))/eps;
// 	     q(2,3) = (uMinusHalfPosFD(2)-uMinusHalfPos(2))/eps;
// 	     auto M2 = JLpos * q.col(3).asDiagonalMatrix();
//
// 	     q(0,6) = (uPlusHalfNegFD(0)-uPlusHalfNeg(0))/eps;
// 	     q(1,6) = (uPlusHalfNegFD(1)-uPlusHalfNeg(1))/eps;
// 	     q(2,6) = (uPlusHalfNegFD(2)-uPlusHalfNeg(2))/eps;
// 	     auto M3 = -JRneg * q.col(6).asDiagonalMatrix();
//
//         auto M = M1 + M2 + M3;
//         for (int k=0; k<numDofPerCell; ++k){
//           for (int l=0; l<numDofPerCell; ++l){
//	         m_J.coeffRef(rowIndex+k, col_l0+l) += M(k,l)*m_hInv;
//           }
//         }
//       }
//
//
//       // Entries for my Cell 
//       if (perturbation[2]!=0){
//         q(0,2) = (uMinusHalfNegFD(0)-uMinusHalfNeg(0))/eps;
//         q(1,2) = (uMinusHalfNegFD(1)-uMinusHalfNeg(1))/eps;
//         q(2,2) = (uMinusHalfNegFD(2)-uMinusHalfNeg(2))/eps;
//         auto M1 = JLneg * q.col(2).asDiagonalMatrix();
//
//         q(0,4) = (uMinusHalfPosFD(0)-uMinusHalfPos(0))/eps;
//         q(1,4) = (uMinusHalfPosFD(1)-uMinusHalfPos(1))/eps;
//         q(2,4) = (uMinusHalfPosFD(2)-uMinusHalfPos(2))/eps;
//         auto M2 = JLpos * q.col(4).asDiagonalMatrix();
//      
//         q(0,7) = (uPlusHalfNegFD(0)-uPlusHalfNeg(0))/eps;
//         q(1,7) = (uPlusHalfNegFD(1)-uPlusHalfNeg(1))/eps;
//         q(2,7) = (uPlusHalfNegFD(2)-uPlusHalfNeg(2))/eps;
//         auto M3 = -JRneg * q.col(7).asDiagonalMatrix();
//     
//         q(0,9) = (uPlusHalfPosFD(0)-uPlusHalfPos(0))/eps;
//         q(1,9) = (uPlusHalfPosFD(1)-uPlusHalfPos(1))/eps;
//         q(2,9) = (uPlusHalfPosFD(2)-uPlusHalfPos(2))/eps;
//         auto M4 = -JRpos * q.col(9).asDiagonalMatrix();
//
//         auto M = M1 + M2 + M3 + M4;
//         for (int k=0; k<numDofPerCell; ++k){
//           for (int l=0; l<numDofPerCell; ++l){
//	         m_J.coeffRef(rowIndex+k, colIndex+l) += M(k,l)*m_hInv;
//           }
//         }
//       }
//
//       // Entries for Cell r0
//       if (perturbation[3]!=0){
//         q(0,5) = (uMinusHalfPosFD(0)-uMinusHalfPos(0))/eps;
//         q(1,5) = (uMinusHalfPosFD(1)-uMinusHalfPos(1))/eps;
//         q(2,5) = (uMinusHalfPosFD(2)-uMinusHalfPos(2))/eps;
//         auto M2 = JLpos * q.col(5).asDiagonalMatrix();
//      
//         q(0,8) = (uPlusHalfNegFD(0)-uPlusHalfNeg(0))/eps;
//         q(1,8) = (uPlusHalfNegFD(1)-uPlusHalfNeg(1))/eps;
//         q(2,8) = (uPlusHalfNegFD(2)-uPlusHalfNeg(2))/eps;
//         auto M3 = -JRneg * q.col(8).asDiagonalMatrix();
//      
//         q(0,10) = (uPlusHalfPosFD(0)-uPlusHalfPos(0))/eps;
//         q(1,10) = (uPlusHalfPosFD(1)-uPlusHalfPos(1))/eps;
//         q(2,10) = (uPlusHalfPosFD(2)-uPlusHalfPos(2))/eps;
//         auto M4 = -JRpos * q.col(10).asDiagonalMatrix();
//
//         auto M = M2 + M3 + M4;
//         for (int k=0; k<numDofPerCell; ++k){
//           for (int l=0; l<numDofPerCell; ++l){
//	         m_J.coeffRef(rowIndex+k, col_r0+l) += M(k,l)*m_hInv;
//           }
//         }
//
//       }
//
//       // Entries for Cell r1
//       if (j==4){
//         q(0,11) = (uPlusHalfPosFD(0)-uPlusHalfPos(0))/eps;
//         q(1,11) = (uPlusHalfPosFD(1)-uPlusHalfPos(1))/eps;
//         q(2,11) = (uPlusHalfPosFD(2)-uPlusHalfPos(2))/eps;
//         auto M = -JRpos * q.col(10).asDiagonalMatrix();
//         for (int k=0; k<numDofPerCell; ++k){
//           for (int l=0; l<numDofPerCell; ++l){
//	         m_J.coeffRef(rowIndex+k, col_r1+l) += M(k,l)*m_hInv;
//           }
//         }
//
//
//       }
//     }
//   }
// };



}}}
#endif
