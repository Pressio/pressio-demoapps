
#ifndef PRESSIODEMOAPPS_EDGE_CELL_JACOBIAN_FIRST_ORDER_ADV_HPP_
#define PRESSIODEMOAPPS_EDGE_CELL_JACOBIAN_FIRST_ORDER_ADV_HPP_

namespace pressiodemoapps{ namespace impladv{

// template<int dim, class IndexT, class GraphType>
// void _get_left_right_cells_indices(IndexT & l0, IndexT & r0,
// 				   const IndexT & smPt,
// 				   const GraphType & graph,
// 				   int axis)
// {
//   if (dim == 1){
//     l0 = graph(smPt, 1);
//     r0 = graph(smPt, 2);
//   }
//   else if(dim==2){
//     l0 = (axis == 1) ? graph(smPt, 1) : graph(smPt, 4);
//     r0 = (axis == 1) ? graph(smPt, 3) : graph(smPt, 2);
//   }
//   else if(dim==3){
//     l0 = (axis == 1) ? graph(smPt, 1) : (axis==2) ? graph(smPt, 4) : graph(smPt, 5);
//     r0 = (axis == 1) ? graph(smPt, 3) : (axis==2) ? graph(smPt, 2) : graph(smPt, 6);
//   }
// }

}}
#endif
