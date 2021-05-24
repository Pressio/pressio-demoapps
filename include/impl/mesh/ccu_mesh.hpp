
#ifndef PRESSIODEMOAPPS_CCU_MESH_HPP_
#define PRESSIODEMOAPPS_CCU_MESH_HPP_

namespace pressiodemoapps{ namespace impl{

template<class scalar_type>
class CellCenteredUniformMeshEigen
{
public:
  using scalar_t = scalar_type;

  // type to use for all indexing, has to be large enough
  // to support indexing fairly large systems
  using index_t  = int32_t;

  using x_t = Eigen::Matrix<scalar_type,-1,1>;
  using y_t = x_t;
  using mesh_graph_t  = Eigen::Matrix<index_t,-1,-1,Eigen::RowMajor>;

  CellCenteredUniformMeshEigen(const std::string & meshDir)
  {
    pressiodemoapps::readMeshInfo(meshDir, dim_,
				  cellDeltas_,
				  stencilSize_,
				  numGptStencilMesh_,
				  numGptSampleMesh_);

    x_.resize(numGptStencilMesh_);
    y_.resize(numGptStencilMesh_);
    pressiodemoapps::readMeshCoordinates(meshDir, x_, y_);

    const auto graphSize = (stencilSize_-1)*dim_;
    graph_.resize(numGptSampleMesh_, graphSize+1);
    pressiodemoapps::readMeshConnectivity(meshDir, graph_);
  }

  const mesh_graph_t & graph() const{
    return graph_;
  }

  const scalar_type dx() const{
    return cellDeltas_[0];
  }

  const scalar_type dxInv() const{
    return cellDeltas_[1];
  }

  const scalar_type dy() const{
    return cellDeltas_[2];
  }

  const scalar_type dyInv() const{
    return cellDeltas_[3];
  }

  index_t stencilMeshSize() const{
    return numGptStencilMesh_;
  }

  index_t sampleMeshSize() const{
    return numGptSampleMesh_;
  }

  const int stencilSize() const {
    return stencilSize_;
  }

  const x_t & viewX() const{
    return x_;
  }

  const y_t & viewY() const{
    return y_;
  }

private:
  int dim_;
  std::array<scalar_type,4> cellDeltas_{};

  index_t numGptStencilMesh_ = {};
  index_t numGptSampleMesh_  = {};

  /*
    graph: contains a list such that
    1 0 3 2 -1

    first col   : contains GIDs of cells where we want velocity
    col 1,2,3,4 : contains GIDs of neighboring cells needed for stencil
    the order of the neighbors is: west, north, east, south
    if needed:
       col 4,5,6,7 : contains GIDs of degree2 neighboring cells needed for stencil
       the order of the neighbors is: west, north, east, south

       col 8,9,10,11 : contains GIDs of degree3 neighboring cells needed for stencil
       the order of the neighbors is: west, north, east, south
  */
  int stencilSize_ = {};
  mesh_graph_t graph_ = {};

  x_t x_ = {};
  y_t y_ = {};
};

}}
#endif
