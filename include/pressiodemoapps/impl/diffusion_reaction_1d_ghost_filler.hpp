
#ifndef PRESSIODEMOAPPS_GHOST_FILLER_DIFFREAC_1D_HPP_
#define PRESSIODEMOAPPS_GHOST_FILLER_DIFFREAC_1D_HPP_

namespace pressiodemoapps{ namespace impldiffreac{

template<class state_t, class mesh_t, class ghost_t>
class GhostFillerProblemA1d
{

public:
  GhostFillerProblemA1d() = delete;
  GhostFillerProblemA1d(const int stencilSize,
			 const state_t & stateIn,
			 const mesh_t & meshIn,
			 ghost_t & ghostLeft,
			 ghost_t & ghostRight)
    : m_stencilSize(stencilSize),
      m_state(stateIn),
      m_meshObj(meshIn),
      m_ghostLeft(ghostLeft),
      m_ghostRight(ghostRight)
  {}

  void operator()()
  {
    if (m_stencilSize!=3){
      throw std::runtime_error("1d-diff-reaction: ghost filler not implemented for stencil !=3 ");
    }

    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    const auto & graph = m_meshObj.graph();
    const auto cellGIDL = graph(0, 0);
    const auto cellGIDR = graph(sampleMeshSize-1, 0);

    m_ghostLeft(0)  = -m_state(cellGIDL);
    m_ghostRight(0) = -m_state(cellGIDR);
  }

private:
  const int m_stencilSize;
  const state_t & m_state;
  const mesh_t & m_meshObj;
  ghost_t & m_ghostLeft;
  ghost_t & m_ghostRight;
};

}}
#endif