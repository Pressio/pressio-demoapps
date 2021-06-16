
#ifndef PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_
#define PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

#include "../stencil_filler.hpp"
#include "../reconstructor_from_stencil.hpp"
#include "./fluxes.hpp"

namespace pressiodemoapps{ namespace ad{ namespace impl{

template<
  class scalar_t,
  class mesh_t,
  class state_t,
  class velo_t
  >
class AdvectionAppT
{
  static_assert
  (std::is_same<scalar_t, typename mesh_t::scalar_t>::value, "");

public:
  using index_t		 = typename mesh_t::index_t;
  using scalar_type	 = scalar_t;
  using state_type	 = state_t;
  using velocity_type	 = state_t;
  using stencil_values_t = state_t;

  static constexpr int numDofPerCell{1};
  static constexpr int dimensionality{1};

public:
  explicit AdvectionAppT(const mesh_t & meshObj,
		      ::pressiodemoapps::Advection1d probEnum,
		      ::pressiodemoapps::ReconstructionType recEnum)
    : m_meshObj(meshObj), m_recEn(recEnum), m_probEn(probEnum)
  {
    m_numDofStencilMesh = m_meshObj.stencilMeshSize();
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize();

    const auto stencilSize = reconstructionTypeToStencilSize(recEnum);
    ::pressiodemoapps::resize(m_stencilVals, stencilSize);
  }

  index_t totalDofSampleMesh() const{ return m_numDofSampleMesh; }
  index_t totalDofStencilMesh() const{ return m_numDofStencilMesh; }

  state_type initialCondition() const
  {
    state_type res(m_numDofStencilMesh);

    const auto & x = m_meshObj.viewX();
    for (int i=0; i<::pressiodemoapps::extent(x,0); ++i){
      res(i) = std::sin(M_PI*x(i));
    }
    return res;
  }

  velocity_type createVelocity() const {
    velocity_type V(m_numDofSampleMesh);
    return V;
  }

  void velocity(const state_type & stateO,
		const scalar_type timeValue,
		velocity_type & veloO) const
  {
    velocityImpl(stateO, timeValue, veloO);
  }

private:
  void velocityImpl(const state_type & U,
		    const scalar_type t,
		    velocity_type & V) const
  {
    scalar_type FL{0};
    scalar_type FR{0};
    scalar_type uMinusHalfNeg{0};
    scalar_type uMinusHalfPos{0};
    scalar_type uPlusHalfNeg {0};
    scalar_type uPlusHalfPos {0};

    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_values_t, state_type, mesh_t, void>;
    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencil<
      scalar_type, stencil_values_t>;

    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    sfiller_t StencilFiller(stencilSize, U, m_meshObj, m_stencilVals);
    rec_fnct_t Reconstructor(m_recEn, m_stencilVals,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

    const auto dxInv = m_meshObj.dxInv();
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
    {
      StencilFiller(smPt);
      Reconstructor.template operator()<numDofPerCell>();
      linAdvRusanovFlux(FL, uMinusHalfNeg, uMinusHalfPos);
      linAdvRusanovFlux(FR, uPlusHalfNeg,  uPlusHalfPos);

      const auto vIndex = smPt*numDofPerCell;
      V(vIndex) = dxInv*(FL - FR);
    }
  }

private:
  const mesh_t & m_meshObj;
  ::pressiodemoapps::Advection1d m_probEn;
  ::pressiodemoapps::ReconstructionType m_recEn;

  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};
  mutable stencil_values_t m_stencilVals;
};

}}}//end namespace
#endif
