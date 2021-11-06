
#ifndef PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_
#define PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

#include "functor_fill_stencil.hpp"
#include "functor_reconstruct_from_state.hpp"

namespace pressiodemoapps{ namespace impladv{

template<typename sc_t>
void linAdvRusanovFlux(sc_t & F,
		       const sc_t & qL,
		       const sc_t & qR)
{
  //constexpr auto oneHalf = static_cast<sc_t>(0.5);
  // sc_t FL;
  // sc_t FR;
  // FL = qL;
  // FR = qR;
  F = qL; // = 0.5*( qL + qR + (qL - qR) );
}

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

  static constexpr int dimensionality{1};
  static constexpr int numDofPerCell{1};

public:
  explicit AdvectionAppT(const mesh_t & meshObj,
		      ::pressiodemoapps::Advection1d probEnum,
		      ::pressiodemoapps::InviscidFluxReconstruction recEnum)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum)
  {
    m_numDofStencilMesh = m_meshObj.stencilMeshSize();
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize();

    const auto stencilSize = reconstructionTypeToStencilSize(recEnum);
    ::pressiodemoapps::resize(m_stencilVals, stencilSize);
  }

  index_t totalDofSampleMesh()  const{ return m_numDofSampleMesh;  }
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

    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromState<
      dimensionality, scalar_type, state_type, mesh_t>;

    rec_fnct_t Reconstructor(m_recEn, U, m_meshObj,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

    const auto dxInv = m_meshObj.dxInv();
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
    {
      Reconstructor.template operator()<numDofPerCell>(smPt);
      linAdvRusanovFlux(FL, uMinusHalfNeg, uMinusHalfPos);
      linAdvRusanovFlux(FR, uPlusHalfNeg,  uPlusHalfPos);
      const auto vIndex = smPt*numDofPerCell;
      V(vIndex) = dxInv*(FL - FR);
    }
  }

private:
  const mesh_t & m_meshObj;
  ::pressiodemoapps::Advection1d m_probEn;
  ::pressiodemoapps::InviscidFluxReconstruction m_recEn;

  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};
  mutable stencil_values_t m_stencilVals;
};

}}//end namespace
#endif
