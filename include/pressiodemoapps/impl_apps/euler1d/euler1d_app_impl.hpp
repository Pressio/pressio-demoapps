
#ifndef PRESSIODEMOAPPS_EULER1D_APP_HPP_
#define PRESSIODEMOAPPS_EULER1D_APP_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

#include "../eulerCommon/energy.hpp"
#include "../eulerCommon/fluxes.hpp"
#include "./initial_condition.hpp"
#include "./ghost_filler.hpp"
#include "../stencil_filler.hpp"
#include "../reconstructor_from_stencil.hpp"
#include "../reconstructor_from_state.hpp"

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<
  class scalar_t,
  class mesh_t,
  class state_t,
  class velo_t,
  class ghost_t
  >
class Euler1dAppT
{
  static_assert
  (std::is_same<scalar_t, typename mesh_t::scalar_t>::value, "");

public:
  using index_t		 = typename mesh_t::index_t;
  using scalar_type	 = scalar_t;
  using state_type	 = state_t;
  using velocity_type	 = state_t;
  using ghost_type	 = ghost_t;
  using stencil_values_t = state_type;
  using flux_t		 = state_type;
  using edge_rec_t	 = state_type;

  static constexpr int dimensionality{1};
  static constexpr index_t numDofPerCell{3};

public:
  Euler1dAppT(const mesh_t & meshObj,
	      ::pressiodemoapps::Euler1d probEnum,
	      ::pressiodemoapps::ReconstructionType recEnum,
	      ::pressiodemoapps::FluxType fluxEnum)
    : m_meshObj(meshObj), m_probEn(probEnum), m_recEn(recEnum), m_fluxEn(fluxEnum)
  {
    m_numDofStencilMesh = m_meshObj.stencilMeshSize() * numDofPerCell;
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize() * numDofPerCell;

    const auto stencilSize = reconstructionTypeToStencilSize(recEnum);
    ::pressiodemoapps::resize(m_stencilVals, numDofPerCell*stencilSize);
    allocateGhosts();
  }

  state_type initialCondition() const
  {
    state_type res(m_numDofStencilMesh);

    if (m_probEn == ::pressiodemoapps::Euler1d::PeriodicSmooth){
      ::pressiodemoapps::ee::euler1dsineInitialCondition(res, m_meshObj, m_gamma);
    }
    else if (m_probEn == ::pressiodemoapps::Euler1d::Sod){
      ::pressiodemoapps::ee::sod1dInitialCondition(res, m_meshObj, m_gamma);
    }
    else if (m_probEn == ::pressiodemoapps::Euler1d::Lax){
      ::pressiodemoapps::ee::lax1dInitialCondition(res, m_meshObj, m_gamma);
    }
    else{
      //nothing
    }

    return res;
  }

  scalar_type gamma()		const{ return m_gamma; }
  index_t totalDofSampleMesh()  const{ return m_numDofSampleMesh; }
  index_t totalDofStencilMesh() const{ return m_numDofStencilMesh; }

  velocity_type createVelocity() const {
    velocity_type V(m_numDofSampleMesh);
    return V;
  }

  void velocity(const state_type & state,
		const scalar_type timeValue,
		velocity_type & V) const
  {
    velocityImpl(state, timeValue, V);
  }

private:
  void velocityImpl(const state_type & U,
		    const scalar_type t,
		    velocity_type & V) const
  {
    flux_t FL(numDofPerCell);
    flux_t FR(numDofPerCell);
    edge_rec_t uMinusHalfNeg(numDofPerCell);
    edge_rec_t uMinusHalfPos(numDofPerCell);
    edge_rec_t uPlusHalfNeg (numDofPerCell);
    edge_rec_t uPlusHalfPos (numDofPerCell);
    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);

    // only need ghosts for specific problems
    if (m_probEn == ::pressiodemoapps::Euler1d::Sod or
	m_probEn == ::pressiodemoapps::Euler1d::Lax)
    {
      using ghost_filler_t  = ::pressiodemoapps::impl::Ghost1dNeumannFiller<
	numDofPerCell, state_type, mesh_t, ghost_t>;
      ghost_filler_t ghF(stencilSize, U, m_meshObj, m_ghostLeft, m_ghostRight);
      ghF();
    }

    using sfiller_t  = ::pressiodemoapps::impl::StencilFiller<
      dimensionality, numDofPerCell, stencil_values_t, state_type, mesh_t, ghost_t>;
    using rec_fnct_t = ::pressiodemoapps::impl::ReconstructorFromStencil<
      edge_rec_t, stencil_values_t>;

    sfiller_t StencilFiller(stencilSize, U, m_meshObj,
			    m_ghostLeft, m_ghostRight, m_stencilVals);
    rec_fnct_t Reconstructor(m_recEn, m_stencilVals,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

    const auto dxInv = m_meshObj.dxInv();
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
      {
	StencilFiller(smPt);
	Reconstructor.template operator()<numDofPerCell>();
	switch(m_fluxEn)
	  {
	  case ::pressiodemoapps::FluxType::Rusanov:
	    eeRusanovFluxThreeDof(FL, uMinusHalfNeg, uMinusHalfPos, m_gamma);
	    eeRusanovFluxThreeDof(FR, uPlusHalfNeg,  uPlusHalfPos,  m_gamma);
	    break;
	  }

	const auto vIndex = smPt*numDofPerCell;
	V(vIndex)   = dxInv*(FL(0) - FR(0));
	V(vIndex+1) = dxInv*(FL(1) - FR(1));
	V(vIndex+2) = dxInv*(FL(2) - FR(2));
      }
  }

private:
  void allocateGhosts()
  {
    /*for stencil = 3, left bd indicated by ||
      ------------------
      |	 0,1,2        ||
      | rho, rho*u, E ||
      |		      ||
      -----------------

      for stencil = 7, at leftBoundary indicated by ||
      --------------------------------------
      |	 6,7,8      | 3,4,5   |  0,1,2    || U | e0 | e1 | e2
      | rho,        | rho,    | rho       ||
      | rho u,	    | rho*u   | rho*u     ||
      | E	    | E       | E         ||
      --------------------------------------
     */

    const auto stencilSize = reconstructionTypeToStencilSize(m_recEn);
    const auto numGhostValues = numDofPerCell*((stencilSize-1)/2);
    ::pressiodemoapps::resize(m_ghostLeft,  numGhostValues);
    ::pressiodemoapps::resize(m_ghostRight, numGhostValues);
  }

private:
  scalar_type m_gamma{1.4};

  const mesh_t & m_meshObj;
  ::pressiodemoapps::Euler1d m_probEn;
  ::pressiodemoapps::ReconstructionType m_recEn;
  ::pressiodemoapps::FluxType m_fluxEn;

  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points. for this problem,
  // the dof = numDofPerCell * number_of_unknown_grid_points
  // SampleMesh_ identifies the velocity/residual locations
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};

  mutable stencil_values_t m_stencilVals = {};
  mutable ghost_t m_ghostLeft;
  mutable ghost_t m_ghostRight;
};

}}}//end namespace
#endif
