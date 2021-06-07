
#ifndef PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_
#define PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_IMPL_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be
// improved later on but we have a starting point.

namespace pressiodemoapps{ namespace ad{ namespace impl{

template<
  class scalar_t,
  class mesh_t,
  class state_t,
  class velo_t
  >
class LinearAdvT
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
  explicit LinearAdvT(const mesh_t & meshObj,
		      ::pressiodemoapps::reconstructionEnum enIn =
		      ::pressiodemoapps::reconstructionEnum::firstOrder)
    : m_recEn(enIn), m_meshObj(meshObj)
  {
    m_numDofStencilMesh = m_meshObj.stencilMeshSize();
    m_numDofSampleMesh  = m_meshObj.sampleMeshSize();

    const auto stencilSize = reconstructionEnumToStencilSize(enIn);
    ::pressiodemoapps::resize(m_stencilVals, stencilSize);
  }

  index_t totalDofSampleMesh() const{ return m_numDofSampleMesh; }
  index_t totalDofStencilMesh() const{ return m_numDofStencilMesh; }

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
    using rec_fnct_t = ::pressiodemoapps::impl::Reconstructor<
      numDofPerCell, scalar_type, stencil_values_t>;

    const auto stencilSize = reconstructionEnumToStencilSize(m_recEn);
    sfiller_t StencilFiller(stencilSize, U, m_meshObj, m_stencilVals);
    rec_fnct_t Reconstructor(m_recEn, m_stencilVals,
			     uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg,  uPlusHalfPos);

    const auto dxInv = m_meshObj.dxInv();
    const auto sampleMeshSize = m_meshObj.sampleMeshSize();
    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
    {
      StencilFiller(smPt);
      Reconstructor();
      linAdvRusanovFlux(FL, uMinusHalfNeg, uMinusHalfPos);
      linAdvRusanovFlux(FR, uPlusHalfNeg,  uPlusHalfPos);

      const auto vIndex = smPt*numDofPerCell;
      V(vIndex) = dxInv*(FL - FR);
    }
  }

private:
  pressiodemoapps::reconstructionEnum m_recEn =
    pressiodemoapps::reconstructionEnum::firstOrder;

  const mesh_t & m_meshObj;
  index_t m_numDofStencilMesh = {};
  index_t m_numDofSampleMesh  = {};
  mutable stencil_values_t m_stencilVals;
};

}}}//end namespace
#endif
