
#ifndef PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_HPP_
#define PRESSIODEMOAPPS_LINEAR_ADVECTION_1D_HPP_

// note that the code below is pretty ugly for now, but it works.
// this is inside impl namespace for a reason, and will need to be 
// improved later on but we have a starting point.

namespace pressiodemoapps{ namespace ad{ namespace impl{

template<class scalar_t, class mesh_t>
class EigenLinearAdv
{
  static_assert
  (std::is_same<scalar_t, typename mesh_t::scalar_t>::value, "");

public:
  using scalar_type	= scalar_t;
  using state_type	= Eigen::Matrix<scalar_t,-1,1>;
  using velocity_type	= Eigen::Matrix<scalar_t,-1,1>;
  using dense_matrix_type = Eigen::Matrix<scalar_t,-1,-1>;
  using ghost_t = Eigen::Matrix<scalar_type,-1,1>;
  using index_t  = typename mesh_t::index_t;

public:
  EigenLinearAdv(const mesh_t & meshObj)
    : meshObj_(meshObj)
  {
    // here we have one dof per cell
    numDofStencilMesh_ = meshObj_.stencilMeshSize();
    numDofSampleMesh_  = meshObj_.sampleMeshSize();

    stencilVals_.resize(meshObj_.stencilSize());
  }

  index_t totalDofSampleMesh() const{
    return numDofSampleMesh_;
  }

  index_t totalDofStencilMesh() const{
    return numDofStencilMesh_;
  }

  velocity_type createVelocity() const {
    velocity_type V(numDofSampleMesh_);
    return V;
  }


  template<class adlist_t, class uindex_t>
  void fillStencil(index_t & smPt,
		   const uindex_t uIndex,
		   const adlist_t & adList,
		   const state_type & U) const
  {
    const auto w0 = adList[1];
    const auto e0 = adList[2];
    const auto w1 = adList[3];
    const auto e1 = adList[4];
    const auto w2 = adList[5];
    const auto e2 = adList[6];
    const auto stencilSize = meshObj_.stencilSize();

    if (stencilSize ==3)
    {
      stencilVals_[0] = U[w0];
      stencilVals_[1] = U[uIndex];
      stencilVals_[2] = U[e0];
    }
    else if (stencilSize ==7){
      stencilVals_[0] = U[w2];
      stencilVals_[1] = U[w1];
      stencilVals_[2] = U[w0];
      stencilVals_[3] = U[uIndex];
      stencilVals_[4] = U[e0];
      stencilVals_[5] = U[e1];
      stencilVals_[6] = U[e2];
    }
  }

  template<class T, class T2>
  void reconstruct(const T2 & stencilVals_,
		   T & uMinusHalfNeg,
		   T & uMinusHalfPos,
		   T & uPlusHalfNeg,
		   T & uPlusHalfPos) const
  {

    const auto stencilSize = meshObj_.stencilSize();
    if (stencilSize==3)
    {
      uMinusHalfNeg = stencilVals_[0];
      uMinusHalfPos = stencilVals_[1];
      uPlusHalfNeg  = stencilVals_[1];
      uPlusHalfPos  = stencilVals_[2];
    }
    else if (stencilSize==7)
    {
      pressiodemoapps::weno5(uMinusHalfNeg, uMinusHalfPos,
			     uPlusHalfNeg, uPlusHalfPos,
			     stencilVals_);
    }
  }

  void velocity(const state_type & U,
		const scalar_type t,
		velocity_type & V) const
  {
    scalar_type FL{0};
    scalar_type FR{0};
    scalar_type uMinusHalfNeg{0};
    scalar_type uMinusHalfPos{0};
    scalar_type uPlusHalfNeg {0};
    scalar_type uPlusHalfPos {0};

    const auto x = meshObj_.viewX();
    const auto dxInv = meshObj_.dxInv();
    const auto & graph = meshObj_.graph();
    const auto sampleMeshSize = meshObj_.sampleMeshSize();

    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
    {
      const auto thisCellAdList = graph.row(smPt);
      const auto cellGID = thisCellAdList(0);

      fillStencil(smPt, cellGID, thisCellAdList, U);

      // do edge reconstruction
      reconstruct(stencilVals_,
		  uMinusHalfNeg, uMinusHalfPos,
		  uPlusHalfNeg, uPlusHalfPos);

      // compute fluxes
      linAdvRusanovFlux(FL, uMinusHalfNeg, uMinusHalfPos);
      linAdvRusanovFlux(FR, uPlusHalfNeg,  uPlusHalfPos);

      V(cellGID) = dxInv*(FL - FR);
    }
  }

private:
  const mesh_t & meshObj_;

  index_t numDofStencilMesh_ = {};
  index_t numDofSampleMesh_  = {};

  mutable std::vector<scalar_type> stencilVals_;
};

}}}//end namespace

#endif



// void fillGhosts(const state_type & U,
// 	      scalar_type t) const
// {
//   const auto dx  = meshObj_.dx();
//   const auto & x = meshObj_.viewX();
//   const auto & graph = meshObj_.graph();
//   const auto sampleMeshSize = meshObj_.sampleMeshSize();
//   const auto stencilSize = meshObj_.stencilSize();

//   if (graph.row(0)[1] != -1){
//     throw std::runtime_error("Adv1d: cell not near BD");
//   }
//   if (graph.row(sampleMeshSize-1)[2] != -1){
//     throw std::runtime_error("Adv1d: cell not near BD");
//   }

//   // assume here that smPt=0 is at left bd and
//   // smPt = sampleMeshSize-1 is at right bd
//   const auto cellGIDL = graph.row(0)[0];
//   const auto Xleft = x(cellGIDL);
//   const auto cellGIDR = graph.row(sampleMeshSize-1)[0];
//   const auto Xright = x(cellGIDR);

//   if (stencilSize == 3)
//     {
//   	ghostLeft_(0)  = analytical(Xleft -dx, t);
//   	ghostRight_(0) = analytical(Xright+dx, t);
//     }
//   else if (stencilSize == 7)
//     {
// 	ghostLeft_(2) = analytical(Xleft-3.*dx, t);
// 	ghostLeft_(1) = analytical(Xleft-2.*dx, t);
// 	ghostLeft_(0) = analytical(Xleft-dx, t);

// 	ghostRight_(0) = analytical(Xright+dx, t);
// 	ghostRight_(1) = analytical(Xright+2.*dx, t);
// 	ghostRight_(2) = analytical(Xright+3.*dx, t);
//     }
// }

// private:
//   void allocateGhosts()
//   {
//     /*
//       for stencil = 3
//       --------  --------
//       |	0   ||  ||  0  |
//       | u   ||  ||  U  |
//       --------  -------

//       for stencil = 7
//       --------------   --------------
//       |	2 | 1 | 0 ||   || 0 | 1 | 2 |
//       | u | u | u ||   || u | u | u |
//       --------------   --------------
//     */

//     const auto stencilSize = meshObj_.stencilSize();
//     const auto numGhostValues = ((stencilSize-1)/2);
//     ghostLeft_.resize(numGhostValues);
//     ghostRight_.resize(numGhostValues);
//   }
