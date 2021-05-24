
#ifndef PRESSIODEMOAPPS_EULER1D_APP_HPP_
#define PRESSIODEMOAPPS_EULER1D_APP_HPP_

/*
  probid = 1 (Sod) taken from: https://www.mdpi.com/2227-7390/6/10/211/pdf
 */

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<class scalar_t, class mesh_t, int probid>
class EigenApp1d
{
  static_assert
  (std::is_same<scalar_t, typename mesh_t::scalar_t>::value, "");

  static_assert
  (probid==0 or probid==1, "Invalid probid");

public:
  using scalar_type	= scalar_t;
  using state_type	= Eigen::Matrix<scalar_t,-1,1>;
  using velocity_type	= Eigen::Matrix<scalar_t,-1,1>;
  using dense_matrix_type = Eigen::Matrix<scalar_t,-1,-1>;
  using ghost_t = Eigen::Matrix<scalar_t,-1,1>;
  using index_t  = typename mesh_t::index_t;

  static constexpr index_t numDofPerCell_{3};

public:
  EigenApp1d(const mesh_t & meshObj)
    : meshObj_(meshObj)
  {
    allocateGhosts();
  }

  // initial condition for Sod
  template<int _probid = probid>
  typename std::enable_if<_probid==1, state_type>::type
  initialCondition() const
  {
    state_type res(numDofStencilMesh_);
    const auto x= meshObj_.viewX();
    std::array<scalar_type, 3> prim;
    for (int i=0; i<x.size(); ++i)
      {
	if (x(i) <= 0.0){
	  prim[0] = 1.0;
	  prim[1] = 0.0;
	  prim[2] = 1.0;
	}

	if (x(i) > 0.0){
	  prim[0] = 0.125;
	  prim[1] = 0.0;
	  prim[2] = 0.1;
	}

	const auto ind = i*3;
	res(ind)   = prim[0];
	res(ind+1) = prim[0]*prim[1];
	res(ind+2) = pressiodemoapps::ee::computeEnergyFromPrimitive(gamma_, prim);
      }
    return res;
  }

  scalar_type gamma() const{
    return gamma_;
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
    const auto stencilSize = meshObj_.stencilSize();

    // indices of the neighboring cells
    const auto w0 = adList[1];
    const auto e0 = adList[2];
    const auto w1 = adList[3];
    const auto e1 = adList[4];
    const auto w2 = adList[5];
    const auto e2 = adList[6];

    // indices of the states from neighboring cells
    const auto w0i = w0*numDofPerCell_;
    const auto e0i = e0*numDofPerCell_;
    const auto w1i = w1*numDofPerCell_;
    const auto e1i = e1*numDofPerCell_;
    const auto w2i = w2*numDofPerCell_;
    const auto e2i = e2*numDofPerCell_;

    if (stencilSize ==3)
    {
      if (w0==-1){
	stencilVals_[0] = ghostLeft_[0];
	stencilVals_[1] = ghostLeft_[1];
	stencilVals_[2] = ghostLeft_[2];
      }else{
	stencilVals_[0] = U[w0i];
	stencilVals_[1] = U[w0i+1];
	stencilVals_[2] = U[w0i+2];
      }

      stencilVals_[3] = U[uIndex];
      stencilVals_[4] = U[uIndex+1];
      stencilVals_[5] = U[uIndex+2];

      if (e0==-1){
	stencilVals_[6] = ghostRight_[0];
	stencilVals_[7] = ghostRight_[1];
	stencilVals_[8] = ghostRight_[2];
      }else{
	stencilVals_[6] = U[e0i];
	stencilVals_[7] = U[e0i+1];
	stencilVals_[8] = U[e0i+2];
      }
    }
    else if (stencilSize ==7)
    {

      if (w2==-1){
	stencilVals_[0] = ghostLeft_[6];
	stencilVals_[1] = ghostLeft_[7];
	stencilVals_[2] = ghostLeft_[8];
      }else{
	stencilVals_[0] = U[w2i];
	stencilVals_[1] = U[w2i+1];
	stencilVals_[2] = U[w2i+2];
      }

      if (w1==-1){
	stencilVals_[3] = ghostLeft_[3];
	stencilVals_[4] = ghostLeft_[4];
	stencilVals_[5] = ghostLeft_[5];
      }else{
	stencilVals_[3] = U[w1i];
	stencilVals_[4] = U[w1i+1];
	stencilVals_[5] = U[w1i+2];
      }

      if (w0==-1){
	stencilVals_[6] = ghostLeft_[0];
	stencilVals_[7] = ghostLeft_[1];
	stencilVals_[8] = ghostLeft_[2];
      }else{
	stencilVals_[6] = U[w0i];
	stencilVals_[7] = U[w0i+1];
	stencilVals_[8] = U[w0i+2];
      }

      stencilVals_[9]  = U[uIndex];
      stencilVals_[10] = U[uIndex+1];
      stencilVals_[11] = U[uIndex+2];

      if (e0==-1){
	stencilVals_[12] = ghostRight_[0];
	stencilVals_[13] = ghostRight_[1];
	stencilVals_[14] = ghostRight_[2];
      }else{
	stencilVals_[12] = U[e0i];
	stencilVals_[13] = U[e0i+1];
	stencilVals_[14] = U[e0i+2];
      }

      if (e1==-1){
	stencilVals_[15] = ghostRight_[3];
	stencilVals_[16] = ghostRight_[4];
	stencilVals_[17] = ghostRight_[5];
      }else{
	stencilVals_[15] = U[e1i];
	stencilVals_[16] = U[e1i+1];
	stencilVals_[17] = U[e1i+2];
      }

      if (e2==-1){
	stencilVals_[18] = ghostRight_[6];
	stencilVals_[19] = ghostRight_[7];
	stencilVals_[20] = ghostRight_[8];
      }else{
	stencilVals_[18] = U[e2i];
	stencilVals_[19] = U[e2i+1];
	stencilVals_[20] = U[e2i+2];
      }
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

    if (stencilSize==3){
      uMinusHalfNeg = {stencilVals_[0], stencilVals_[1], stencilVals_[2]};
      uMinusHalfPos = {stencilVals_[3], stencilVals_[4], stencilVals_[5]};
      uPlusHalfNeg  = {stencilVals_[3], stencilVals_[4], stencilVals_[5]};
      uPlusHalfPos  = {stencilVals_[6], stencilVals_[7], stencilVals_[8]};
    }
    else if (stencilSize==7)
    {

      for (int dof=0; dof<numDofPerCell_; ++dof)
      {
	std::array<scalar_type,7> mys;
	int start = dof;
	for (int k=0; k<7; ++k){
	  mys[k] = stencilVals_[start + 3*k];
	}
	pressiodemoapps::weno5(uMinusHalfNeg[dof], uMinusHalfPos[dof],
			       uPlusHalfNeg[dof],  uPlusHalfPos[dof], mys);
      }

    }
  }

  void fillGhosts(const state_type & U,
		  scalar_type t) const
  {
    std::array<scalar_type, 4> prim = {0,0,0};

    const auto dx  = meshObj_.dx();
    const auto x  = meshObj_.viewX();
    const auto & graph = meshObj_.graph();
    const auto sampleMeshSize = meshObj_.sampleMeshSize();
    const auto stencilSize = meshObj_.stencilSize();

    // assume here that smPt=0 is at left bd and
    // smPt = sampleMeshSize-1 is at right bd
    const auto cellGIDL = graph.row(0)[0];
    const auto myXL = x(cellGIDL);
    if (graph.row(0)[1] != -1){
      throw std::runtime_error("Point not on boundary, something wrong");
    }

    const auto cellGIDR = graph.row(sampleMeshSize-1)[0];
    const auto myXR = x(cellGIDR);
    if (graph.row(sampleMeshSize-1)[2] != -1){
      throw std::runtime_error("Point not on boundary, something wrong");
    }

    // will need to fix this for the sample mesh
    // bcause it is not for sure we have cells right next to bd
    // that we can use to fill ghosts

    ghostLeft_(0) = U[0];
    ghostLeft_(1) = U[1];
    ghostLeft_(2) = U[2];

    ghostRight_(0) = U[cellGIDR*numDofPerCell_];
    ghostRight_(1) = U[cellGIDR*numDofPerCell_+1];
    ghostRight_(2) = U[cellGIDR*numDofPerCell_+2];

    if (stencilSize==7)
      {
	ghostLeft_(6) = U[2*numDofPerCell_];
	ghostLeft_(7) = U[2*numDofPerCell_+1];
	ghostLeft_(8) = U[2*numDofPerCell_+2];

	ghostLeft_(3) = U[1*numDofPerCell_];
	ghostLeft_(4) = U[1*numDofPerCell_+1];
	ghostLeft_(5) = U[1*numDofPerCell_+2];

	auto cid = cellGIDR-1;
	ghostRight_(3) = U[cid*numDofPerCell_];
	ghostRight_(4) = U[cid*numDofPerCell_+1];
	ghostRight_(5) = U[cid*numDofPerCell_+2];

	cid = cellGIDR-2;
	ghostRight_(6) = U[cid*numDofPerCell_];
	ghostRight_(7) = U[cid*numDofPerCell_+1];
	ghostRight_(8) = U[cid*numDofPerCell_+2];
      }
  }

  void velocity(const state_type & U,
		const scalar_type t,
		velocity_type & V) const
  {
    using arr_t = std::array<scalar_type, numDofPerCell_>;
    arr_t FL = {0,0,0};
    arr_t FR = {0,0,0};
    arr_t uMinusHalfNeg = {0,0,0};
    arr_t uMinusHalfPos = {0,0,0};
    arr_t uPlusHalfNeg  = {0,0,0};
    arr_t uPlusHalfPos  = {0,0,0};

    const auto dx  = meshObj_.dx();
    const auto dxInv  = meshObj_.dxInv();
    const auto & graph = meshObj_.graph();
    const auto sampleMeshSize = meshObj_.sampleMeshSize();
    const auto stencilSize = meshObj_.stencilSize();

    // only need ghosts when doing non periodic problem
    if (probid!=0){
      fillGhosts(U, t);
    }

    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
      {
	// ajacency list of this cell
	const auto thisCellAdList = graph.row(smPt);
	// gID of this cell
	const auto & cellGID = thisCellAdList(0);
	// given a cell, compute index in V of the correspond first dof
	const auto vIndex = smPt*numDofPerCell_;
	// given a cell, compute index in u of the correspond first dof
	const auto uIndex = cellGID*numDofPerCell_;

	fillStencil(smPt, uIndex, thisCellAdList, U);

	reconstruct(stencilVals_,
		    uMinusHalfNeg, uMinusHalfPos,
		    uPlusHalfNeg, uPlusHalfPos);


	eeRusanovFlux(FL, uMinusHalfNeg, uMinusHalfPos, gamma_);
	eeRusanovFlux(FR, uPlusHalfNeg,  uPlusHalfPos,  gamma_);

	V(vIndex)   = dxInv*(FL[0] - FR[0]);
	V(vIndex+1) = dxInv*(FL[1] - FR[1]);
	V(vIndex+2) = dxInv*(FL[2] - FR[2]);
      }
  }

private:
  void allocateGhosts()
  {
    numDofStencilMesh_ = meshObj_.stencilMeshSize() * numDofPerCell_;
    numDofSampleMesh_  = meshObj_.sampleMeshSize() * numDofPerCell_;

    const auto stencilSize = meshObj_.stencilSize();
    stencilVals_.resize(numDofPerCell_ * stencilSize);

    /*
      for stencil = 3, left bd indicated by ||:
      ------------------
      |	 0,1,2        ||
      | rho, rho*u, E ||
      |		      ||
      -----------------

      for stencil = 7, at leftBoundary indicated by ||:
      --------------------------------------
      |	 6,7,8      | 3,4,5   |  0,1,2    ||
      | rho,        | rho,    | rho       ||
      | rho u,	    | rho*u   | rho*u     ||
      | E	    | E       | E         ||
      --------------------------------------
     */

    const auto numGhostValues = numDofPerCell_*((stencilSize-1)/2);
    ghostLeft_.resize(numGhostValues);
    ghostRight_.resize(numGhostValues);
  }

private:
  const scalar_type gamma_ = (5.+2.)/5.;
  const scalar_type gammaMinusOne_ = gamma_ - 1.;
  const scalar_type gammaMinusOneInv_ = 1./gammaMinusOne_;
  const scalar_type gammaMinusOneDiv16_ = gammaMinusOne_/(8. * gamma_ * M_PI * M_PI);

  mutable index_t uWestIndex_  = {};
  mutable index_t uEastIndex_  = {};

  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points. for this problem,
  // the dof = numDofPerCell_ * number_of_unknown_grid_points
  // SampleMesh_ identifies the velocity/residual locations
  index_t numDofStencilMesh_ = {};
  index_t numDofSampleMesh_  = {};

  const std::array<scalar_type, 2> normalX_{1, 0};

  mutable ghost_t ghostLeft_;
  mutable ghost_t ghostRight_;

  mutable std::vector<scalar_type> stencilVals_;
  const mesh_t & meshObj_;
};

}}}//end namespace

#endif
