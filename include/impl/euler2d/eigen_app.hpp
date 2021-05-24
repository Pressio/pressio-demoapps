
#ifndef PRESSIODEMOAPPS_EULER2D_HPP_
#define PRESSIODEMOAPPS_EULER2D_HPP_

/*
  problemId = 0: generic problem with periodic BC so ghosts are not used
  problemId = 1: Sedov
  problemId = 2: Riemann
 */

namespace pressiodemoapps{ namespace ee{ namespace impl{

template<class scalar_t, class mesh_t, int problemId>
class EigenApp2d
{
  static_assert
  (std::is_same<scalar_t, typename mesh_t::scalar_t>::value, "");

  static_assert
  (problemId==0 or problemId==1 or problemId==2,
   "Invalid problemId");

public:
  using scalar_type	= scalar_t;
  using state_type	= Eigen::Matrix<scalar_t,-1,1>;
  using velocity_type	= Eigen::Matrix<scalar_t,-1,1>;
  using dense_matrix_type = Eigen::Matrix<scalar_t,-1,-1>;
  using ghost_t = Eigen::Matrix<scalar_t,-1,-1,Eigen::RowMajor>;
  using index_t  = typename mesh_t::index_t;

  static constexpr index_t numDofPerCell_{4};

public:
  explicit EigenApp2d(const mesh_t & meshObj,
		      int icId = 1)
    : meshObj_(meshObj), icId_(icId)
  {
    numDofStencilMesh_ = meshObj_.stencilMeshSize() * numDofPerCell_;
    numDofSampleMesh_  = meshObj_.sampleMeshSize() * numDofPerCell_;

    const auto stencilSize = meshObj_.stencilSize();
    stencilValsX_.resize(numDofPerCell_ * stencilSize);
    stencilValsY_.resize(numDofPerCell_ * stencilSize);

    allocateGhosts();
  }

  state_type initialCondition() const
  {
    state_type IC(numDofStencilMesh_);

    if (problemId == 1){
      sedov2dInitialCondition(IC, meshObj_, numDofPerCell_, gamma_);
    }
    else if (problemId == 2 and icId_ == 1){
      riemann2dInitialCondition1(IC, meshObj_, numDofPerCell_, gamma_);
    }
    else if (problemId == 2 and icId_ == 2){
      riemann2dInitialCondition2(IC, meshObj_, numDofPerCell_, gamma_);
    }

    return IC;
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

  template<class adlist_t, class uindex_t, class gh_t, class st_t>
  void fillStencil3(index_t & smPt,
		    const uindex_t uIndex,
		    const adlist_t & adList,
		    const state_type & U,
		    const gh_t & ghL,
		    const gh_t & ghR,
		    st_t & stencilVals,
		    int dim) const
  {
    const auto l0 = (dim == 1) ? adList[1] : adList[4];
    const auto r0 = (dim == 1) ? adList[3] : adList[2];

    if (l0 == -1){
      stencilVals[0] = ghL(smPt, 0);
      stencilVals[1] = ghL(smPt, 1);
      stencilVals[2] = ghL(smPt, 2);
      stencilVals[3] = ghL(smPt, 3);
    }else{
      const auto index = l0*numDofPerCell_;
      stencilVals[0]  = U[index];
      stencilVals[1]  = U[index+1];
      stencilVals[2] = U[index+2];
      stencilVals[3] = U[index+3];
    }

    stencilVals[4] = U[uIndex];
    stencilVals[5] = U[uIndex+1];
    stencilVals[6] = U[uIndex+2];
    stencilVals[7] = U[uIndex+3];

    if (r0 == -1){
      stencilVals[8]  = ghR(smPt, 0);
      stencilVals[9]  = ghR(smPt, 1);
      stencilVals[10] = ghR(smPt, 2);
      stencilVals[11] = ghR(smPt, 3);
    }else{
      const auto index = r0*numDofPerCell_;
      stencilVals[8]  = U[index];
      stencilVals[9]  = U[index+1];
      stencilVals[10] = U[index+2];
      stencilVals[11] = U[index+3];
    }
  }

  template<class adlist_t, class uindex_t, class gh_t, class st_t>
  void fillStencil7(index_t & smPt,
		    const uindex_t uIndex,
		    const adlist_t & adList,
		    const state_type & U,
		    const gh_t & ghL,
		    const gh_t & ghR,
		    st_t & stencilVals,
		    int dim) const
  {
    const auto l0 = (dim == 1) ? adList[1] : adList[4];
    const auto r0 = (dim == 1) ? adList[3] : adList[2];
    const auto l1 = (dim == 1) ? adList[5] : adList[8];
    const auto r1 = (dim == 1) ? adList[7] : adList[6];
    const auto l2 = (dim == 1) ? adList[9] : adList[12];
    const auto r2 = (dim == 1) ? adList[11] : adList[10];

    if (l2 == -1){
      stencilVals[0] = ghL(smPt, 8);
      stencilVals[1] = ghL(smPt, 9);
      stencilVals[2] = ghL(smPt, 10);
      stencilVals[3] = ghL(smPt, 11);
    }else{
      const auto index = l2*numDofPerCell_;
      stencilVals[0] = U[index];
      stencilVals[1] = U[index+1];
      stencilVals[2] = U[index+2];
      stencilVals[3] = U[index+3];
    }

    if (l1 == -1){
      stencilVals[4] = ghL(smPt, 4);
      stencilVals[5] = ghL(smPt, 5);
      stencilVals[6] = ghL(smPt, 6);
      stencilVals[7] = ghL(smPt, 7);
    }else{
      const auto index = l1*numDofPerCell_;
      stencilVals[4] = U[index];
      stencilVals[5] = U[index+1];
      stencilVals[6] = U[index+2];
      stencilVals[7] = U[index+3];
    }

    if (l0 == -1){
      stencilVals[8]  = ghL(smPt, 0);
      stencilVals[9]  = ghL(smPt, 1);
      stencilVals[10] = ghL(smPt, 2);
      stencilVals[11] = ghL(smPt, 3);
    }else{
      const auto index = l0*numDofPerCell_;
      stencilVals[8]  = U[index];
      stencilVals[9]  = U[index+1];
      stencilVals[10] = U[index+2];
      stencilVals[11] = U[index+3];
    }

    stencilVals[12] = U[uIndex];
    stencilVals[13] = U[uIndex+1];
    stencilVals[14] = U[uIndex+2];
    stencilVals[15] = U[uIndex+3];

    if (r0 == -1){
      stencilVals[16] = ghR(smPt, 0);
      stencilVals[17] = ghR(smPt, 1);
      stencilVals[18] = ghR(smPt, 2);
      stencilVals[19] = ghR(smPt, 3);
    }else{
      const auto index = r0*numDofPerCell_;
      stencilVals[16] = U[index];
      stencilVals[17] = U[index+1];
      stencilVals[18] = U[index+2];
      stencilVals[19] = U[index+3];
    }

    if (r1 == -1){
      stencilVals[20] = ghR(smPt, 4);
      stencilVals[21] = ghR(smPt, 5);
      stencilVals[22] = ghR(smPt, 6);
      stencilVals[23] = ghR(smPt, 7);
    }else{
      const auto index = r1*numDofPerCell_;
      stencilVals[20] = U[index];
      stencilVals[21] = U[index+1];
      stencilVals[22] = U[index+2];
      stencilVals[23] = U[index+3];
    }

    if (r2 == -1){
      stencilVals[24] = ghR(smPt, 8);
      stencilVals[25] = ghR(smPt, 9);
      stencilVals[26] = ghR(smPt, 10);
      stencilVals[27] = ghR(smPt, 11);
    }else{
      const auto index = r2*numDofPerCell_;
      stencilVals[24] = U[index];
      stencilVals[25] = U[index+1];
      stencilVals[26] = U[index+2];
      stencilVals[27] = U[index+3];
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
      uMinusHalfNeg = {stencilVals_[0], stencilVals_[1], stencilVals_[2], stencilVals_[3]};
      uMinusHalfPos = {stencilVals_[4], stencilVals_[5], stencilVals_[6], stencilVals_[7]};
      uPlusHalfNeg = uMinusHalfPos;
      uPlusHalfPos  = {stencilVals_[8], stencilVals_[9], stencilVals_[10], stencilVals_[11]};
    }
    else if (stencilSize==7)
    {
      for (int dof=0; dof<4; ++dof){
	std::array<scalar_type,7> mys;
	int start = dof;
	for (int k=0; k<7; ++k){
	  mys[k] = stencilVals_[start + 4*k];
	}
	pressiodemoapps::weno5(uMinusHalfNeg[dof], uMinusHalfPos[dof],
			       uPlusHalfNeg[dof],  uPlusHalfPos[dof],
			       mys);
      }
    }
  }

  void velocity(const state_type & U,
		const scalar_type t,
		velocity_type & V) const
  {
    using arr_t = std::array<scalar_type, numDofPerCell_>;
    arr_t FL = {0,0,0,0};
    arr_t FR = {0,0,0,0};
    arr_t FU = {0,0,0,0};
    arr_t FD = {0,0,0,0};
    arr_t uMinusHalfNeg = {0,0,0,0};
    arr_t uMinusHalfPos = {0,0,0,0};
    arr_t uPlusHalfNeg  = {0,0,0,0};
    arr_t uPlusHalfPos  = {0,0,0,0};

    const auto dx  = meshObj_.dx();
    const auto dxInv  = meshObj_.dxInv();
    const auto dyInv  = meshObj_.dyInv();
    const auto & graph = meshObj_.graph();
    const auto sampleMeshSize = meshObj_.sampleMeshSize();
    const auto stencilSize = meshObj_.stencilSize();

    fillGhosts(U, t);

    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
      {
	const auto thisCellAdList = graph.row(smPt);
	const auto & cellGID = thisCellAdList(0);
	// given a cell, compute index in V of the correspond first dof
	const auto vIndex = smPt*numDofPerCell_;
	// given a cell, compute index in u of the correspond first dof
	const auto uIndex = cellGID*numDofPerCell_;

	if (stencilSize==3){
	  fillStencil3(smPt, uIndex, thisCellAdList, U,
		       ghostLeft_,   ghostRight_, stencilValsX_, 1);
	  fillStencil3(smPt, uIndex, thisCellAdList, U,
		       ghostBottom_, ghostTop_, stencilValsY_, 2);
	}else{
	  fillStencil7(smPt, uIndex, thisCellAdList, U,
		       ghostLeft_,   ghostRight_, stencilValsX_, 1);
	  fillStencil7(smPt, uIndex, thisCellAdList, U,
		       ghostBottom_, ghostTop_, stencilValsY_, 2);
	}

	//*** X direction ***
	reconstruct(stencilValsX_,
		    uMinusHalfNeg, uMinusHalfPos,
		    uPlusHalfNeg, uPlusHalfPos);
	eeRusanovFlux(FL, uMinusHalfNeg, uMinusHalfPos, normalX_, gamma_);
	eeRusanovFlux(FR, uPlusHalfNeg,  uPlusHalfPos,  normalX_, gamma_);

	//*** Y direction ***
	reconstruct(stencilValsY_,
		    uMinusHalfNeg, uMinusHalfPos,
		    uPlusHalfNeg, uPlusHalfPos);
	eeRusanovFlux(FD, uMinusHalfNeg, uMinusHalfPos, normalY_, gamma_);
	eeRusanovFlux(FU, uPlusHalfNeg,  uPlusHalfPos,  normalY_, gamma_);

	V(vIndex)   = dxInv*(FL[0] - FR[0]) + dyInv*(FD[0] - FU[0]);
	V(vIndex+1) = dxInv*(FL[1] - FR[1]) + dyInv*(FD[1] - FU[1]);
	V(vIndex+2) = dxInv*(FL[2] - FR[2]) + dyInv*(FD[2] - FU[2]);
	V(vIndex+3) = dxInv*(FL[3] - FR[3]) + dyInv*(FD[3] - FU[3]);
      }
  }


  template<int _probid = problemId>
  typename std::enable_if<_probid==0>::type
  fillGhosts(const state_type & U, scalar_type t) const
  {
    // noop
  }

  /*
    for Sedov and Riemann, use neumann conditions
    Currently this is done using zeroth-order, but later one
    we should change to using higher-order (if it makes sense)
   */
  template<int _probid = problemId>
  typename std::enable_if<_probid==1 or _probid==2>::type
  fillGhosts(const state_type & U, scalar_type t) const
  {
    std::array<scalar_type, 4> prim = {0,0,0,0};

    const auto dx  = meshObj_.dx();
    const auto dy  = meshObj_.dy();
    const auto x  = meshObj_.viewX();
    const auto y  = meshObj_.viewY();
    const auto & graph = meshObj_.graph();
    const auto sampleMeshSize = meshObj_.sampleMeshSize();
    const auto stencilSize = meshObj_.stencilSize();

    for (index_t smPt=0; smPt < sampleMeshSize; ++smPt)
    {
      const auto & thisCellAdList = graph.row(smPt);
      const auto & cellGID = thisCellAdList(0);
      const auto myX = x(cellGID);
      const auto myY = y(cellGID);

      const auto w0 = thisCellAdList[1];
      const auto n0 = thisCellAdList[2];
      const auto e0 = thisCellAdList[3];
      const auto s0 = thisCellAdList[4];

      if (w0 == -1){
	const auto ind = cellGID*numDofPerCell_;
	ghostLeft_(smPt, 0) = U[ind];
	ghostLeft_(smPt, 1) = U[ind+1];
	ghostLeft_(smPt, 2) = U[ind+2];
	ghostLeft_(smPt, 3) = U[ind+3];
      }

      if (n0 == -1){
	const auto ind = cellGID*numDofPerCell_;
	ghostTop_(smPt, 0) = U[ind];
	ghostTop_(smPt, 1) = U[ind+1];
	ghostTop_(smPt, 2) = U[ind+2];
	ghostTop_(smPt, 3) = U[ind+3];
      }

      if (e0 == -1){
	const auto ind = cellGID*numDofPerCell_;
	ghostRight_(smPt, 0) = U[ind];
	ghostRight_(smPt, 1) = U[ind+1];
	ghostRight_(smPt, 2) = U[ind+2];
	ghostRight_(smPt, 3) = U[ind+3];
      }

      if (s0 == -1){
	const auto ind = cellGID*numDofPerCell_;
	ghostBottom_(smPt, 0) = U[ind];
	ghostBottom_(smPt, 1) = U[ind+1];
	ghostBottom_(smPt, 2) = U[ind+2];
	ghostBottom_(smPt, 3) = U[ind+3];
      }

      if (stencilSize == 7)
      {
	const auto w1 = thisCellAdList[5];
	const auto n1 = thisCellAdList[6];
	const auto e1 = thisCellAdList[7];
	const auto s1 = thisCellAdList[8];
	const auto w2 = thisCellAdList[9];
	const auto n2 = thisCellAdList[10];
	const auto e2 = thisCellAdList[11];
	const auto s2 = thisCellAdList[12];

	if (w1 == -1){
	  const auto ind = e0*numDofPerCell_;
	  ghostLeft_(smPt, 4) = U[ind];
	  ghostLeft_(smPt, 5) = U[ind+1];
	  ghostLeft_(smPt, 6) = U[ind+2];
	  ghostLeft_(smPt, 7) = U[ind+3];
	}

	if (w2 == -1){
	  const auto ind = e1*numDofPerCell_;
	  ghostLeft_(smPt, 8)  = U[ind];
	  ghostLeft_(smPt, 9)  = U[ind+1];
	  ghostLeft_(smPt, 10) = U[ind+2];
	  ghostLeft_(smPt, 11) = U[ind+3];
	}

	if (n1 == -1){
	  const auto ind = s0*numDofPerCell_;
	  ghostTop_(smPt, 4) = U[ind];
	  ghostTop_(smPt, 5) = U[ind+1];
	  ghostTop_(smPt, 6) = U[ind+2];
	  ghostTop_(smPt, 7) = U[ind+3];
	}

	if (n2 == -1){
	  const auto ind = s1*numDofPerCell_;
	  ghostTop_(smPt, 8) = U[ind];
	  ghostTop_(smPt, 9) = U[ind+1];
	  ghostTop_(smPt, 10) = U[ind+2];
	  ghostTop_(smPt, 11) = U[ind+3];
	}

	if (e1 == -1){
	  const auto ind = w0*numDofPerCell_;
	  ghostRight_(smPt, 4) = U[ind];
	  ghostRight_(smPt, 5) = U[ind+1];
	  ghostRight_(smPt, 6) = U[ind+2];
	  ghostRight_(smPt, 7) = U[ind+3];
	}

	if (e2 == -1){
	  const auto ind = w1*numDofPerCell_;
	  ghostRight_(smPt, 8) = U[ind];
	  ghostRight_(smPt, 9) = U[ind+1];
	  ghostRight_(smPt, 10) = U[ind+2];
	  ghostRight_(smPt, 11) = U[ind+3];
	}

	if (s1 == -1){
	  const auto ind = n0*numDofPerCell_;
	  ghostBottom_(smPt, 4) = U[ind];
	  ghostBottom_(smPt, 5) = U[ind+1];
	  ghostBottom_(smPt, 6) = U[ind+2];
	  ghostBottom_(smPt, 7) = U[ind+3];
	}

	if (s2 == -1){
	  const auto ind = n1*numDofPerCell_;
	  ghostBottom_(smPt, 8) = U[ind];
	  ghostBottom_(smPt, 9) = U[ind+1];
	  ghostBottom_(smPt, 10) = U[ind+2];
	  ghostBottom_(smPt, 11) = U[ind+3];
	}

      }
    }
  }

private:
  void allocateGhosts()
  {
    /*
      for stencil = 3, at leftBoundary:
      ---------------
      |	 0,1,2,3   ||
      | rho,       ||
      | rho u,	   ||
      | rho v,     ||
      | E	   ||
      ---------------

      for stencil = 7, at leftBoundary:
      --------------------------------------
      |	 8,9,10,11  | 4,5,6,7 |  0,1,2,3  ||
      |	     	    |
      | rho,        | rho,    | rho       ||
      | rho u,	    | rho*u   | rho*u     ||
      | rho v,      | rho*v   | rho*v     ||
      | E	    | E       | E         ||
      --------------------------------------

     */

    const auto stencilSize = meshObj_.stencilSize();
    const auto sampleMeshSize = meshObj_.sampleMeshSize();
    const auto numGhostValues = numDofPerCell_*((stencilSize-1)/2);
    ghostLeft_.resize(sampleMeshSize,   numGhostValues);
    ghostTop_.resize(sampleMeshSize,    numGhostValues);
    ghostRight_.resize(sampleMeshSize,  numGhostValues);
    ghostBottom_.resize(sampleMeshSize, numGhostValues);
  }

private:
  const scalar_type gamma_ = (5.+2.)/5.;
  const scalar_type gammaMinusOne_ = gamma_ - 1.;
  const scalar_type gammaMinusOneInv_ = 1./(gamma_ - 1.);
  const scalar_type gammaMinusOneDiv16_ = gammaMinusOne_/(8. * gamma_ * M_PI * M_PI);

  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points. for this problem,
  // the dof = numDofPerCell_ * number_of_unknown_grid_points
  // SampleMesh_ identifies the velocity/residual locations
  index_t numDofStencilMesh_ = {};
  index_t numDofSampleMesh_  = {};

  const std::array<scalar_type, 2> normalX_{1, 0};
  const std::array<scalar_type, 2> normalY_{0, 1};

  mutable ghost_t ghostLeft_;
  mutable ghost_t ghostTop_;
  mutable ghost_t ghostRight_;
  mutable ghost_t ghostBottom_;

  mutable std::vector<scalar_type> stencilValsX_;
  mutable std::vector<scalar_type> stencilValsY_;

  const mesh_t & meshObj_;

  // which initial condition to use
  int icId_ = 1;
};

}}}//end namespace

#endif
