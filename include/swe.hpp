
#ifndef PRESSIO_DEMOAPPS_SWE_HPP_
#define PRESSIO_DEMOAPPS_SWE_HPP_

#include "Eigen/Core"
#include "Eigen/Sparse"
#include "swe_fluxes.hpp"
#include <array>

namespace pressiodemoapps{

class Swe
{
  using eigVec = Eigen::VectorXd;

public:
  using scalar_type	= double;
  using state_type	= eigVec;
  using velocity_type	= eigVec;
  using dense_matrix_type = Eigen::MatrixXd;

  // type to use for all indexing, has to be large enough
  // to support indexing fairly large systems
  using index_t  = int32_t;

  // number of dofs for each cell
  static constexpr index_t numDofPerCell_{3};

  // // type for the adjacency list for each mesh graph
  // // we use an array of 5 indices since all graph nodes have
  // // the same number of connections
  // using node_al_t = std::array<index_t,5>;

  // type to represent connectivity
  using mesh_graph_t  = Eigen::Matrix<index_t,-1,-1,Eigen::RowMajor>;

public:
  template<typename T>
  Swe(const std::string & meshDir,
        const T & params)
    : g_(params[0]),
      mu_ic_(params[1]),
      mu_f_(params[2])
  {
    this->readMesh(meshDir);
  }

  index_t totalDofSampleMesh() const{
    return numDofSampleMesh_;
  }

  index_t totalDofStencilMesh() const{
    return numDofStencilMesh_;
  }

  const eigVec & viewX() const{
    return x_;
  }

  const eigVec & viewY() const{
    return y_;
  }

  // void writeCoordinatesToFile(const std::string & fileout) const{
  //   std::ofstream file; file.open(fileout);
  //   for (index_t i=0; i<(index_t) x_.size(); i++){
  //     file << std::setprecision(14)
	 //   << x_(i) << " " << y_(i)
	 //   << " \n";
  //   }
  //   file.close();
  // }

  state_type getGaussianIC(scalar_type mu) const
  {
    state_type U0(numDofStencilMesh_);
    for (index_t pt=0; pt < numGptStencilMesh_; ++pt)
      {
	const index_t id = pt*numDofPerCell_;
	const auto eX  = pow( x_(pt) - 1.75, 2);
	const auto eY  = pow( y_(pt) - 1.5, 2);
	const auto eX2 = pow( x_(pt) - 2.5, 2);
	const auto eY2 = pow( y_(pt) - 1.25, 2);
	U0(id) = 1. + mu_ic_*exp(-(eX+eY)/0.5) + mu_ic_*exp(-(eX2+eY2)/0.5);
	U0(id+1) = 0.;
	U0(id+2) = 0.;
      }
    return U0;
  }

  velocity_type createVelocity() const {
    velocity_type V(numDofSampleMesh_);
    return V;
  }

  void velocity(const state_type & U,
                const scalar_type /*t*/,
                velocity_type & V) const
  {
    std::array<scalar_type,3> FL = {0,0,0};
    std::array<scalar_type,3> FR = {0,0,0};
    std::array<scalar_type,3> FU = {0,0,0};
    std::array<scalar_type,3> FD = {0,0,0};
    std::array<scalar_type,3> forcing = {0,0,0};

    /* the graph contains the connectivity for all the cells where I need
     * to copute the velocity. I can loop over it. */
    for (index_t smPt=0; smPt < numGptSampleMesh_; ++smPt)
      {
	// ajacency list of this cell
	const auto & thisCellAdList = graph_.row(smPt);
	// gID of this cell
	const auto & cellGID = thisCellAdList(0);
	// given a cell, compute index in V of the correspond first dof
	const auto vIndex = smPt*numDofPerCell_;
	// given a cell, compute index in u of the correspond first dof
	const auto uIndex = cellGID*numDofPerCell_;

	if (stencilSize_ == 3){
	  // the gids of the neighboring cells (we assume 2nd ordered)
	  const auto westCellGid  = thisCellAdList[1];
	  const auto northCellGid = thisCellAdList[2];
	  const auto eastCellGid  = thisCellAdList[3];
	  const auto southCellGid = thisCellAdList[4];

	  uWestIndex_  = westCellGid*numDofPerCell_;
	  uNorthIndex_ = northCellGid*numDofPerCell_;
	  uEastIndex_  = eastCellGid*numDofPerCell_;
	  uSouthIndex_ = southCellGid*numDofPerCell_;
	  sweRusanovFluxFullStateIn(FL, U, uWestIndex_,  uIndex,	normalX_, g_);
	  sweRusanovFluxFullStateIn(FR, U, uIndex,	  uEastIndex_,	normalX_, g_);
	  sweRusanovFluxFullStateIn(FD, U, uSouthIndex_, uIndex,	normalY_, g_);
	  sweRusanovFluxFullStateIn(FU, U, uIndex,	  uNorthIndex_, normalY_, g_);

	  forcing[0] = 0;
	  forcing[1] = mu_f_*U(uIndex+2);
	  forcing[2] = mu_f_*U(uIndex+1);
	  V(vIndex)   = -dxInv_*(FR[0] - FL[0]) - dyInv_*(FU[0] - FD[0]) + forcing[0];
	  V(vIndex+1) = -dxInv_*(FR[1] - FL[1]) - dyInv_*(FU[1] - FD[1]) + forcing[1];
	  V(vIndex+2) = -dxInv_*(FR[2] - FL[2]) - dyInv_*(FU[2] - FD[2]) + forcing[2];
	}
	else if(stencilSize_==7)
	{
	  // const auto w0  = thisCellAdList[1]*numDofPerCell_;
	  // const auto n0  = thisCellAdList[2]*numDofPerCell_;
	  // const auto e0  = thisCellAdList[3]*numDofPerCell_;
	  // const auto s0  = thisCellAdList[4]*numDofPerCell_;
	  // const auto w1  = thisCellAdList[5]*numDofPerCell_;
	  // const auto n1  = thisCellAdList[6]*numDofPerCell_;
	  // const auto e1  = thisCellAdList[7]*numDofPerCell_;
	  // const auto s1  = thisCellAdList[8]*numDofPerCell_;
	  // const auto w2  = thisCellAdList[9]*numDofPerCell_;
	  // const auto n2  = thisCellAdList[10]*numDofPerCell_;
	  // const auto e2  = thisCellAdList[11]*numDofPerCell_;
	  // const auto s2  = thisCellAdList[12]*numDofPerCell_;

	  // // q = (h, hu, hv);

	  // std::array<scalar_type,3> uMinusHalfNeg = {};
	  // std::array<scalar_type,3> uMinusHalfPos = {};
	  // std::array<scalar_type,3> uPlusHalfNeg = {};
	  // std::array<scalar_type,3> uPlusHalfPos = {};

	  // // i-1/2^{-,+}
	  // weno5(uMinusHalfNeg, uMinusHalfPos, U, w2, w1, w0, uIndex, e0, e1);
	  // LF(FL, uMinusHalfNeg, uMinusHalfPos, normalX_, g_);
	  // // i+1/2^{-,+}
	  // weno5(uPlusHalfNeg, uPlusHalfPos, U, w1, w0, uIndex, e0, e1, e2);
	  // LF(FR, uPlusHalfNeg, uPlusHalfPos, normalX_, g_);

	  // // j-1/2^{-,+}
	  // weno5(uMinusHalfNeg, uMinusHalfPos, U, s2, s1, s0, uIndex, n0, n1);
	  // LF(FD, uMinusHalfNeg, uMinusHalfPos, normalY_, g_);
	  // // j+1/2^{-,+}
	  // weno5(uPlusHalfNeg, uPlusHalfPos, U, s1, s0, uIndex, n0, n1, n2);
	  // LF(FU, uPlusHalfNeg, uPlusHalfPos, normalY_, g_);

	  // forcing[0] = 0;
	  // forcing[1] = mu_f_*U(uIndex+2);
	  // forcing[2] = mu_f_*U(uIndex+1);
	  // V(vIndex)   = -dxInv_*(FR[0] - FL[0]) - dyInv_*(FU[0] - FD[0]) + forcing[0];
	  // V(vIndex+1) = -dxInv_*(FR[1] - FL[1]) - dyInv_*(FU[1] - FD[1]) + forcing[1];
	  // V(vIndex+2) = -dxInv_*(FR[2] - FL[2]) - dyInv_*(FU[2] - FD[2]) + forcing[2];
	}
      }
  }

private:
  void readMeshInfo(const std::string & meshDir)
  {
    const auto inFile   = meshDir+"/info.dat";
    std::ifstream foundFile(inFile);
    if(!foundFile){
      std::cout << "file not found " << inFile << std::endl;
      exit(EXIT_FAILURE);
    }

    constexpr auto one = static_cast<scalar_type>(1);
    std::ifstream source( inFile, std::ios_base::in);
    std::string line;
    while (std::getline(source, line) )
      {
	std::istringstream ss(line);
	std::string colVal;
	ss >> colVal;

	if (colVal == "dx"){
	  ss >> colVal;
	  dx_ = std::stod(colVal);
	  dxInv_ = one/dx_;
	  std::cout << "dx = " << dx_ << "\n";
	}

	else if (colVal == "dy"){
	  ss >> colVal;
	  dy_ = std::stod(colVal);
	  dyInv_ = one/dy_;
	  std::cout << "dy = " << dy_ << "\n";
	}

	else if (colVal == "sampleMeshSize"){
	  ss >> colVal;
	  numGptSampleMesh_ = std::stoi(colVal);
	  numDofSampleMesh_ = numGptSampleMesh_ * numDofPerCell_;
	  std::cout << "numGptSampleMesh = " << numGptSampleMesh_ << "\n";
	  std::cout << "numDofSampleMesh = " << numDofSampleMesh_ << "\n";
	}

	else if (colVal == "stencilMeshSize"){
	  ss >> colVal;
	  numGptStencilMesh_ = std::stoi(colVal);
	  numDofStencilMesh_   = numGptStencilMesh_ * numDofPerCell_;
	  std::cout << "numGptStencilMesh = " << numGptStencilMesh_ << "\n";
	  std::cout << "numDofStencilmesh = " << numDofStencilMesh_ << "\n";
	}

	else if (colVal == "stencilSize"){
	  ss >> colVal;
	  stencilSize_ = std::stoi(colVal);
	  std::cout << "stencilSize = " << stencilSize_ << "\n";
	}
      }
    source.close();
  }

  void readMeshCoordinates(const std::string & meshDir)
  {
    const auto inFile   = meshDir+"/coordinates.dat";
    std::ifstream foundFile(inFile);
    if(!foundFile){
      std::cout << "file not found " << inFile << "\n";
      exit(EXIT_FAILURE);
    }

    x_.resize(numGptStencilMesh_);
    y_.resize(numGptStencilMesh_);
    std::ifstream source( inFile, std::ios_base::in);
    std::string line;
    while (std::getline(source, line) )
      {
	std::istringstream ss(line);
	std::string colVal;
	ss >> colVal;
	auto thisGid = std::stoi(colVal);
	ss >> colVal; x_[thisGid] = std::stod(colVal);
	ss >> colVal; y_[thisGid] = std::stod(colVal);
      }
    source.close();
  }

  void readMeshConnectivity(const std::string & meshDir)
  {
    const auto inFile   = meshDir+"/connectivity.dat";
    std::ifstream foundFile(inFile);
    if(!foundFile){
      std::cout << "file not found " << inFile << "\n";
      exit(EXIT_FAILURE);
    }

    const auto graphSize = (stencilSize_-1)*2;
    graph_.resize(numGptSampleMesh_, graphSize+1);
    std::ifstream source(inFile, std::ios_base::in);
    std::string line;
    index_t count = 0;
    while (std::getline(source, line) )
      {
	std::istringstream ss(line);
	std::string colVal;
	ss >> colVal;
	auto thisGid = std::stoi(colVal);
	graph_(count, 0) = thisGid;
	// store others on same row
	for (auto i=1; i<=graph_.cols()-1; ++i){
	  ss >> colVal; thisGid = std::stoi(colVal);
	  graph_(count,i) = thisGid;
	}
	count++;
      }
    source.close();
  }

  void readMesh(const std::string & meshDir)
  {
    readMeshInfo(meshDir);
    readMeshCoordinates(meshDir);
    readMeshConnectivity(meshDir);
  }

private:
  scalar_type dx_{};
  scalar_type dy_{};
  scalar_type dxInv_{};
  scalar_type dyInv_{};

  /*
    graph: contains a list such that
    1 0 3 2 -1

    first col   : contains GIDs of cells where we want velocity
    col 1,2,3,4 : contains GIDs of neighboring cells needed for stencil
    the order of the neighbors is: west, north, east, south
    col 4,5,6,7 : contains GIDs of degree2 neighboring cells needed for stencil
    the order of the neighbors is: west, north, east, south
  */
  int stencilSize_ = {};
  mesh_graph_t graph_ = {};

  mutable index_t uWestIndex_  = {};
  mutable index_t uNorthIndex_ = {};
  mutable index_t uEastIndex_  = {};
  mutable index_t uSouthIndex_ = {};

  // note that dof refers to the degress of freedom,
  // which is NOT same as grid points. for this problem,
  // the dof = numDofPerCell_ * number_of_unknown_grid_points
  // SampleMesh_ identifies the velocity/residual locations
  index_t numGptStencilMesh_ = {};
  index_t numDofStencilMesh_ = {};
  index_t numGptSampleMesh_  = {};
  index_t numDofSampleMesh_  = {};

  // x,y define the coords of all cells centers
  eigVec x_ = {};
  eigVec y_ = {};

  const std::array<scalar_type, 2> normalX_{1, 0};
  const std::array<scalar_type, 2> normalY_{0, 1};

  scalar_type g_ = 6.;
  scalar_type mu_ic_ = 0.125;
  scalar_type mu_f_ = 0.2;
};

}//end namespace

#endif
