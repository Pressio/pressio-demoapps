#include <vector>
#include "pressiodemoapps/mesh.hpp"
#include "pressiodemoapps/gradient.hpp"
#include <iomanip>

template<class T>
T sinef(const T & x, const T & y){
  return std::sin(M_PI * x * y);
}

template<class T>
T gradsinex(const T & x, const T & y){
  return y*M_PI*std::cos(M_PI * x * y);
}

template<class T>
T gradsiney(const T & x, const T & y){
  return x*M_PI*std::cos(M_PI * x * y);
}


template<class MeshType>
std::pair<double, double> do_test_api_A(const MeshType& mesh, const std::string & ss)
{
  /* this test the api where we pass a function one dof per cell */

  namespace pda  = pressiodemoapps;
  std::cout << "Gradient test for " << ss << '\n';

  // function to compute grad of
  const auto & x = mesh.viewX();
  const auto & y = mesh.viewY();
  Eigen::VectorXd f(mesh.stencilMeshSize());
  for (int i=0; i<x.size(); ++i){
    f(i) = sinef(x(i), y(i));
  }

  // evaluate gradients
  pda::GradientEvaluator grads(mesh);
  grads(f);

  const auto & G = mesh.graph();
  const auto & rowsCellsOnBD = mesh.graphRowsOfCellsStrictlyOnBd();
  std::ofstream file; file.open("pda_result" + ss + ".txt");

  auto toFile = [&](int cellGID, auto const & faceIn, auto gold){
      file << faceIn.centerCoordinates[0] << " "
	   << faceIn.centerCoordinates[1] << " "
	   << faceIn.normalGradient << " "
	   << gold << " "
	   << faceIn.normalDirection
	   << std::endl;
  };

  double errorGradX = {};
  double errorGradY = {};
  int count = 0;
  for (auto rowInd : rowsCellsOnBD){
    const bool bL = mesh.cellHasLeftFaceOnBoundary2d(rowInd);
    const bool bF = mesh.cellHasFrontFaceOnBoundary2d(rowInd);
    const bool bR = mesh.cellHasRightFaceOnBoundary2d(rowInd);
    const bool bB = mesh.cellHasBackFaceOnBoundary2d(rowInd);
    const int cellGID = G(rowInd, 0);

    if (bL){
      auto & face = grads.queryFace(cellGID, pda::FacePosition::Left);
      const auto gold = gradsinex(face.centerCoordinates[0], face.centerCoordinates[1]);
      toFile(cellGID, face, gold);
      const auto err  = face.normalGradient - gold;
      errorGradX += err*err;
      count++;
    }
    if (bF){
      auto & face = grads.queryFace(cellGID, pda::FacePosition::Front);
      const auto gold = gradsiney(face.centerCoordinates[0], face.centerCoordinates[1]);
      toFile(cellGID, face, gold);
      const auto err  = face.normalGradient - gold;
      errorGradY += err*err;
      count++;
    }
    if (bR){
      auto & face = grads.queryFace(cellGID, pda::FacePosition::Right);
      const auto gold = gradsinex(face.centerCoordinates[0], face.centerCoordinates[1]);
      toFile(cellGID, face, gold);
      const auto err  = face.normalGradient - gold;
      errorGradX += err*err;
      count++;
    }
    if (bB){
      auto & face = grads.queryFace(cellGID, pda::FacePosition::Back);
      const auto gold = gradsiney(face.centerCoordinates[0], face.centerCoordinates[1]);
      toFile(cellGID, face, gold);
      const auto err  = face.normalGradient - gold;
      errorGradY += err*err;
      count++;
    }
  }
  file.close();

  return std::make_pair(std::sqrt(errorGradX/count), std::sqrt(errorGradY/count));
}


template<std::size_t MaxNumDofPerCell, class MeshType>
auto do_test_api_B(const MeshType& mesh,
		   const std::string & ss,
		   int numDofPerCell)
{
  /* test the api where we pass a function that has multiple dof per cell
     and here for simplicity we just mimic this by simply repeating
     the sine function multiple times
  */

  namespace pda  = pressiodemoapps;
  std::cout << "Gradient test for " << ss << ' '
	    << " MaxNumDofPerCell = " << MaxNumDofPerCell << ' '
	    << " numDofPerCell = " << numDofPerCell << '\n';

  // function to compute grad of
  const auto & x = mesh.viewX();
  const auto & y = mesh.viewY();
  Eigen::VectorXd f(mesh.stencilMeshSize()*numDofPerCell);
  for (int i=0; i<x.size(); ++i){
    double value = sinef(x(i), y(i));
    for (int j=0; j<numDofPerCell; ++j){
      int fi = i*numDofPerCell;
      f(fi+j) = value;
    }
  }

  // evaluate gradients
  pda::GradientEvaluator<MeshType, MaxNumDofPerCell> grads(mesh);
  grads(f, numDofPerCell);

  const auto & G = mesh.graph();
  const auto & rowsCellsOnBD = mesh.graphRowsOfCellsStrictlyOnBd();
  std::ofstream file; file.open("pda_result" + ss + ".txt");

  auto toFile = [&](int cellGID, auto const & faceIn, auto gold){
      file << faceIn.centerCoordinates[0] << " "
	   << faceIn.centerCoordinates[1] << "\n";
      for (int j=0; j<numDofPerCell; ++j){
	file << faceIn.normalGradient[j] << " ";
      }
      file << gold << " " << faceIn.normalDirection << "\n";
  };

  std::array<double, MaxNumDofPerCell> errorGradX = {};
  std::array<double, MaxNumDofPerCell> errorGradY = {};
  int count = 0;
  for (auto rowInd : rowsCellsOnBD){
    const bool bL = mesh.cellHasLeftFaceOnBoundary2d(rowInd);
    const bool bF = mesh.cellHasFrontFaceOnBoundary2d(rowInd);
    const bool bR = mesh.cellHasRightFaceOnBoundary2d(rowInd);
    const bool bB = mesh.cellHasBackFaceOnBoundary2d(rowInd);
    const int cellGID = G(rowInd, 0);

    if (bL){
      auto & face = grads.queryFace(cellGID, pda::FacePosition::Left);
      const auto gold = gradsinex(face.centerCoordinates[0], face.centerCoordinates[1]);
      toFile(cellGID, face, gold);
      assert(face.normalGradient.size() == MaxNumDofPerCell);
      for (int j=0; j<numDofPerCell; ++j){
	const auto err = face.normalGradient[j] - gold;
	errorGradX[j] += err*err;
      }
      count++;
    }

    if (bF){
      auto & face = grads.queryFace(cellGID, pda::FacePosition::Front);
      const auto gold = gradsiney(face.centerCoordinates[0], face.centerCoordinates[1]);
      toFile(cellGID, face, gold);
      assert(face.normalGradient.size() == MaxNumDofPerCell);
      for (int j=0; j<numDofPerCell; ++j){
	const auto err  = face.normalGradient[j] - gold;
	errorGradY[j] += err*err;
      }
      count++;
    }

    if (bR){
      auto & face = grads.queryFace(cellGID, pda::FacePosition::Right);
      const auto gold = gradsinex(face.centerCoordinates[0], face.centerCoordinates[1]);
      toFile(cellGID, face, gold);
      assert(face.normalGradient.size() == MaxNumDofPerCell);
      for (int j=0; j<numDofPerCell; ++j){
	const auto err  = face.normalGradient[j] - gold;
	errorGradX[j] += err*err;
      }
      count++;
    }

    if (bB){
      auto & face = grads.queryFace(cellGID, pda::FacePosition::Back);
      const auto gold = gradsiney(face.centerCoordinates[0], face.centerCoordinates[1]);
      toFile(cellGID, face, gold);
      assert(face.normalGradient.size() == MaxNumDofPerCell);
      for (int j=0; j<numDofPerCell; ++j){
	const auto err  = face.normalGradient[j] - gold;
	errorGradY[j] += err*err;
      }
      count++;
    }
  }
  file.close();

  for (auto & it : errorGradX){ it = std::sqrt(it/count); }
  for (auto & it : errorGradY){ it = std::sqrt(it/count); }
  return std::make_pair(errorGradX, errorGradY);
}

struct RMSE{
  double rmseGX;
  double rmseGY;
};

bool test_driver_A(const std::unordered_map<std::string, RMSE> & goldData,
		   const std::string & meshstring,
		   int ss)
{
  namespace pda  = pressiodemoapps;
  const std::string caseID = meshstring + "_s" + std::to_string(ss);
  auto mesh = pda::load_cellcentered_uniform_mesh_eigen(caseID);
  const auto [eX, eY] = do_test_api_A(mesh, caseID);
  std::cout << "  eX = " << eX << " eY = " << eY << '\n';
  std::cout << "\n";

  if (std::abs(eX - goldData.at(caseID).rmseGX) > 1e-6){ return false; }
  if (std::abs(eY - goldData.at(caseID).rmseGY) > 1e-6){ return false; }
  return true;
}

template<std::size_t MaxNumDofPerCell>
bool test_driver_B(const std::unordered_map<std::string, RMSE> & goldData,
		   const std::string & meshstring,
		   int ss,
		   int numDofPerCell)
{
  namespace pda  = pressiodemoapps;
  const std::string caseID = meshstring + "_s" + std::to_string(ss);
  auto mesh = pda::load_cellcentered_uniform_mesh_eigen(caseID);
  const auto [eX, eY] = do_test_api_B<MaxNumDofPerCell>(mesh, caseID, numDofPerCell);
  for (int j=0; j<numDofPerCell; ++j){
    std::cout << "  j = " << j
	      << " eX = " << eX[j]
	      << " eY = " << eY[j]
	      << '\n';
    if (std::abs(eX[j] - goldData.at(caseID).rmseGX) > 1e-6){ return false; }
    if (std::abs(eY[j] - goldData.at(caseID).rmseGY) > 1e-6){ return false; }
  }
  std::cout << "\n";

  return true;
}

int main()
{

  std::unordered_map<std::string, RMSE> goldData =
    { {"fullmesh_s3",   {0.118044,  0.082173}},
      {"samplemesh_s3", {0.134392,  0.0759626}},
      {"fullmesh_s5",   {0.0512737, 0.0298508}},
      {"samplemesh_s5", {0.0682185, 0.0307574}},
      {"fullmesh_s7",   {0.0512737, 0.0298508}},
      {"samplemesh_s7", {0.0682185, 0.0307574}}
    };
  // NOTE: the case for s5 and s7 is the same because at the time this test
  // is written, for both cases the stencil used for the gradient is
  // a 3 pt one sided FD which is supported when s5 and s7.
  // Whereas for s3 the FD stencil for gradient is 2pt stencil
  // which, indeed, yields higher error.

  std::vector<bool> flags;
  for (auto ss : {3}){//,5,7}){
    flags.push_back(test_driver_A(goldData, "fullmesh",   ss));
    flags.push_back(test_driver_A(goldData, "samplemesh", ss));

    {
      constexpr int maxNumDofPerCell = 3;
      flags.push_back(test_driver_B<maxNumDofPerCell>(goldData, "fullmesh", ss, 1));
      flags.push_back(test_driver_B<maxNumDofPerCell>(goldData, "fullmesh", ss, 2));
      flags.push_back(test_driver_B<maxNumDofPerCell>(goldData, "fullmesh", ss, 3));
    }

    {
      constexpr int maxNumDofPerCell = 5;
      flags.push_back(test_driver_B<maxNumDofPerCell>(goldData, "fullmesh", ss, 1));
      flags.push_back(test_driver_B<maxNumDofPerCell>(goldData, "fullmesh", ss, 2));
      flags.push_back(test_driver_B<maxNumDofPerCell>(goldData, "fullmesh", ss, 3));
      flags.push_back(test_driver_B<maxNumDofPerCell>(goldData, "fullmesh", ss, 4));
      flags.push_back(test_driver_B<maxNumDofPerCell>(goldData, "fullmesh", ss, 5));
    }
  }

  if (std::any_of(flags.begin(), flags.end(), [](auto v){ return (v==false); })){
    std::puts("FAILED"); return 0;
  }

  std::puts("PASS");
  return 0;
}
