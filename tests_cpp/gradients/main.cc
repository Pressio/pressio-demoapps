#include <vector>
#include "pressiodemoapps/mesh.hpp"
#include "pressiodemoapps/gradient.hpp"
#include <iomanip>

// constexpr int grad_x = 1;
// constexpr int grad_y = 2;

// struct Face{
//   int parentCellGID;
//   int parentCellGraphRow;
//   double x;
//   double y;
//   double gradValue;
//   int gradDirection;
// };

// template<class MeshType>
// auto one_sided_fd_normal_gradient_at_boundary_faces(const MeshType & mesh,
// 						    const Eigen::VectorXd & f,
// 						    int order)
// {
//   namespace pda  = pressiodemoapps;
//   const int nDofPerCell = 1;
//   const int dofShift = 0;

//   const int ss = mesh.stencilSize();
//   const auto & x = mesh.viewX();
//   const auto & y = mesh.viewY();
//   const auto dx  = mesh.dx();
//   const auto dy  = mesh.dy();

//   const auto & rowsForLoop = mesh.graphRowsOfCellsStrictlyOnBd();
//   const auto & G = mesh.graph();
//   std::vector<Face> result;
//   for (auto rowInd : rowsForLoop)
//   {
//     const bool bdL = mesh.cellHasLeftFaceOnBoundary2d(rowInd);
//     const bool bdF = mesh.cellHasFrontFaceOnBoundary2d(rowInd);
//     const bool bdR = mesh.cellHasRightFaceOnBoundary2d(rowInd);
//     const bool bdB = mesh.cellHasBackFaceOnBoundary2d(rowInd);
//     std::cout << "row = " << rowInd << " " << G(rowInd,0) << std::endl;

//     const auto graphRow = G.row(rowInd);
//     const int cellGID   = graphRow(0);
//     const auto cellX    = x(cellGID);
//     const auto cellY    = y(cellGID);

//     const auto SkewRight  = (order == 1) ? pda::GradFdMode::ForwardTwoPt  : pda::GradFdMode::ForwardThreePt;
//     const auto SkewLeft   = (order == 1) ? pda::GradFdMode::BackwardTwoPt : pda::GradFdMode::BackwardThreePt;
//     Face face;
//     face.parentCellGID = cellGID;
//     face.parentCellGraphRow = rowInd;

//     if (bdL){
//       face.x = cellX-dx*0.5;
//       face.y = cellY;
//       face.gradValue = pda::impl::face_normal_gradient_for_cell_centered_function_2d(f, graphRow, 'x',
// 										     SkewRight, dx, nDofPerCell, dofShift);
//       face.gradDirection = grad_x;
//       result.push_back(face);
//     }

//     if (bdF){
//       face.x = cellX;
//       face.y = cellY+dy*0.5;
//       face.gradValue = pda::impl::face_normal_gradient_for_cell_centered_function_2d(f, graphRow, 'y',
// 										     SkewLeft, dy, nDofPerCell, dofShift);
//       face.gradDirection = grad_y;
//       result.push_back(face);
//     }

//     if (bdR){
//       face.x = cellX+dx*0.5;
//       face.y = cellY;
//       face.gradValue = pda::impl::face_normal_gradient_for_cell_centered_function_2d(f, graphRow, 'x',
// 										     SkewLeft, dx, nDofPerCell, dofShift);
//       face.gradDirection = grad_x;
//       result.push_back(face);
//     }

//     if (bdB){
//       face.x = cellX;
//       face.y = cellY-dy*0.5;
//       face.gradValue = pda::impl::face_normal_gradient_for_cell_centered_function_2d(f, graphRow, 'y',
// 										     SkewRight, dy, nDofPerCell, dofShift);
//       face.gradDirection = grad_y;
//       result.push_back(face);
//     }
//   }

//   return result;
// }


//   // const auto fd_result1 = one_sided_fd_normal_gradient_at_boundary_faces(mesh, f, 1);
//   // write_result(fd_result1, "fd_results_1.txt");
//   // const auto fd_result2 = one_sided_fd_normal_gradient_at_boundary_faces(mesh, f, 2);
//   // write_result(fd_result2, "fd_results_2.txt");


// void write_result(const std::vector<Face> & M,
// 		  const std::string & filename)
// {
//   std::ofstream file; file.open(filename);
//   for (auto & it : M){
//     file << it.parentCellGID << " "
// 	 << it.x << " " << it.y << " "
// 	 << it.gradValue << " "
// 	 << it.gradDirection
// 	 << std::endl;
//   }
//   file.close();
// };

template<class T, class StateType>
void sinef(const T & x, const T & y, StateType & f){
  for (int i=0; i<x.size(); ++i){
    f(i) = std::sin(M_PI * x(i) * y(i));
  }
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
std::pair<double, double> do_test_A(const MeshType& mesh, const std::string & ss)
{
  namespace pda  = pressiodemoapps;
  std::cout << "Doing gradient test for " << ss << '\n';

  // function to compute grad of
  const auto & x = mesh.viewX();
  const auto & y = mesh.viewY();
  Eigen::VectorXd f(mesh.stencilMeshSize());
  sinef(x, y, f);

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

int main()
{
  namespace pda  = pressiodemoapps;

  std::unordered_map<std::string, std::array<double,2>> goldData =
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

  for (auto ss : {3,5,7}){
    {
      const std::string caseID = "fullmesh_s" + std::to_string(ss);
      auto mesh = pda::load_cellcentered_uniform_mesh_eigen(caseID);
      const auto [eX, eY] = do_test_A(mesh, caseID);
      std::cout << "eX = " << eX << " eY = " << eY << '\n';
      if (std::abs(eX - goldData[caseID][0]) > 1e-6){
	std::puts("FAILED");
	return 0;
      }
      if (std::abs(eY - goldData[caseID][1]) > 1e-6){
	std::puts("FAILED");
	return 0;
      }
    }

    {
      const std::string caseID = "samplemesh_s" + std::to_string(ss);
      auto mesh2 = pda::load_cellcentered_uniform_mesh_eigen(caseID);
      const auto [eX, eY] = do_test_A(mesh2, caseID);
      std::cout << "eX = " << eX << " eY = " << eY << '\n';
      if (std::abs(eX - goldData[caseID][0]) > 1e-6){
	std::puts("FAILED");
	return 0;
      }
      if (std::abs(eY - goldData[caseID][1]) > 1e-6){
	std::puts("FAILED");
	return 0;
      }
    }
  }

  //std::cout << "FINISHED\n";
  // std::puts("FAILED");
  std::puts("PASS");

  return 0;
}
