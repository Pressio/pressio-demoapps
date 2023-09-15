#include <vector>
#include "pressiodemoapps/mesh.hpp"
#include "pressiodemoapps/gradient.hpp"
#include <iomanip>

template<class T, class StateType>
void sinef(const T & x, const T & y, StateType & f)
{
  for (int i=0; i<x.size(); ++i){
    f(i) = std::sin(M_PI * x(i) * y(i));
  }
}

constexpr int grad_x = 0;
constexpr int grad_y = 1;

struct Face{
  int parentCellGID;
  int parentCellGraphRow;
  double x;
  double y;
  double gradValue;
  int gradDirection;
};

template<class MeshType>
auto one_sided_fd_normal_gradient_at_boundary_faces(const MeshType & mesh,
						    const Eigen::VectorXd & f,
						    int order)
{
  namespace pda  = pressiodemoapps;

  const int ss = mesh.stencilSize();
  const auto & x = mesh.viewX();
  const auto & y = mesh.viewY();
  const auto dx  = mesh.dx();
  const auto dy  = mesh.dy();

  const auto & rowsForLoop = mesh.graphRowsOfCellsStrictlyOnBd();
  const auto & G = mesh.graph();
  std::vector<Face> result;
  for (auto rowInd : rowsForLoop)
  {
    const bool bdL = mesh.cellHasLeftFaceOnBoundary2d(rowInd);
    const bool bdF = mesh.cellHasFrontFaceOnBoundary2d(rowInd);
    const bool bdR = mesh.cellHasRightFaceOnBoundary2d(rowInd);
    const bool bdB = mesh.cellHasBackFaceOnBoundary2d(rowInd);
    std::cout << "row = " << rowInd << " " << G(rowInd,0) << std::endl;

    const auto graphRow = G.row(rowInd);
    const int cellGID   = graphRow(0);
    const auto cellX    = x(cellGID);
    const auto cellY    = y(cellGID);

    const auto SkewRight  = (order == 1) ? pda::GradFdMode::ForwardTwoPt  : pda::GradFdMode::ForwardThreePt;
    const auto SkewLeft   = (order == 1) ? pda::GradFdMode::BackwardTwoPt : pda::GradFdMode::BackwardThreePt;
    Face face;
    face.parentCellGID = cellGID;
    face.parentCellGraphRow = rowInd;

    if (bdL){
      face.x = cellX-dx*0.5;
      face.y = cellY;
      face.gradValue = pda::impl::fd_face_normal_gradient_for_cell_centered_state_2d(f, graphRow, 'x',
										     SkewRight, dx);
      face.gradDirection = grad_x;
      result.push_back(face);
    }

    // if (bdF){
    //   face.x = cellX;
    //   face.y = cellY+dy*0.5;
    //   face.gradValue = pda::impl::fd_face_normal_gradient_for_cell_centered_state_2d(f, graphRow, 'y',
    // 										     SkewLeft, dy);
    //   face.gradDirection = grad_y;
    //   result.push_back(face);
    // }

    if (bdR){
      face.x = cellX+dx*0.5;
      face.y = cellY;
      face.gradValue = pda::impl::fd_face_normal_gradient_for_cell_centered_state_2d(f, graphRow, 'x',
										     SkewLeft, dx);
      face.gradDirection = grad_x;
      result.push_back(face);
    }

    // if (bdB){
    //   face.x = cellX;
    //   face.y = cellY-dy*0.5;
    //   face.gradValue = pda::impl::fd_face_normal_gradient_for_cell_centered_state_2d(f, graphRow, 'y',
    // 										     SkewRight, dy);
    //   face.gradDirection = grad_y;
    //   result.push_back(face);
    // }
  }

  return result;
}

void write_result(const std::vector<Face> & M,
		  const std::string & filename)
{
  std::ofstream file; file.open(filename);
  for (auto & it : M){
    file << it.parentCellGID << " "
	 << it.parentCellGraphRow << " "
	 << it.x << " " << it.y << " "
	 << it.gradValue << " "
	 << it.gradDirection
	 << std::endl;
  }
  file.close();
};

int main()
{
  namespace pda  = pressiodemoapps;

  auto mesh = pda::load_cellcentered_uniform_mesh_eigen(".");

  const auto & x = mesh.viewX();
  const auto & y = mesh.viewY();
  Eigen::VectorXd f(mesh.stencilMeshSize());
  sinef(x, y, f);

  const auto fd_result1 = one_sided_fd_normal_gradient_at_boundary_faces(mesh, f, 1);
  write_result(fd_result1, "fd_results_1.txt");

  // const auto fd_result2 = one_sided_fd_normal_gradient_at_boundary_faces(mesh, f, 2);
  // write_result(fd_result2, "fd_results_2.txt");

  pda::GradientEvaluator grads(mesh);
  grads.compute(f);
  const auto & G = mesh.graph();
  const auto & rowsCellsOnBD = mesh.graphRowsOfCellsStrictlyOnBd();
  std::ofstream file; file.open("pda_result.txt");
  for (auto rowInd : rowsCellsOnBD){
    const bool bdL = mesh.cellHasLeftFaceOnBoundary2d(rowInd);
    const bool bdF = mesh.cellHasFrontFaceOnBoundary2d(rowInd);
    const bool bdR = mesh.cellHasRightFaceOnBoundary2d(rowInd);
    const bool bdB = mesh.cellHasBackFaceOnBoundary2d(rowInd);
    const int cellGID = G(rowInd, 0);

    pda::FacePosition fp ={};

    // const auto fp = bdL ? pda::FacePosition::Left
    //   : bdF ? pda::FacePosition::Front
    //   : bdR ? pda::FacePosition::Right
    //   : pda::FacePosition::Back;
    int gradDirection = (bdL || bdR ) ? grad_x : grad_y;

    if (bdR){
      auto & f = grads.queryFace(cellGID, fp);
      std::cout << cellGID << " "
		<< f.centerCoords[0] << " "
		<< f.centerCoords[1]
		<< '\n';
      file << cellGID << " "
	   << 0 << " "
	   << f.centerCoords[0] << " "
	   << f.centerCoords[1] << " "
	   << f.normalGradValue[0] << " "
	   << gradDirection
	   << std::endl;
    }
  }
  file.close();

  // std::puts("FAILED");
  // std::puts("PASS");

  return 0;
}
