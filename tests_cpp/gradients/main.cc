#include <vector>
#include "pressiodemoapps/mesh.hpp"
#include <iomanip>

    // double x0 = cellX;
    // double x1 = cellX + dx;
    // double y1 = cellY;
    // double f0 = f(row[0]);
    // double f1 = f(row[3]);
    // double alfa = -(f0-f1)/dx;
    // double beta = f0 - alfa*x0;
    // double fleft  = alfa * bdFaceX + beta;
    // gradVals.push_back( (-fleft + f1)/(2.*dx) );

template<class MeshType>
bool cell_has_left_face_on_boundary(const MeshType & mesh,
				    int graphRow){
  // true if the first order neighboring cell index is -1
  const auto & G = mesh.graph();
  return (G(graphRow,1)==-1);
}

template<class MeshType>
bool cell_has_front_face_on_boundary(const MeshType & mesh,
				     int graphRow){
  const auto & G = mesh.graph();
  return (G(graphRow,2)==-1);
}

template<class MeshType>
bool cell_has_right_face_on_boundary(const MeshType & mesh,
				     int graphRow){
  const auto & G = mesh.graph();
  return (G(graphRow,3)==-1);
}

template<class MeshType>
bool cell_has_back_face_on_boundary(const MeshType & mesh,
				     int graphRow){
  const auto & G = mesh.graph();
  return (G(graphRow,4)==-1);
}

template<class MeshType>
bool cell_is_strictly_next_to_boundary(const MeshType & mesh,
				       int graphRow)
{
  const auto bL = cell_has_left_face_on_boundary(mesh, graphRow);
  const auto bF = cell_has_front_face_on_boundary(mesh, graphRow);
  const auto bR = cell_has_right_face_on_boundary(mesh, graphRow);
  const auto bB = cell_has_back_face_on_boundary(mesh, graphRow);
  return bL || bF || bR || bB;
}

template<class T, class StateType>
void sinef(const T & x, const T & y, StateType & f)
{
  for (int i=0; i<x.size(); ++i){
    f(i) = std::sin(M_PI * x(i) * y(i));
  }
}

enum class Skew{Right2, Left2, Center3, Right3, Left3, Center5};

template<class Connect>
double fd_gradient_x(Skew direction,
		     int stencilSize,
		     const Connect & conn,
		     const double h,
		     const Eigen::VectorXd & f)
{
  if (stencilSize != 3 && stencilSize!=5){
    throw std::runtime_error("invalid stencil size");
  }

  // https://web.media.mit.edu/~crtaylor/calculator.html
  switch(direction){
  case Skew::Right2:
    return (-f(conn[0]) + f(conn[3]))/h;
  case Skew::Left2:
    return ( f(conn[0]) - f(conn[1]))/h;
  case Skew::Center3:
    return (-f(conn[1]) + f(conn[3]))/(2*h);

  case Skew::Right3:
    return (-2.*f(conn[0]) + 3.*f(conn[3]) - 1.*f(conn[7]))/h;
  case Skew::Left3:
    return ( 2.*f(conn[0]) - 3.*f(conn[1]) + 1.*f(conn[5]))/h;
  case Skew::Center5:
    return (-f(conn[5]) - 8.*f(conn[1]) + 8.*f(conn[3]) - f(conn[7]))/(12.*h);

  default:
    return 0;
  }
}

template<class Connect>
double fd_gradient_y(Skew direction,
		     int stencilSize,
		     const Connect & conn,
		     const double h,
		     const Eigen::VectorXd & f)
{
  if (stencilSize != 3 && stencilSize!=5){
    throw std::runtime_error("invalid stencil size");
  }

  // https://web.media.mit.edu/~crtaylor/calculator.html
  switch(direction){
  case Skew::Right2:
    return (-f(conn[0]) + f(conn[2]))/h;
  case Skew::Left2:
    return ( f(conn[0]) - f(conn[4]))/h;
  case Skew::Center3:
    return (-f(conn[4]) + f(conn[2]))/(2*h);

  case Skew::Right3:
    return (-2.*f(conn[0]) + 3.*f(conn[2]) - 1.*f(conn[6]))/h;
  case Skew::Left3:
    return ( 2.*f(conn[0]) - 3.*f(conn[4]) + 1.*f(conn[8]))/h;
  case Skew::Center5:
    return (-f(conn[8]) - 8.*f(conn[4]) + 8.*f(conn[2]) - f(conn[6]))/(12.*h);

  default:
    return 0;
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

  const int ss = mesh.stencilSize();
  const auto & x = mesh.viewX();
  const auto & y = mesh.viewY();
  const auto dx  = mesh.dx();
  const auto dy  = mesh.dy();

  const auto & G = mesh.graph();
  std::vector<int> rowsForLoop;
  for (auto it : mesh.graphRowsOfCellsNearBd()){
    const bool b = cell_is_strictly_next_to_boundary(mesh, it);
    if (b) rowsForLoop.push_back(it);
  }

  std::vector<Face> result;
  for (auto rowInd : rowsForLoop)
  {
    const bool bdL = cell_has_left_face_on_boundary(mesh,  rowInd);
    const bool bdF = cell_has_front_face_on_boundary(mesh, rowInd);
    const bool bdR = cell_has_right_face_on_boundary(mesh, rowInd);
    const bool bdB = cell_has_back_face_on_boundary(mesh,  rowInd);
    std::cout << "row = " << rowInd << " " << G(rowInd,0) << std::endl;

    const auto graphRow = G.row(rowInd);
    const int cellGID   = graphRow(0);
    const auto cellX    = x(cellGID);
    const auto cellY    = y(cellGID);

    const auto SkewRight  = (order == 1) ? Skew::Right2  : Skew::Right3;
    const auto SkewLeft   = (order == 1) ? Skew::Left2   : Skew::Left3;

    Face face;
    face.parentCellGID = cellGID;
    face.parentCellGraphRow = rowInd;

    if (bdL){
      face.x = cellX-dx*0.5;
      face.y = cellY;
      face.gradValue = fd_gradient_x(SkewRight, ss, graphRow, dx, f);
      face.gradDirection = grad_x;
      result.push_back(face);

      if (bdB){
	face.x = cellX;
	face.y = cellY-dy*0.5;
	face.gradValue = fd_gradient_y(SkewRight, ss, graphRow, dy, f);
	face.gradDirection = grad_y;
	result.push_back(face);
      }

      if (bdF){
	face.x = cellX;
	face.y = cellY+dy*0.5;
	face.gradValue = fd_gradient_y(SkewLeft, ss, graphRow, dy, f);
	face.gradDirection = grad_y;
	result.push_back(face);
      }
    }

    if (bdF && !bdL && !bdR){
      face.x = cellX;
      face.y = cellY+dy*0.5;
      face.gradValue = fd_gradient_y(SkewLeft, ss, graphRow, dy, f);
      face.gradDirection = grad_y;
      result.push_back(face);
    }

    if (bdR){
      face.x = cellX+dx*0.5;
      face.y = cellY;
      face.gradValue = fd_gradient_x(SkewLeft, ss, graphRow, dx, f);
      face.gradDirection = grad_x;
      result.push_back(face);

      if (bdB){
	face.x = cellX;
	face.y = cellY-dy*0.5;
	face.gradValue = fd_gradient_y(SkewRight, ss, graphRow, dy, f);
	face.gradDirection = grad_y;
	result.push_back(face);
      }

      if (bdF){
	face.x = cellX;
	face.y = cellY+dy*0.5;
	face.gradValue = fd_gradient_y(SkewLeft, ss, graphRow, dy, f);
	face.gradDirection = grad_y;
	result.push_back(face);
      }
    }

    if (bdB && !bdL && !bdR){
      face.x = cellX;
      face.y = cellY-dy*0.5;
      face.gradValue = fd_gradient_y(SkewRight, ss, graphRow, dy, f);
      face.gradDirection = grad_y;
      result.push_back(face);
    }
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
};
//file << std::setprecision(14);

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

  const auto fd_result2 = one_sided_fd_normal_gradient_at_boundary_faces(mesh, f, 2);
  write_result(fd_result2, "fd_results_2.txt");

  // std::vector< std::array<double, 3> > fd_grads;
  // const auto & G = mesh.graph();
  // for (int i=0; i<G.rows(); ++i){
  //   one_sided_gradient_fd(i, mesh, fd_grads);
  //   for (int j=1; j<G.cols(); ++j){
  //     std::cout << G(i,j ) << " ";
  //   }
  //   std::cout << '\n';
  // }

  // write_vector_to_ascii(std::string("grad_x.txt"), gradXs);
  // write_vector_to_ascii(std::string("grad_y.txt"), gradYs);
  // write_vector_to_ascii(std::string("grad_v.txt"), gradVals);

  // std::puts("FAILED");
  // std::puts("PASS");

  return 0;
}
