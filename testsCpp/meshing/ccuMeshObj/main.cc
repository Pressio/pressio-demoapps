
#include "mesh.hpp"
#include <iomanip>
#include <vector>

std::vector<double> goldX = {
  0.10000000000000, 0.30000000000000 ,
  0.50000000000000, 0.70000000000000 ,
  0.90000000000000, 0.10000000000000 ,
  0.30000000000000, 0.50000000000000 ,
  0.70000000000000, 0.90000000000000 ,
  0.10000000000000, 0.30000000000000 ,
  0.50000000000000, 0.70000000000000 ,
  0.90000000000000, 0.10000000000000 ,
  0.30000000000000, 0.50000000000000 ,
  0.70000000000000, 0.90000000000000 ,
  0.10000000000000, 0.30000000000000 ,
  0.50000000000000, 0.70000000000000 , 0.90000000000000
};

std::vector<double> goldY = {
  0.1, 0.1, 0.1, 0.1, 0.1,
  0.3, 0.3, 0.3, 0.3, 0.3,
  0.5, 0.5, 0.5, 0.5, 0.5,
  0.7, 0.7, 0.7, 0.7, 0.7,
  0.9, 0.9, 0.9, 0.9, 0.9
};

int main(int argc, char *argv[])
{
  using scalar_t    = double;
  using mesh_t      = pressiodemoapps::CellCenteredUniformMeshEigen<scalar_t>;
  mesh_t meshO(".");

  std::array<bool, 11> v= {};

  v[0] = meshO.dimensionality() == 2;
  v[1] = meshO.stencilMeshSize() == 25;
  v[2] = meshO.sampleMeshSize() == 25;
  v[3] = meshO.stencilSize() == 3;
  v[4] = meshO.graph().rows() == 25;
  v[5] = meshO.graph().cols() == 5;
  v[6] = meshO.dx() == 0.2;
  v[7] = meshO.dxInv() == 5.;
  v[8] = meshO.dy() == 0.2;
  v[9] = meshO.dyInv() == 5.;

  const auto & x = meshO.viewX();
  const auto & y = meshO.viewY();
  for (int i=0; i<x.size(); ++i){
    const auto ex = std::abs( x(i) - goldX[i] );
    const auto ey = std::abs( y(i) - goldY[i] );
    if (ex > 1e-13 or ey > 1e-13) {
      std::puts("FAILED");
      return 0;
    }
  }

  std::puts("PASS");
  return 0;
}
