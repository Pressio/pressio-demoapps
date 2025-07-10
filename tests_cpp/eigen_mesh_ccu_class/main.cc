
#include "pressiodemoapps/mesh.hpp"
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
  0.50000000000000, 0.70000000000000 ,
  0.90000000000000
};

std::vector<double> goldY = {
  0.1, 0.1, 0.1, 0.1, 0.1,
  0.3, 0.3, 0.3, 0.3, 0.3,
  0.5, 0.5, 0.5, 0.5, 0.5,
  0.7, 0.7, 0.7, 0.7, 0.7,
  0.9, 0.9, 0.9, 0.9, 0.9
};

auto test_default_constr()
{
  namespace pda = pressiodemoapps;
  pda::cellcentered_uniform_mesh_eigen_type meshObj;

  std::vector<bool> v;
  v.push_back( meshObj.dimensionality() == 0 );
  v.push_back( meshObj.stencilMeshSize() == 0 );
  v.push_back( meshObj.sampleMeshSize() == 0 );
  v.push_back( meshObj.stencilSize() == 0 );
  v.push_back( meshObj.graph().rows() == 0 );
  v.push_back( meshObj.graph().cols() == 0 );
  v.push_back( meshObj.dx() == 0 );
  v.push_back( meshObj.dxInv() == 0 );
  v.push_back( meshObj.dy() == 0 );
  v.push_back( meshObj.dyInv() == 0 );
  v.push_back( meshObj.dz() == 0 );
  v.push_back( meshObj.dzInv() == 0 );
  v.push_back( meshObj.numCellsInner()==0);
  v.push_back( meshObj.numCellsNearBd()==0);

  const auto g1 = meshObj.graphRowsOfCellsAwayFromBd();
  const auto g2 = meshObj.graphRowsOfCellsNearBd();
  v.push_back( g1.size()==0);
  v.push_back( g2.size()==0);

  const auto & x = meshObj.viewX();
  const auto & y = meshObj.viewY();
  const auto & z = meshObj.viewZ();
  v.push_back( x.size()==0 );
  v.push_back( y.size()==0 );
  v.push_back( z.size()==0 );

  return v;
}

template<class T>
auto check_nontrivial(const T & meshObj)
{
  std::vector<bool> v;
  v.push_back( meshObj.dimensionality() == 2 );
  v.push_back( meshObj.stencilMeshSize() == 25 );
  v.push_back( meshObj.sampleMeshSize() == 25 );
  v.push_back( meshObj.stencilSize() == 3  );
  v.push_back( meshObj.graph().rows() == 25 );
  v.push_back( meshObj.graph().cols() == 5 );
  v.push_back( meshObj.dx() == 0.2 );
  v.push_back( meshObj.dxInv() == 5. );
  v.push_back( meshObj.dy() == 0.2 );
  v.push_back( meshObj.dyInv() == 5. );
  v.push_back( meshObj.dz() == 0 );

  const auto & x = meshObj.viewX();
  const auto & y = meshObj.viewY();
  const auto & z = meshObj.viewZ();
  for (int i=0; i<x.size(); ++i){
    const auto ex = std::abs( x(i) - goldX[i] );
    const auto ey = std::abs( y(i) - goldY[i] );
    const auto ez = std::abs( z(i) - 0.0 );
    if (ex > 1e-13 or ey > 1e-13 or ez > 1e-13) {
      v.push_back(false);
      return v;
    }
  }
  return v;
}

auto test_nontrivial()
{
  const auto meshObj = pressiodemoapps::load_cellcentered_uniform_mesh_eigen(".");
  return check_nontrivial(meshObj);
}

auto test_assign()
{
  namespace pda = pressiodemoapps;
  pda::cellcentered_uniform_mesh_eigen_type meshObj;
  const auto meshObj2 = pressiodemoapps::load_cellcentered_uniform_mesh_eigen(".");
  meshObj = meshObj2;
  return check_nontrivial(meshObj);
}

auto test_move()
{
  namespace pda = pressiodemoapps;
  pda::cellcentered_uniform_mesh_eigen_type meshObj;
  const auto meshObj2 = pressiodemoapps::load_cellcentered_uniform_mesh_eigen(".");
  meshObj = std::move(meshObj2);
  return check_nontrivial(meshObj);
}

bool should_pass(std::vector<bool> const & r){
  return std::any_of(r.cbegin(), r.cend(), [](bool v){ return v==false;} );
}

int main()
{
  if ( should_pass(test_default_constr()) ){
    std::puts("FAILED");
    return 0;
  }

  if ( should_pass(test_nontrivial()) ){
    std::puts("FAILED");
    return 0;
  }

  if ( should_pass(test_assign()) ){
    std::puts("FAILED");
    return 0;
  }

  if ( should_pass(test_move()) ){
    std::puts("FAILED");
    return 0;
  }

  std::puts("PASS");
  return 0;
}
