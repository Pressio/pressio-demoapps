
#include "pressiodemoapps/mesh.hpp"
#include <iomanip>
#include <vector>

int main(int argc, char *argv[])
{
  const auto meshObj = pressiodemoapps::loadCellCenterUniformMeshEigen(".");
  using int_t = typename decltype(meshObj)::index_t;

  std::vector<int_t> goldInn = {6,7,8,11,12,13,16,17,18};
  const auto & inn = meshObj.graphRowsOfCellsAwayFromBd();
  if (goldInn != inn ){ std::puts("FAILED"); return 0; }

  std::vector<int_t> goldBd = {0,1,2,3,4,5,9,10,14,15,19,20,21,22,23,24};
  const auto & nearBd = meshObj.graphRowsOfCellsNearBd();
  if (goldBd != nearBd ){ std::puts("FAILED"); return 0; }

  auto s1 = meshObj.numCellsInner();
  auto s2 = meshObj.numCellsBd();
  if (s1!=9 or s2!=16){
  std::puts("FAILED");
  }

  std::puts("PASS");
  return 0;
}
