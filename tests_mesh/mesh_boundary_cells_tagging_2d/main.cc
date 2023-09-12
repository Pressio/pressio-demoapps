#include "pressiodemoapps/mesh.hpp"

int main()
{
  namespace pda = pressiodemoapps;

  const std::vector<int> goldCellsGIDsNearBoundariesS3({0,1,2,3,4,5,6,7,
      8,15,16,23,24,31,32,39,40,47,48,55,56,63,64,65,66,67,68,69,70,71});

  const std::vector<int> goldCellsGIDsNearBoundariesS5({0,1,2,3,4,5,6,7,
      8,9,10,11,12,13,14,15,16,17,22,23,24,25,30,31,32,33,
      38,39,40,41,46,47,48,49,54,55,56,57,58,59,60,61,62,
      63,64,65,66,67,68,69,70,71});

  const std::vector<int> goldCellsGIDsNearBoundariesS7({0,1,2,3,4,5,6,7,
      8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,
      30,31,32,33,34,37,38,39,40,41,42,45,46,47,48,49,50,51,52,53,54,55,
      56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71});

  const std::vector<int> goldCellsGIDsStrictlyAtBoundaries({0,1,2,3,4,5,6,7,
      8,15,16,23,24,31,32,39,40,47,48,55,56,63,64,65,66,67,68,69,70,71});

  for (int ss=3; ss<=7; ss+=2){
    const std::string meshPath = "./fullmesh_s" + std::to_string(ss);
    const auto mesh = pda::load_cellcentered_uniform_mesh_eigen(meshPath);
    const auto & G = mesh.graph();

    const auto & rowsCellsStrictlyOnBd = mesh.graphRowsOfCellsStrictlyOnBd();
    std::vector<int> strictlyOnBdGIDs;
    for (auto rowInd : rowsCellsStrictlyOnBd){ strictlyOnBdGIDs.push_back(G(rowInd, 0)); }
    const bool flag1 = std::equal(strictlyOnBdGIDs.begin(), strictlyOnBdGIDs.end(),
				  goldCellsGIDsStrictlyAtBoundaries.begin());

    const auto & rowsCellsNearBd = mesh.graphRowsOfCellsNearBd();
    std::vector<int> nearBdGIDs;
    for (auto rowInd : rowsCellsNearBd){ nearBdGIDs.push_back(G(rowInd, 0)); }
    auto & goldNearBd = (ss==3) ? goldCellsGIDsNearBoundariesS3 :
      (ss==5) ? goldCellsGIDsNearBoundariesS5 : goldCellsGIDsNearBoundariesS7;
    const bool flag2 = std::equal(nearBdGIDs.begin(), nearBdGIDs.end(),
				  goldNearBd.begin());

    if (flag1 && flag2){ std::puts("PASS"); }
    else{ std::puts("FAILED"); }
  }

  return 0;
}
