#include "pressiodemoapps/mesh.hpp"

int main()
{
  namespace pda = pressiodemoapps;

  // gold for full mesh
  const std::vector<int> goldCellsGIDsNearBoundaries_fullS3({0,1,2,3,4,5,6,7,
      8,15,16,23,24,31,32,39,40,47,48,55,56,63,64,65,66,67,68,69,70,71});

  const std::vector<int> goldCellsGIDsNearBoundaries_fullS5({0,1,2,3,4,5,6,7,
      8,9,10,11,12,13,14,15,16,17,22,23,24,25,30,31,32,33,
      38,39,40,41,46,47,48,49,54,55,56,57,58,59,60,61,62,
      63,64,65,66,67,68,69,70,71});

  const std::vector<int> goldCellsGIDsNearBoundaries_fullS7({0,1,2,3,4,5,6,7,
      8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,29,
      30,31,32,33,34,37,38,39,40,41,42,45,46,47,48,49,50,51,52,53,54,55,
      56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71});

  const std::vector<int> goldCellsGIDsStrictlyAtBoundaries_full({0,1,2,3,4,5,6,7,
      8,15,16,23,24,31,32,39,40,47,48,55,56,63,64,65,66,67,68,69,70,71});

  // gold for sample mesh
  const std::vector<int> goldCellsGIDsNearBoundaries_sampleS3({0,3,6,7,18,22,29,35,40,42,45});
  const std::vector<int> goldCellsGIDsStrictlyAtBoundaries_sampleS3({0,3,6,7,18,22,29,35,40,42,45});

  const std::vector<int> goldCellsGIDsNearBoundaries_sampleS5({0,5,6,11,27,32,46,47,48,53});
  const std::vector<int> goldCellsGIDsStrictlyAtBoundaries_sampleS5({0,5,6,32,46,47,48,53});

  const std::vector<int> goldCellsGIDsNearBoundaries_sampleS7({0,6,7,12,29,37,41,52,53,54,60});
  const std::vector<int> goldCellsGIDsStrictlyAtBoundaries_sampleS7({0,6,7,37,52,53,54,60});

  for (auto sIt : {std::string("fullmesh"), std::string("samplemesh")}){
    for (auto ss : {3, 5, 7}){
      std::cout << sIt << " " << ss << '\n';

      const std::string meshPath = "./" + sIt + "_s" + std::to_string(ss);
      const auto mesh = pda::load_cellcentered_uniform_mesh_eigen(meshPath);
      const auto & G = mesh.graph();

      // first check cells strictly on BD
      const auto & rowsCellsStrictlyOnBd = mesh.graphRowsOfCellsStrictlyOnBd();
      std::vector<int> strictlyOnBdGIDs;
      for (auto rowInd : rowsCellsStrictlyOnBd){ strictlyOnBdGIDs.push_back(G(rowInd, 0)); }
      const std::vector<int> * goldStrictlyOnBd;
      if (sIt == "fullmesh"){
	goldStrictlyOnBd = &goldCellsGIDsStrictlyAtBoundaries_full;
      }
      else{
	goldStrictlyOnBd = (ss==3) ? &goldCellsGIDsStrictlyAtBoundaries_sampleS3 :
	  (ss==5) ? &goldCellsGIDsStrictlyAtBoundaries_sampleS5 : &goldCellsGIDsStrictlyAtBoundaries_sampleS7;
      }
      const bool flag1 = (strictlyOnBdGIDs.size() == goldStrictlyOnBd->size());
      const bool flag2 = std::equal(strictlyOnBdGIDs.begin(), strictlyOnBdGIDs.end(), goldStrictlyOnBd->begin());

      // check cells "near" BD
      const auto & rowsCellsNearBd = mesh.graphRowsOfCellsNearBd();
      std::vector<int> nearBdGIDs;
      for (auto rowInd : rowsCellsNearBd){ nearBdGIDs.push_back(G(rowInd, 0)); }
      std::for_each(nearBdGIDs.begin(), nearBdGIDs.end(), [](auto v){ std::cout << v << " "; });

      const std::vector<int> * goldNearBd;
      if (sIt == "fullmesh"){
	goldNearBd = (ss==3) ? &goldCellsGIDsNearBoundaries_fullS3 :
	  (ss==5) ? &goldCellsGIDsNearBoundaries_fullS5 : &goldCellsGIDsNearBoundaries_fullS7;
      }
      else{
	goldNearBd = (ss==3) ? &goldCellsGIDsNearBoundaries_sampleS3 :
	  (ss==5) ? &goldCellsGIDsNearBoundaries_sampleS5 : &goldCellsGIDsNearBoundaries_sampleS7;
      }
      const bool flag3 = (nearBdGIDs.size() == goldNearBd->size());
      const bool flag4 = std::equal(nearBdGIDs.begin(), nearBdGIDs.end(), goldNearBd->begin());

      std::cout << "flag1 = " << flag1 << " , "
		<< "flag2 = " << flag2 << " , "
		<< "flag3 = " << flag3 << " , "
		<< "flag4 = " << flag4 << '\n';

      if (flag1 && flag2 && flag3 && flag4){ std::puts("PASS"); }
      else{ std::puts("FAILED"); return 0; }
    }
  }

  return 0;
}
