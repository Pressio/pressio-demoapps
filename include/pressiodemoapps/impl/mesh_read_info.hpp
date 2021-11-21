
#ifndef PRESSIODEMOAPPS_MESH_READ_INFO_HPP_
#define PRESSIODEMOAPPS_MESH_READ_INFO_HPP_

namespace pressiodemoapps{ namespace impl{

template<class sc_t, class int_t>
void read_mesh_info(const std::string & meshDir,
		    int & dim,
		    std::array<sc_t, 3> & cellDeltas,
		    std::array<sc_t, 3> & cellDeltasInv,
		    int & stencilSize,
		    int_t & numGptStencilMesh,
		    int_t & numGptSampleMesh)
{
  const auto inFile   = meshDir+"/info.dat";
  std::ifstream foundFile(inFile);
  if(!foundFile){
    std::cout << "file not found " << inFile << std::endl;
    exit(EXIT_FAILURE);
  }

  constexpr auto one = static_cast<sc_t>(1);
  std::ifstream source( inFile, std::ios_base::in);
  std::string line;
  while (std::getline(source, line) )
    {
      std::istringstream ss(line);
      std::string colVal;
      ss >> colVal;

      if (colVal == "dim"){
	ss >> colVal;
	dim = std::stoi(colVal);
      }

      else if (colVal == "dx"){
	ss >> colVal;
	cellDeltas[0]    = std::stod(colVal);
	cellDeltasInv[0] = one/cellDeltas[0];
      }

      else if (colVal == "dy"){
	ss >> colVal;
	cellDeltas[1]    = std::stod(colVal);
	cellDeltasInv[1] = one/cellDeltas[1];
      }

      else if (colVal == "dz"){
	ss >> colVal;
	cellDeltas[2]    = std::stod(colVal);
	cellDeltasInv[2] = one/cellDeltas[2];
      }

      else if (colVal == "sampleMeshSize"){
	ss >> colVal;
	numGptSampleMesh = std::stoi(colVal);
      }

      else if (colVal == "stencilMeshSize"){
	ss >> colVal;
	numGptStencilMesh = std::stoi(colVal);
      }

      else if (colVal == "stencilSize"){
	ss >> colVal;
	stencilSize = std::stoi(colVal);
      }

    }
  source.close();
}

}}
#endif
