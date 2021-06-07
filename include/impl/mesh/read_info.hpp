
#ifndef PRESSIODEMOAPPS_MESH_READ_INFO_HPP_
#define PRESSIODEMOAPPS_MESH_READ_INFO_HPP_

namespace pressiodemoapps{ namespace impl{

template<class sc_t, class int_t>
void readMeshInfo(const std::string & meshDir,
		  int & dim,
		  std::array<sc_t, 4> & cellDeltas,
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
	//std::cout << "dim = " << dim << "\n";
      }

      else if (colVal == "dx"){
	ss >> colVal;
	cellDeltas[0] = std::stod(colVal);
	cellDeltas[1] = one/cellDeltas[0];
	//std::cout << "dx = " << cellDeltas[0] << "\n";
      }

      else if (colVal == "dy"){
	ss >> colVal;
	cellDeltas[2] = std::stod(colVal);
	cellDeltas[3] = one/cellDeltas[2];
	//std::cout << "dy = " << cellDeltas[2] << "\n";
      }

      else if (colVal == "sampleMeshSize"){
	ss >> colVal;
	numGptSampleMesh = std::stoi(colVal);
	//std::cout << "numGptSampleMesh = " << numGptSampleMesh << "\n";
      }

      else if (colVal == "stencilMeshSize"){
	ss >> colVal;
	numGptStencilMesh = std::stoi(colVal);
	//std::cout << "numGptStencilMesh = " << numGptStencilMesh << "\n";
      }

      else if (colVal == "stencilSize"){
	ss >> colVal;
	stencilSize = std::stoi(colVal);
	//std::cout << "stencilSize = " << stencilSize << "\n";
      }
    }
  source.close();
}

}}
#endif
