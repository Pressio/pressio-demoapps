
#ifndef PRESSIODEMOAPPS_MESH_READ_CONNECTIVITY_HPP_
#define PRESSIODEMOAPPS_MESH_READ_CONNECTIVITY_HPP_

namespace pressiodemoapps{ namespace impl{

template<class g_t>
void read_mesh_connectivity(const std::string & meshDir,
			    g_t & graph,
			    int numCols)
{
  const auto inFile   = meshDir+"/connectivity.dat";
  std::ifstream foundFile(inFile);
  if(!foundFile){
    std::cout << "file not found " << inFile << "\n";
    exit(EXIT_FAILURE);
  }

  std::ifstream source(inFile, std::ios_base::in);
  std::string line;
  std::size_t count = 0;
  while (std::getline(source, line) )
    {
      std::istringstream ss(line);
      std::string colVal;
      ss >> colVal;
      graph(count, 0) = std::stoi(colVal);

      for (auto i=1; i<=numCols-1; ++i){
	ss >> colVal;
	graph(count,i) = std::stoi(colVal);
      }
      count++;
    }
  source.close();
}

}}
#endif
