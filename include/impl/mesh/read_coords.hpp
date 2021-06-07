
#ifndef PRESSIODEMOAPPS_MESH_READ_COORDS_HPP_
#define PRESSIODEMOAPPS_MESH_READ_COORDS_HPP_

namespace pressiodemoapps{ namespace impl{

template<class T>
void readMeshCoordinates(const std::string & meshDir,
			 T& x, T& y)
{
  const auto inFile   = meshDir+"/coordinates.dat";
  std::ifstream foundFile(inFile);
  if(!foundFile){
    std::cout << "file not found " << inFile << "\n";
    exit(EXIT_FAILURE);
  }

  std::ifstream source(inFile, std::ios_base::in);
  std::string line;
  while (std::getline(source, line) )
    {
      std::istringstream ss(line);
      std::string colVal;
      ss >> colVal;
      auto thisGid = std::stoi(colVal);
      ss >> colVal; x(thisGid) = std::stod(colVal);
      ss >> colVal; y(thisGid) = std::stod(colVal);
    }
  source.close();
}

}}
#endif
