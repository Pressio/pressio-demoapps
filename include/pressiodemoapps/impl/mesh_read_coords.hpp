
#ifndef PRESSIODEMOAPPS_MESH_READ_COORDS_HPP_
#define PRESSIODEMOAPPS_MESH_READ_COORDS_HPP_

namespace pressiodemoapps{ namespace impl{

template<class x_t, class y_t, class z_t>
void readMeshCoordinates(const std::string & meshDir,
			 int dim,
			 x_t& x,
			 y_t& y,
			 z_t& z)
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

      if (dim<=2){
	z(thisGid) = 0.0;
      }
      else{
	ss >> colVal; z(thisGid) = std::stod(colVal);
      }
    }
  source.close();
}

}}
#endif
