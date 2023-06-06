/*
//@HEADER
// ************************************************************************
//
// mesh_read_info.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
		    int_t & numGptSampleMesh,
                    std::array<int_t, 3> & meshDims)
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

      else if (colVal == "nx"){
	ss >> colVal;
	meshDims[0] = std::stoi(colVal);
      }

      else if (colVal == "ny"){
	ss >> colVal;
	meshDims[1] = std::stoi(colVal);
      }

      else if (colVal == "nz"){
	ss >> colVal;
	meshDims[2] = std::stoi(colVal);
      }

    }
  source.close();
}

}}
#endif
