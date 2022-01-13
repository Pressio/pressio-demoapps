
import numpy as np
import sys, os
from argparse import ArgumentParser

if __name__== "__main__":

  parser = ArgumentParser()
  parser.add_argument("--numdofspercell",
                      dest="numDofsPerCell", type=int)  
  args = parser.parse_args()

  if args.numDofsPerCell <=0 :
    sys.exit("numDofsPerCell must be a positive integer")

  # read samples mesh gids
  sampleMeshGids  = np.loadtxt("sample_mesh_gids.dat", dtype=int)
  stencilMeshGids = np.loadtxt("stencil_mesh_gids.dat", dtype=int)
  print("sample  mesh size: ", len(sampleMeshGids))
  print("stencil mesh size: ", len(stencilMeshGids))

  # read full velo
  fullVelo = np.loadtxt("./full/velo.txt")

  # read full jacob
  fullJacobian = np.loadtxt("./full/jacobian.txt")
  print("full J shape = ", fullJacobian.shape)

  # read sample mesh velo
  sampleVelo = np.loadtxt("velo.txt")

  # read sample mesh jac
  sampleJacobian = np.loadtxt("jacobian.txt")
  print("sample J shape = ", sampleJacobian.shape)

  goldVelo = []
  goldJacob= []
  if args.numDofsPerCell == 1:
    for i in sampleMeshGids:
      goldVelo.append(fullVelo[i])
      goldJacob.append(fullJacobian[i,stencilMeshGids])

  elif args.numDofsPerCell == 2:
    for i in sampleMeshGids:
      for k in range(2):
        goldVelo.append(fullVelo[i*2+k])
        values = []
        for j in stencilMeshGids:
          values.append( fullJacobian[i*2+k,j*2] )
          values.append( fullJacobian[i*2+k,j*2+1] )
        goldJacob.append(values)

  elif args.numDofsPerCell == 3:
    for i in sampleMeshGids:
      for k in range(3):
        goldVelo.append(fullVelo[i*3+k])
        values = []
        for j in stencilMeshGids:
          values.append( fullJacobian[i*3+k,j*3] )
          values.append( fullJacobian[i*3+k,j*3+1] )
          values.append( fullJacobian[i*3+k,j*3+2] )
        goldJacob.append(values)

  elif args.numDofsPerCell == 4:
    for i in sampleMeshGids:
      for k in range(4):
        goldVelo.append(fullVelo[i*4+k])
        values = []
        for j in stencilMeshGids:
          values.append( fullJacobian[i*4+k,j*4] )
          values.append( fullJacobian[i*4+k,j*4+1] )
          values.append( fullJacobian[i*4+k,j*4+2] )
          values.append( fullJacobian[i*4+k,j*4+3] )
        goldJacob.append(values)

  elif args.numDofsPerCell == 5:
    for i in sampleMeshGids:
      for k in range(5):
        goldVelo.append(fullVelo[i*5+k])
        values = []
        for j in stencilMeshGids:
          values.append( fullJacobian[i*5+k,j*5] )
          values.append( fullJacobian[i*5+k,j*5+1] )
          values.append( fullJacobian[i*5+k,j*5+2] )
          values.append( fullJacobian[i*5+k,j*5+3] )
          values.append( fullJacobian[i*5+k,j*5+4] )
        goldJacob.append(values)

  goldVelo  = np.array(goldVelo)
  goldJacob = np.array(goldJacob)
  print("masked J shape = ", goldJacob.shape)

  assert(np.allclose(goldVelo.shape, sampleVelo.shape))
  assert(np.isnan(sampleVelo).all() == False)
  assert(np.isnan(fullVelo).all() == False)
  assert(np.allclose(sampleVelo, goldVelo,rtol=1e-8, atol=1e-10))

  assert(np.allclose(goldJacob.shape, sampleJacobian.shape))
  assert(np.isnan(sampleJacobian).all() == False)
  assert(np.isnan(fullJacobian).all() == False)
  assert(np.allclose(sampleJacobian, goldJacob, rtol=1e-8, atol=1e-10))
