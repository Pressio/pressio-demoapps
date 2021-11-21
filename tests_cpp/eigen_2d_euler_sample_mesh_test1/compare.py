
import numpy as np
import sys, os

if __name__== "__main__":

  # read samples mesh gids
  sampleMeshGids  = np.loadtxt("sample_mesh_gids.dat", dtype=int)
  stencilMeshGids = np.loadtxt("stencil_mesh_gids.dat", dtype=int)
  print("sample  mesh size: ", len(sampleMeshGids))
  print("stencil mesh size: ", len(stencilMeshGids))

  # read full velo
  fv = np.loadtxt("./full/velo.txt")

  # read full velo
  fullJ = np.loadtxt("./full/jacobian.txt")
  print("full J shape = ", fullJ.shape)

  # read sample mesh velo
  sv = np.loadtxt("velo.txt")

  # read sample mesh jac
  sjac = np.loadtxt("jacobian.txt")
  print("sample J shape = ", sjac.shape)

  maskedVelo = []
  maskedJacob= []
  for i in sampleMeshGids:
    for k in range(4):
      maskedVelo.append(fv[i*4+k])

      values = []
      for j in stencilMeshGids:
        values.append( fullJ[i*4+k,j*4] )
        values.append( fullJ[i*4+k,j*4+1] )
        values.append( fullJ[i*4+k,j*4+2] )
        values.append( fullJ[i*4+k,j*4+3] )
      maskedJacob.append(values)

  maskedVelo  = np.array(maskedVelo)
  maskedJacob = np.array(maskedJacob)
  print("masked J shape = ", maskedJacob.shape)

  assert(np.allclose(maskedVelo.shape, sv.shape))
  assert(np.isnan(sv).all() == False)
  assert(np.isnan(fv).all() == False)
  assert(np.allclose(sv, maskedVelo,rtol=1e-8, atol=1e-10))

  assert(np.allclose(maskedJacob.shape, sjac.shape))
  assert(np.isnan(sjac).all() == False)
  assert(np.isnan(fullJ).all() == False)
  assert(np.allclose(sjac, maskedJacob, rtol=1e-8, atol=1e-10))
