
import numpy as np
import sys, os

if __name__== "__main__":

  # read samples mesh gids
  smgids = np.loadtxt("sample_mesh_gids.dat", dtype=int)
  print(smgids)

  # read full velo
  fv = np.loadtxt("./full/velo.txt")

  # read sample mesh velo
  sv = np.loadtxt("velo.txt")

  maskedVelo = []
  for i in smgids:
    maskedVelo.append(fv[i*3])
    maskedVelo.append(fv[i*3+1])
    maskedVelo.append(fv[i*3+2])
  maskedVelo = np.array(maskedVelo)

  assert(np.allclose(maskedVelo.shape, sv.shape))
  assert(np.isnan(sv).all() == False)
  assert(np.isnan(fv).all() == False)
  assert(np.allclose(sv, maskedVelo,rtol=1e-8, atol=1e-10))

  # import matplotlib.pyplot as plt
  # plt.plot(maskedVelo, 'o')
  # plt.plot(sv, '*')
  # plt.show()
