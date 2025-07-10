
import numpy as np
import sys, os

if __name__== "__main__":
  nx=20
  ny=20
  fomTotDofs = nx*ny*2

  D = np.fromfile("burgers2d_solution.bin")
  nt = int(np.size(D)/fomTotDofs)
  D = np.reshape(D, (nt, fomTotDofs))
  D = D[-1, :]
  np.savetxt("field.txt", D)

  goldD = np.loadtxt("gold.txt")
  assert(np.allclose(D.shape, goldD.shape))
  assert(np.isnan(D).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(D, goldD,rtol=1e-10, atol=1e-12))
