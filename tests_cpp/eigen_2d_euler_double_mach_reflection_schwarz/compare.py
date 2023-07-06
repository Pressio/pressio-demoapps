
import numpy as np
import sys, os

if __name__== "__main__":
  nx=50
  ny=200
  fomTotDofs = nx*ny*4

  D = np.fromfile("doubleMach2d_solution.bin")
  nt = int(np.size(D)/fomTotDofs)
  D = np.reshape(D, (nt, fomTotDofs))
  D = D[-1, :]
  D = np.reshape(D, (nx*ny, 4))
  rho = D[:,0]
  np.savetxt("rho.txt", rho)

  goldD = np.loadtxt("rho_gold.txt")
  assert(np.allclose(rho.shape, goldD.shape))
  assert(np.isnan(rho).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(rho, goldD,rtol=1e-10, atol=1e-12))
