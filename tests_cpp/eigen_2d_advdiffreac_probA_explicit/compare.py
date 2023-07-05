
import numpy as np
import sys, os

if __name__== "__main__":
  nx=18
  ny=18
  fomTotDofs = nx*ny

  D = np.fromfile("advdiffreac_2d_solution.bin")
  nt = int(np.size(D)/fomTotDofs)
  D = np.reshape(D, (nt, fomTotDofs))
  D = D[-1, :]
  D = np.reshape(D, (nx*ny, 1))
  h = D[:,0]
  np.savetxt("field.txt", h)

  goldD = np.loadtxt("gold.txt")
  assert(np.allclose(h.shape, goldD.shape))
  assert(np.isnan(h).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(h, goldD,rtol=1e-10, atol=1e-12))
