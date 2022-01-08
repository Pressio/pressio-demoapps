
import numpy as np
import sys, os

if __name__== "__main__":
  nx=50
  ny=50
  totDofs = nx*ny

  D = np.fromfile("diffusion_reaction_2d_solution.bin")
  nt = int(np.size(D)/totDofs)
  D = np.reshape(D, (nt, totDofs))
  y = D[-1, :]
  np.savetxt("field.txt", y)

  goldD = np.loadtxt("gold.txt")
  assert(np.allclose(y.shape, goldD.shape))
  assert(np.isnan(y).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(y, goldD,rtol=1e-10, atol=1e-12))
