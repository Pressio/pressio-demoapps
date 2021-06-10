
import numpy as np
import sys, os

if __name__== "__main__":
  nx=100
  fomTotDofs = nx*3

  D = np.fromfile("sod1d_solution.bin")
  nt = int(np.size(D)/fomTotDofs)
  D = np.reshape(D, (nt, fomTotDofs))
  D = D[-1, :]
  np.savetxt("gigi.txt", D)
  goldD = np.loadtxt("gold.txt")
  # for i in range(fomTotDofs):
  #   print(np.abs(D[i]-goldD[i]))

  assert(np.allclose(D.shape, goldD.shape))
  assert(np.isnan(D).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(D, goldD,rtol=1e-9, atol=1e-11))
