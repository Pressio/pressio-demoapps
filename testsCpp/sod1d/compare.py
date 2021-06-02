
import numpy as np
import sys, os

if __name__== "__main__":
  st = int(sys.argv[1])

  nx=100
  fomTotDofs = nx*3

  D = np.fromfile("solution.bin")
  nt = int(np.size(D)/fomTotDofs)
  D = np.reshape(D, (nt, fomTotDofs))
  D = D[-1, :]
  np.savetxt("gigi.txt", D)

  goldD = np.loadtxt("golds"+str(st)+".txt")
  assert(np.allclose(D.shape, goldD.shape))
  assert(np.isnan(D).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(D, goldD,rtol=1e-8, atol=1e-10))
