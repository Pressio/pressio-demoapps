
import numpy as np
import sys, os

gamma = (5.+2.)/5.

def computePressure(rho, u, v, E):
  vel = u**2 + v**2
  return (gamma - 1.) * (E - rho*vel*0.5)

if __name__== "__main__":
  nx=50
  ny=50
  fomTotDofs = nx*ny*4

  D = np.fromfile("sedov2d_solution.bin")
  nt = int(np.size(D)/fomTotDofs)
  D = np.reshape(D, (nt, fomTotDofs))
  D = D[-1, :]
  D = np.reshape(D, (nx*ny, 4))
  rho = D[:,0]
  u   = D[:,1]/rho
  v   = D[:,2]/rho
  p   = computePressure(rho, u, v, D[:,3])
  np.savetxt("gigi.txt", p)

  goldD = np.loadtxt("p_gold.txt")
  assert(np.allclose(p.shape, goldD.shape))
  assert(np.isnan(p).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(p, goldD))#,rtol=1e-10, atol=1e-12))
