
import numpy as np
import sys, os, re

gamma = 5./3.

def computePressure(rho, u, v, E):
  vel = u**2 + v**2
  return (gamma - 1.) * (E - rho*vel*0.5)

def extractN(ns):
  reg = re.compile(r''+ns+'.+')
  file1 = open('info.dat', 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return int(strings.group().split()[1])

if __name__== "__main__":
  nx = extractN('nx')
  ny = extractN('ny')
  fomTotDofs = nx*ny*4

  D = np.fromfile("rt2d_solution.bin")
  nt = int(np.size(D)/fomTotDofs)
  D = np.reshape(D, (nt, fomTotDofs))
  D = D[-1, :]
  D = np.reshape(D, (nx*ny, 4))
  rho = D[:,0]
  u   = D[:,1]/rho
  v   = D[:,2]/rho
  p   = computePressure(rho, u, v, D[:,3])
  # np.savetxt("rho.txt", rho)
  # np.savetxt("p.txt", p)

  goldR = np.loadtxt("rho_gold.txt")
  assert(np.allclose(rho.shape, goldR.shape))
  assert(np.isnan(rho).all() == False)
  assert(np.isnan(goldR).all() == False)
  assert(np.allclose(rho, goldR,rtol=1e-10, atol=1e-12))

  goldP = np.loadtxt("p_gold.txt")
  assert(np.allclose(p.shape, goldP.shape))
  assert(np.isnan(p).all() == False)
  assert(np.isnan(goldP).all() == False)
  assert(np.allclose(p, goldP,rtol=1e-10, atol=1e-12))
