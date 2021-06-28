
import numpy as np
import sys, os, re, math
from numpy import linalg as LA

gamma = 1.4
gammaMinusOneInv = 1./(gamma-1.)

def computePressure(rho, u, v, w, E):
  return (gamma - 1.) * (E - rho*(u**2+v**2+w**2)*0.5)

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
  nz = extractN('nz')
  numCells   = nx*ny*nz
  fomTotDofs = numCells*5

  D = np.fromfile("3d_sedov_sym_solution.bin")
  nt = int(np.size(D)/fomTotDofs)
  D = np.reshape(D, (nt, fomTotDofs))
  D = D[-1, :]
  np.savetxt("comp_state.txt", D, fmt="%.14f")

  goldD = np.loadtxt("gold_state.txt")
  assert(np.allclose(D.shape, goldD.shape))
  assert(np.isnan(D).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(D, goldD,rtol=1e-8, atol=1e-10))