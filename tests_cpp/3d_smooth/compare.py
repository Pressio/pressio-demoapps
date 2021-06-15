
import numpy as np
import sys, os, re, math
from numpy import linalg as LA

gamma = 1.4
gammaMinusOneInv = 1./(gamma-1.)

def sol(xin,yin,zin, t):
  return 1. + 0.2*np.sin(np.pi*(xin+yin+zin-t*3.))

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
  st = int(sys.argv[1])
  nx = extractN('nx')
  ny = extractN('ny')
  nz = extractN('nz')
  numCells   = nx*ny*nz
  fomTotDofs = numCells*5

  D = np.fromfile("solution.bin")
  nt = int(np.size(D)/fomTotDofs)
  D = np.reshape(D, (nt, fomTotDofs))
  D = D[-1, :]
  D = np.reshape(D, (numCells, 5))
  rho = D[:,0]
  np.savetxt("rho_comp.txt", rho)

  coo = np.loadtxt("coordinates.dat", dtype=float)
  x_fom = coo[:,1]
  y_fom = coo[:,2]
  z_fom = coo[:,3]

  gold = sol(x_fom, y_fom, z_fom, 2.)
  err = LA.norm(rho-gold, np.inf)
  print(err)

  if(st==3):
    assert(math.isclose(err, 0.19999999613178898))
