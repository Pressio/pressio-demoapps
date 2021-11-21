
import numpy as np
import sys, os, math
from numpy import linalg as LA
import matplotlib.pyplot as plt

gamma = (5.+2.)/5.

def analytical(xin, yin, t):
  return 1. + 0.2*np.sin(np.pi*(xin+yin-t*2.))

def computePressure(rho, u, v, E):
  vel = u**2 + v**2
  return (gamma - 1.) * (E - rho*vel*0.5)

if __name__== "__main__":
  st = int(sys.argv[1])

  nx=25
  ny=25
  dt=0.001
  finalTime = 2.
  fomTotDofs = nx*ny*4

  fomCoords = np.loadtxt('coordinates.dat', dtype=float)
  x_fom, y_fom = fomCoords[:,1], fomCoords[:,2]

  D = np.fromfile("smooth2d_solution.bin")
  nt = int(np.size(D)/fomTotDofs)
  D0 = np.reshape(D, (nt, fomTotDofs))
  print(D0.shape)
  D = D0[-1, :]
  D = np.reshape(D, (nx*ny, 4))
  rho = D[:,0]
  u   = D[:,1]/rho
  v   = D[:,2]/rho
  p   = computePressure(rho, u, v, D[:,3])

  goldRho = analytical(x_fom, y_fom, finalTime)
  err = LA.norm(rho-goldRho, np.inf)
  print("Linf error wrt analytical = ", err)

  if(st==3):
    assert(math.isclose(err, 0.19629187623144784))
  if(st==5):
    assert(math.isclose(err, 0.09177511987157294))
  elif(st==7):
    assert(math.isclose(err, 0.0007232573399452713))
