#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from numpy import linalg as LA
import re

gamma_ = (5.+2.)/5.
gammaMinusOne_ = gamma_ - 1.
gammaMinusOneInv_ = 1./(gamma_ - 1.)
gammaMinusOneDiv16_ = gammaMinusOne_/(8. * gamma_ * np.pi * np.pi)

def computePressure(rho, u, v, E):
  vel = u**2 + v**2
  return (gamma_ - 1.) * (E - rho*vel*0.5)

def extractN(ns):
  reg = re.compile(r''+ns+'.+')
  file1 = open('info.dat', 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return int(strings.group().split()[1])

##########################
if __name__== "__main__":
##########################
  nx = extractN('nx')
  ny = extractN('ny')
  print(nx, ny)
  fomTotDofs = nx*ny*4

  fomCoords = np.loadtxt('coordinates.dat', dtype=float)
  x_fom, y_fom = fomCoords[:,1], fomCoords[:,2]
  x_fom = np.reshape(x_fom, (ny,nx))
  y_fom = np.reshape(y_fom, (ny,nx))

  fomTestD = np.fromfile("riemann2d_solution.bin")
  nt = int(np.size(fomTestD)/fomTotDofs)
  print("fomTest: nt = ", nt)
  fomTestD = np.reshape(fomTestD, (nt, fomTotDofs))

  cm1 = plt.cm.get_cmap('cividis')
  fig = plt.figure(1)
  for i in range(nt-1, nt):
    fomS = fomTestD[i,:]
    fomS = np.reshape(fomS, (nx*ny, 4))
    rho = fomS[:,0]
    u   = fomS[:,1]/rho
    v   = fomS[:,2]/rho
    p   = computePressure(rho, u, v, fomS[:,3])
    rho1 = np.reshape(rho, (nx,ny))
    print(np.min(rho1), np.max(rho1))
    p1 = np.reshape(p, (nx,ny))

    plt.clf()
    ax = plt.gca()
    h = plt.contourf(x_fom, y_fom, rho1, 18)
    ax.set_aspect(aspect=1.)
    plt.colorbar()
    # ax.set_xlim([0., 0.5])
    # ax.set_ylim([0., 0.5])
    # plt.pause(0.001)
    fig.savefig("density.png", format="png",
      bbox_inches='tight', dpi=300, transparent=True)
    plt.show()
