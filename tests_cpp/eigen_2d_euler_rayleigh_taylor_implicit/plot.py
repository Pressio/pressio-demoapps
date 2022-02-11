#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from numpy import linalg as LA
import re

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

  fomTestD = np.fromfile("rt2d_solution.bin")
  nt = int(np.size(fomTestD)/fomTotDofs)
  print("fomTest: nt = ", nt)
  fomTestD = np.reshape(fomTestD, (nt, fomTotDofs))

  # cm1 = plt.cm.get_cmap('cividis')
  fig = plt.figure(1)
  for i in range(nt-1, nt):
    fomS = fomTestD[i,:]
    fomS = np.reshape(fomS, (nx*ny, 4))
    rho = fomS[:,0]
    rho1 = np.reshape(rho, (ny,nx))
    print(np.min(rho1), np.max(rho1))
    # p1 = np.reshape(p, (ny,nx))

    plt.clf()
    ax = plt.gca()
    h = plt.contourf(x_fom, y_fom, rho1, 20)
    #, vmin=0.95, vmax=2.55)
    ax.set_aspect(aspect=1.)
    plt.colorbar()
    # plt.pause(0.001)
    # fig.savefig("density.png", format="png",
    #   bbox_inches='tight', dpi=300, transparent=True)
    plt.show()
