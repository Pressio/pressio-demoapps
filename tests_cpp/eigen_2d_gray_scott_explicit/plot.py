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
  fomTotDofs = nx*ny*2

  fomCoords = np.loadtxt('coordinates.dat', dtype=float)
  x_fom, y_fom = fomCoords[:,1], fomCoords[:,2]
  x_fom = np.reshape(x_fom, (ny,nx))
  y_fom = np.reshape(y_fom, (ny,nx))

  fomTestD = np.fromfile("gs_2d_solution.bin")
  nt = int(np.size(fomTestD)/fomTotDofs)
  print("fomTest: nt = ", nt)
  fomTestD = np.reshape(fomTestD, (nt, fomTotDofs))

  fig = plt.figure(1)
  for i in range(nt-1, nt):#0, nt):
    fomS = fomTestD[i,:]
    fomS = np.reshape(fomS, (nx*ny, 2))
    A = np.reshape(fomS[:,0], (ny,nx))
    B = np.reshape(fomS[:,1], (ny,nx))
    print(np.min(A), np.max(A))
    print(np.min(B), np.max(B))

    plt.clf()
    plt.subplot(121)
    plt.contourf(x_fom, y_fom, A, 25)
    ax = plt.gca()
    ax.set_aspect(aspect=1.)
    ax.set_xlabel("x", fontsize=12)
    ax.set_ylabel("y", fontsize=12)
    ax.set_title("A", fontsize=15)

    plt.subplot(122)
    plt.contourf(x_fom, y_fom, B, 25)
    ax = plt.gca()
    ax.set_aspect(aspect=1.)
    ax.set_xlabel("x", fontsize=12)
    ax.set_title("B", fontsize=15)

    # plt.colorbar()
    # plt.pause(0.001)
    plt.show()

  #fig.savefig("solution.png", format="png", bbox_inches='tight', dpi=450)
