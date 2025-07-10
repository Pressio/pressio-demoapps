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
  fomTotDofs = nx*ny*3

  fomCoords = np.loadtxt('coordinates.dat', dtype=float)
  x_fom, y_fom = fomCoords[:,1], fomCoords[:,2]
  r = np.sqrt(x_fom**2 + y_fom**2)
  zeroIndx = np.argmin(np.abs(r))
  x_fom = np.reshape(x_fom, (ny,nx))
  y_fom = np.reshape(y_fom, (ny,nx))

  fomTestD = np.fromfile("swe_slipWall2d_solution.bin")
  nt = int(np.size(fomTestD)/fomTotDofs)
  print("fomTest: nt = ", nt)
  fomTestD = np.reshape(fomTestD, (nt, fomTotDofs))
  t = np.linspace(0,7.5,nt)
  fig = plt.figure(1)
  heightAtCenter = np.zeros(nt)
  for i in range(0, nt):#nt-1, nt):
    fomS = fomTestD[i,:]
    fomS = np.reshape(fomS, (nx*ny, 3))
    heightAtCenter[i] = fomS[:,0].flatten()[zeroIndx]

  plt.clf()
  ax = plt.gca()
  hf = plt.plot(t,heightAtCenter)
  plt.xlabel(r'$t$')
  plt.ylabel(r'$h(0,0,t)$')
  plt.show()
