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
  print(nx)
  fomTotDofs = nx*3

  fomCoords = np.loadtxt('coordinates.dat', dtype=float)
  x_fom = fomCoords[:,1]

  fomTestD = np.fromfile("sod1d_solution.bin")
  nt = int(np.size(fomTestD)/fomTotDofs)
  print("fomTest: nt = ", nt)
  fomTestD = np.reshape(fomTestD, (nt, fomTotDofs))

  fig = plt.figure(1)
  rho0 = np.reshape(fomTestD[0,:], (nx, 3))[:,0]
  rho1 = np.reshape(fomTestD[int(nt/3),:], (nx, 3))[:,0]
  rho2 = np.reshape(fomTestD[nt-1,:], (nx, 3))[:,0]
  plt.plot(x_fom, rho0, '-r')
  plt.plot(x_fom, rho1, '-g')
  plt.plot(x_fom, rho2, '-b')
  plt.show()
