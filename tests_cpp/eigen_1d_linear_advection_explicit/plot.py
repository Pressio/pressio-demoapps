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
  fomTotDofs = nx

  fomCoords = np.loadtxt('coordinates.dat', dtype=float)
  x_fom = fomCoords[:,1]

  data = np.fromfile("sod1d_solution.bin")
  nt = int(np.size(data)/fomTotDofs)
  print("fomTest: nt = ", nt)
  data = np.reshape(data, (nt, fomTotDofs))

  fig = plt.figure(1)
  y0 = np.reshape(data[0,:], (nx, 1))[:,0]
  y1 = np.reshape(data[int(nt/2),:], (nx, 1))[:,0]
  y2 = np.reshape(data[nt-1,:], (nx, 1))[:,0]
  plt.plot(x_fom, y0, '-r', label='Initial condition')
  plt.plot(x_fom, y1, '-g', label='Half simulation time')
  plt.plot(x_fom[0::8], y2[0::8], 'ob', label='Final time')
  plt.xlabel("x", fontsize=12)
  plt.ylabel("Solution", fontsize=12)
  plt.legend()
  fig.savefig("solution.png", format="png", bbox_inches='tight', dpi=450)
  plt.show()
