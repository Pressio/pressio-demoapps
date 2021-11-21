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

  fomTestD = np.fromfile("1d_diffreac_solution.bin")
  nt = int(np.size(fomTestD)/fomTotDofs)
  print("fomTest: nt = ", nt)
  fomTestD = np.reshape(fomTestD, (nt, fomTotDofs))

  fig = plt.figure(1)
  S0 = fomTestD[0,:]
  S1 = fomTestD[int(nt/3),:]
  S2 = fomTestD[nt-1,:]
  print(np.min(S2), np.max(S2))
  plt.plot(x_fom, S0, '-r')
  plt.plot(x_fom, S1, '-g')
  plt.plot(x_fom, S2, '-b')
  plt.show()
