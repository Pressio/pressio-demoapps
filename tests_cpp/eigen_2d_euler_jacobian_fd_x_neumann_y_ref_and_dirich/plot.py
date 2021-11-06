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

  coords = np.loadtxt('coordinates.dat', dtype=float)
  x,y = coords[:,1], coords[:,2]
  x = np.reshape(x, (ny,nx))
  y = np.reshape(y, (ny,nx))

  D = np.loadtxt("IC.txt")
  rho = D[0::4]
  xm  = D[1::4]
  ym  = D[2::4]
  en  = D[3::4]

  for i in range(4):
    F = np.reshape(D[i::4], (ny,nx))
    print(np.min(F), np.max(F))

    fig = plt.figure(i)
    ax = plt.gca()
    h = plt.contourf(x, y, F)
    ax.set_aspect(aspect=1.)
    plt.colorbar()
  plt.show()
