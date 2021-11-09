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

  coords = np.loadtxt('coordinates.dat', dtype=float)
  x = coords[:,1]

  y = np.loadtxt("IC.txt")
  fig = plt.figure(1)
  plt.plot(x, y, '-r')
  plt.show()
