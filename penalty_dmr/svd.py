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
  # figure out num of cells
  nx = extractN('nx')
  ny = extractN('ny')
  print(nx, ny)
  # total dofs
  fomTotDofs = nx*ny*4

  fomSnaps = np.fromfile("fom_dmr_solution.bin")
  # number of snapshots
  nt = int(np.size(fomSnaps)/fomTotDofs)
  print("fom: num of snapshots = ", nt)
  fomSnaps = np.reshape(fomSnaps, (nt, fomTotDofs))

  # do svd here
  # ...
