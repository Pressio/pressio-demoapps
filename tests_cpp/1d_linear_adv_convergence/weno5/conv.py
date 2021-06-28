#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from numpy import linalg as LA
import re

def sol(xin, t):
  return np.sin(np.pi*(xin-t))

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
  errors = []
  for nx in [10, 20, 40, 80]:
    D = np.fromfile("1d_linear_adv_convergence_weno5_sol"+str(nx)+".bin")
    nt = int(np.size(D)/nx)
    print("fomTest: nt = ", nt)
    D = np.reshape(D, (nt, nx))

    x = np.loadtxt('coords'+str(nx)+'.dat', dtype=float)[:,1]

    anY = sol(x, 2.)
    errors.append( LA.norm(D[-1,:]-anY, np.inf) )

  goldErr = [0.047441409401830836, \
             0.0025252539521077866, \
             8.82400025065122e-05, \
             2.828379264019354e-06]
  print(goldErr)
  print(errors)
  for i in range(1,4):
    print( np.log2(errors[i-1]/errors[i]) )

  assert(np.allclose(goldErr, errors, rtol=1e-10, atol=1e-12))
