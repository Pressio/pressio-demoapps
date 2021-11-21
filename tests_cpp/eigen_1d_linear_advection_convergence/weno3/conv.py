#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from numpy import linalg as LA
import re

# https://link.springer.com/content/pdf/10.1007/s12591-019-00508-5.pdf

def sol2(xin, t):
  pixmt = np.pi*(xin-t)
  return np.sin(pixmt)**4

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
  meshes = [80, 160, 320, 640]

  errors = []
  for nx in meshes:
    D = np.fromfile("eigen_1d_linear_adv_convergence_weno3_sol"+str(nx)+".bin")
    nt = int(np.size(D)/nx)
    print("fomTest: nt = ", nt)
    D = np.reshape(D, (nt, nx))

    x = np.loadtxt('coords'+str(nx)+'.dat', dtype=float)[:,1]

    anY = sol2(x, 2.)
    errors.append( LA.norm(D[-1,:]-anY, np.inf) )

  goldErr = [0.12288453364448215, \
             0.0497108502484469, \
             0.01649554220169369, \
             0.00374869100239660]
  print(goldErr)

  print(errors)
  for i in range(1,len(meshes)):
    print( np.log2(errors[i-1]/errors[i]) )

  assert(np.allclose(goldErr, errors, rtol=1e-10, atol=1e-12))
