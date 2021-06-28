#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from numpy import linalg as LA
import re

gamma_ = (5.+2.)/5.
gammaMinusOne_ = gamma_ - 1.
gammaMinusOneDiv16_ = gammaMinusOne_/(8. * gamma_ * np.pi * np.pi)

def sol(xin, t):
  rho = 1. + 0.2*np.sin(np.pi*(xin-t))
  u = np.ones(len(rho))
  p = np.ones(len(rho))
  return rho,u,p

def computeEnergy(prim):
  usq = prim[0]*prim[0]
  vsq = prim[1]*prim[1]
  return prim[3]*gammaMinusOneInv_ + 0.5*prim[0]*(usq+vsq)

def computePressure(rho, u, E):
  return (gamma_ - 1.) * (E - rho*u*u*0.5)

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
  for nx in [10,20,40,80,160,320]:
    fomTotDofs = nx*3

    D = np.fromfile("1d_euler_convergence_weno3_sol"+str(nx)+".bin")
    nt = int(np.size(D)/fomTotDofs)
    print("fomTest: nt = ", nt)
    D = np.reshape(D, (nt, fomTotDofs))

    D = D[-1, :]
    D2 = np.reshape(D, (nx, 3))
    rho = D2[:,0]
    u   = D2[:,1]/rho
    p   = computePressure(rho, u, D2[:,2])

    coo = np.loadtxt("coords"+str(nx)+".dat", dtype=float)
    x_fom = coo[:,1]

    rho1,u1,p1 = sol(x_fom, 2.)
    errors.append( LA.norm(rho-rho1, np.inf) )

    # fig = plt.figure(1)
    # ax = plt.gca()
    # plt.plot(x_fom, rho, '-')
    # plt.plot(x_fom, rho1, 'o')
    # ax.set_xlim([-1, 1])
    # #ax.set_ylim([-1, 1])
    # plt.show()

  goldErr = [0.15448981799473138, \
             0.07049282248026001, \
             0.027574028896264702, \
             0.009011655383856065, \
             0.0018161054795283738, \
             0.00014573402681006264]
  print(goldErr)
  print(errors)

  assert(np.allclose(goldErr, errors, rtol=1e-10, atol=1e-12))
  for i in range(1,6):
    print( np.log2(errors[i-1]/errors[i]) )
