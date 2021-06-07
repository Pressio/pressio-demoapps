#!/usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from numpy import linalg as LA
import re

gamma_ = (5.+2.)/5.
gammaMinusOne_ = gamma_ - 1.
gammaMinusOneInv_ = 1./(gamma_ - 1.)
gammaMinusOneDiv16_ = gammaMinusOne_/(8. * gamma_ * np.pi * np.pi)

def computePressure(rho, u, v, E):
  vel = u**2 + v**2
  return (gamma_ - 1.) * (E - rho*vel*0.5)

def computeMach(rho, u, v, p):
  return np.sqrt((u**2+v**2))/np.sqrt(gamma_*p/rho)

mycmap = cm.get_cmap('gist_rainbow_r', 20)

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
  fomTotDofs = nx*ny*4
  levels = 10 #np.linspace(0.0, 1.8, 40)

  fomCoords = np.loadtxt('coordinates.dat', dtype=float)
  x_fom, y_fom = fomCoords[:,1], fomCoords[:,2]
  x_fom = np.reshape(x_fom, (ny,nx))
  y_fom = np.reshape(y_fom, (ny,nx))

  fomTestD = np.fromfile("solution.bin")
  nt = int(np.size(fomTestD)/fomTotDofs)
  print("fomTest: nt = ", nt)
  fomTestD = np.reshape(fomTestD, (nt, fomTotDofs))

  # fomS = fomTestD[-1,:]
  # fomS = np.reshape(fomS, (nx*ny, 4))
  # rho = fomS[:,0]
  # u   = fomS[:,1]/rho
  # v   = fomS[:,2]/rho
  # p   = computePressure(rho, u, v, fomS[:,3])
  # np.savetxt("p.txt", p)

  # rho1 = np.reshape(rho, (nx,ny))
  # u1 = np.reshape(u, (nx,ny))
  # v1 = np.reshape(v, (nx,ny))
  # p1 = np.reshape(p, (nx,ny))
  # print("\n")
  # print("rho: ", np.min(rho), np.max(rho))
  # print("u:   ", np.min(u), np.max(u))
  # print("v:   ", np.min(v), np.max(v))
  # print("p:   ", np.min(p), np.max(p))

  fig = plt.figure(1)
  for i in range(0,1):#np.shape(fomTestD)[0]):
    fomS = fomTestD[i,:]
    fomS = np.reshape(fomS, (nx*ny, 4))
    rho = fomS[:,0]
    u   = fomS[:,1]/rho
    v   = fomS[:,2]/rho
    p   = computePressure(rho, u, v, fomS[:,3])

    rho1 = np.reshape(rho, (nx,ny))
    p1 = np.reshape(p, (nx,ny))

    plt.clf()
    ax = plt.gca()
    h = plt.contourf(x_fom, y_fom, p1, cmap=mycmap)
    ax.set_aspect(aspect=1.)
    plt.colorbar()
    ax.set_xlim([-0.5, .5])
    ax.set_ylim([-0.5, .5])
    # plt.pause(0.001)
    plt.show()


  # fig = plt.figure(1)
  # # ax = fig.gca(projection='3d')
  # # h = ax.plot_surface(x_fom,y_fom, p1, cmap=cm.jet)
  # ax = plt.gca()
  # h = plt.contourf(x_fom, y_fom, p1, cmap=mycmap)
  # ax.set_aspect(aspect=1.)
  # plt.colorbar()
  # # ax.set_xlim([-10, 10])
  # # ax.set_ylim([-10, 10])
  # ax.set_xlim([-0.5, .5])
  # ax.set_ylim([-0.5, .5])

  # # ----------------------------------------------------
  # rho,u,v,p = sol(x_fom, y_fom, ptime)
  # # print("rho: ", np.min(rho), np.max(rho))
  # # print("u:   ", np.min(u), np.max(u))
  # # print("v:   ", np.min(v), np.max(v))
  # # print("p:   ", np.min(p), np.max(p))

  # print(LA.norm(rho1-rho))
  # print(LA.norm(rho1-rho, 1))
  # print(LA.norm(rho1-rho, np.inf))

  # fig = plt.figure(2)
  # ax = fig.gca(projection='3d')
  # h = ax.plot_surface(x_fom,y_fom, rho, cmap=cm.jet)
  # # ax = plt.gca()
  # # h2 = plt.contourf(x_fom, y_fom, rho, levels, cmap=mycmap)
  # # plt.colorbar()
  # # ax.set_aspect(aspect=1.)
  # ax.set_xlim([0, 2])
  # ax.set_ylim([0, 2])
  # plt.show()

