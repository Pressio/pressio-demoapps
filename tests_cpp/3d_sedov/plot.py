#!/usr/bin/env python3

import plotly as py
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from numpy import linalg as LA
import re
from mayavi import mlab 

gamma_ = (5.+2.)/5.
gammaMinusOne_ = gamma_ - 1.
gammaMinusOneInv_ = 1./(gamma_ - 1.)
gammaMinusOneDiv16_ = gammaMinusOne_/(8. * gamma_ * np.pi * np.pi)

def computePressure(rho, u, v, w, E):
  vel = u**2 + v**2 + w**2
  return (gamma_ - 1.) * (E - rho*vel*0.5)

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

  # X, Y, Z = np.mgrid[-1:1:30j, -1:1:30j, -1:1:30j]
  # values =    np.sin(np.pi*X) * np.cos(np.pi*Z) * np.sin(np.pi*Y)

  # fig = go.Figure(data=go.Volume(
  #     x=X.flatten(),
  #     y=Y.flatten(),
  #     z=Z.flatten(),
  #     value=values.flatten(),
  #     isomin=-0.1,
  #     isomax=0.8,
  #     opacity=0.1, # needs to be small to see through all surfaces
  #     surface_count=21, # needs to be a large number for good volume rendering
  #     ))
  # fig.show()
  # # sys.exit(1)

  nx = extractN('nx')
  ny = extractN('ny')
  nz = extractN('nz')
  numCells   = nx*ny*nz
  fomTotDofs = numCells*5
  levels = 10 #np.linspace(0.0, 1.8, 40)

  coords = np.loadtxt('coordinates.dat', dtype=float)
  x, y, z = coords[:,1], coords[:,2], coords[:,3] 
  # # x_fom = np.reshape(x_fom, (ny,nx))
  # # y_fom = np.reshape(y_fom, (ny,nx))

  fomTestD = np.fromfile("sedov3dsym_solution.bin")
  nt = int(np.size(fomTestD)/fomTotDofs)
  print("fomTest: nt = ", nt)
  fomTestD = np.reshape(fomTestD, (nt, fomTotDofs))

  # X, Y, Z = np.mgrid[0:0.5:25j, 0:0.5:25j, 0:0.5:25j]
  # targetSol= fomTestD[-1, :] 
  # rho0 = targetSol[0::5]
  # u0   = targetSol[1::5]
  # v0   = targetSol[2::5]
  # w0   = targetSol[3::5]
  # E0   = targetSol[4::5]
  # p0   = computePressure(rho0, u0, v0, w0, E0)
  # fig = go.Figure(data=go.Volume(
  #     x=X.flatten(),
  #     y=Y.flatten(),
  #     z=Z.flatten(),
  #     value=p0.flatten(),
  #     isomin=np.min(p0),
  #     isomax=np.max(p0),
  #     # slices_y=dict(show=True, locations=[0]),
  #     opacity=0.2, # needs to be small to see through all surfaces
  #     surface_count=41 # needs to be a large number for good volume rendering
  #     ))
  # fig.show()


  for i in range(nt-1,nt):#nt,50):
    fig = plt.figure()
    ax = plt.axes(projection ="3d")
    fomS = fomTestD[i,:]
    fomS = np.reshape(fomS, (numCells, 5))
    rho = fomS[:,0]
    u   = fomS[:,1]/rho
    v   = fomS[:,2]/rho
    w   = fomS[:,3]/rho
    p   = computePressure(rho, u, v, w, fomS[:,4])
    print(np.min(p), np.max(p))
    # p1 = np.reshape(p, (nx,ny))

    # plt.clf()
    ax.scatter3D(x, y, z, c=p, marker='o', s=5)
    # ax = plt.gca()
    # h = plt.contourf(x_fom, y_fom, p1)
    # ax.set_aspect(aspect=1.)
    # plt.colorbar()
    # ax.set_xlim([0., 1.])
    # ax.set_ylim([0., 0.2])
    # ax.set_zlim([0., 0.2])
    # plt.pause(0.001) 
    # plt.close()
    plt.show()
