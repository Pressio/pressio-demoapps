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
  nx = extractN('nx')
  ny = extractN('ny')
  nz = extractN('nz')
  numCells   = nx*ny*nz
  fomTotDofs = numCells*5
  print(nx)

  coords = np.loadtxt('coordinates.dat', dtype=float)
  x, y, z = coords[:,1], coords[:,2], coords[:,3] 

  data = np.fromfile("3d_sedov_sym_solution.bin")
  nt = int(np.size(data)/fomTotDofs)
  print("fomTest: nt = ", nt)
  data = np.reshape(data, (nt, fomTotDofs))

  # X, Y, Z = np.mgrid[0:0.5:20j, 0:0.5:20j, 0:0.5:20j]
  # targetSol= data[-1, :] 
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

  for i in range(nt-1,nt):
    print(i)
    fig = plt.figure()
    ax = plt.axes(projection ="3d")
    currState = data[i,:]
    currState = np.reshape(currState, (numCells, 5))
    rho = currState[:,0]
    u   = currState[:,1]/rho
    v   = currState[:,2]/rho
    w   = currState[:,3]/rho
    p   = computePressure(rho, u, v, w, currState[:,4])
    print(np.min(p), np.max(p))

    ax.scatter3D(x, y, z, c=p, marker='o', s=25)
    # ax.set_xlim([0., 0.4])
    # ax.set_ylim([0., 0.4])
    # ax.set_zlim([0., 0.4])
    # plt.pause(0.001) 
    plt.show()


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
