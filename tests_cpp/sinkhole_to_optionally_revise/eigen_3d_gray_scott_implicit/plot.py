#!/usr/bin/env python3

import plotly as py
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from numpy import linalg as LA
import re
from mayavi import mlab

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
  numDofPerCell = 2
  numCells   = nx*ny*nz
  fomTotDofs = numCells*numDofPerCell
  print(nx)

  coords = np.loadtxt('coordinates.dat', dtype=float)
  x, y, z = coords[:,1], coords[:,2], coords[:,3]

  data = np.fromfile("gs_3d_solution.bin")
  nt = int(np.size(data)/fomTotDofs)
  print("fomTest: nt = ", nt)
  data = np.reshape(data, (nt, fomTotDofs))

  X, Y, Z = np.mgrid[-1:1:64j, -1:1:64j, -1:1:64j]
  targetSol= data[-1, :]
  u0   = targetSol[0::numDofPerCell]
  v0   = targetSol[1::numDofPerCell]
  print(np.min(u0), np.max(u0))
  print(np.min(v0), np.max(v0))

  fig = go.Figure(data=go.Volume(
      x=X.flatten(),
      y=Y.flatten(),
      z=Z.flatten(),
      value=v0.flatten(),
      isomin=np.min(v0),
      isomax=np.max(v0),
      # slices_y=dict(show=True, locations=[0]),
      opacity=0.2, # needs to be small to see through all surfaces
      surface_count=40 # needs to be a large number for good volume rendering
      ))
  fig.show()

  # for i in range(0, 1): #nt-1,nt):
  #   print(i)
  #   fig = plt.figure()
  #   ax = plt.axes(projection ="3d")
  #   currState = data[i,:]
  #   currState = np.reshape(currState, (numCells, 2))
  #   A = currState[:,0]
  #   B = currState[:,1]

  #   ax.scatter3D(x, y, z, c=B, marker='o', s=25)
  #   # ax.set_xlim([0., 0.4])
  #   # ax.set_ylim([0., 0.4])
  #   # ax.set_zlim([0., 0.4])
  #   # plt.pause(0.001)
  #   plt.show()


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
