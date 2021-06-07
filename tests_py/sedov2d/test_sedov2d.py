
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
import time
import matplotlib.pyplot as plt
from numpy import linalg as LA
from pressiodemoapps.enums import *
from pressiodemoapps.mesh import *
from pressiodemoapps.euler import *

# def computePressure(rho, u, v, E):
#   gamma_ = (5.+2.)/5.
#   vel = u**2 + v**2
#   return (gamma_ - 1.) * (E - rho*vel*0.5)

# def test_create_vel():
#   meshPath = str(file_path)
#   meshO    = loadCellCenterUniformMesh(meshPath)
#   appObj   = createEuler2dProblem(meshO, reconstructWith.fifthOrderWeno, euler2d.sedov, 0)
#   v = appObj.createVelocity()
#   print(v.shape)
#   assert(v.shape[0] == 10000)

def test_eval_vel():
  meshPath = str(file_path)
  meshObj  = loadCellCenterUniformMesh(meshPath)
  appObj = createEuler2dProblem(meshObj, reconstructWith.fifthOrderWeno, euler2d.sedov, 0)
  yn = appObj.initialCondition()
  # yn_add = yn.__array_interface__['data'][0]
  # print("yn: ", yn_add)

  dt = 0.0001
  Nsteps = int(0.05/dt)
  timec = 0.0
  v = appObj.createVelocity()
  # v_add = v.__array_interface__['data'][0]
  # print("v: ", v_add)
  # y1 = yn.copy()
  # y2 = yn.copy()
  start = time.time()
  for step in range(1, Nsteps+1):

    # v = np.sin(yn)
    # for smPt in range(40000):
    #   cellGID = smPt*4
    #   uIndex = cellGID
    #   vIndex = smPt*4
    #   v[vIndex]   = np.sin(yn[uIndex])
    #   v[vIndex+1] = np.sin(yn[uIndex+1])
    #   v[vIndex+2] = np.sin(yn[uIndex+2])
    #   v[vIndex+3] = np.sin(yn[uIndex+3])

    # if step % 200 == 0: print("step = ", step)
    appObj.velocity(yn, timec, v)
    # y1[:] = yn + dt * v

    # appObj.velocity(y1, timec, v)
    # y2[:] = (3./4.)*yn + 0.25*y1 + 0.25*dt*v

    # appObj.velocity(y2, timec, v)
    # yn[:] = (1./3.)*(yn + 2.*y2 + 2.*dt*v)

    # timec = step * dt

  end = time.time()
  print(end - start)

  #nx=50
  #ny=50
  # fomCoords = np.loadtxt(meshPath+"/coordinates.dat", dtype=float)
  # x_fom, y_fom = fomCoords[:,1], fomCoords[:,2]
  # x_fom = np.reshape(x_fom, (ny,nx))
  # y_fom = np.reshape(y_fom, (ny,nx))
  # fig = plt.figure(1)
  #fomS = np.reshape(yn, (nx*ny, 4))
  #rho = fomS[:,0]
  #u   = fomS[:,1]/rho
  #v   = fomS[:,2]/rho
  #p   = computePressure(rho, u, v, fomS[:,3])
  # rho1 = np.reshape(rho, (nx,ny))
  # p1 = np.reshape(p, (nx,ny))
  # from matplotlib import cm
  # mycmap = cm.get_cmap('gist_rainbow_r', 20)
  # fig = plt.figure(1)
  # ax = plt.gca()
  # h = plt.contourf(x_fom, y_fom, p1, cmap=mycmap)
  # ax.set_aspect(aspect=1.)
  # plt.colorbar()
  # ax.set_xlim([-0.5, .5])
  # ax.set_ylim([-0.5, .5])
  # plt.show()

  # nx=50
  # ny=50
  # fomS = np.reshape(yn, (nx*ny, 4))
  # rho = fomS[:,0]
  # u   = fomS[:,1]/rho
  # v   = fomS[:,2]/rho
  # p   = computePressure(rho, u, v, fomS[:,3])
  # gold = np.loadtxt(str(file_path)+"/p_gold.txt")
  # error = LA.norm(p-gold)
  # print(error)
  # assert( error < 1e-13 )

if __name__ == '__main__':
  # test_create_vel()
  test_eval_vel()
