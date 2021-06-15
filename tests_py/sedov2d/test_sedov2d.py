
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
import time
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import cm
import pressiodemoapps as pda

def computePressure(rho, u, v, E):
  gamma_ = (5.+2.)/5.
  vel = u**2 + v**2
  return (gamma_ - 1.) * (E - rho*vel*0.5)

def makePlot(meshPath, yn):
  nx=50
  ny=50
  fomCoords = np.loadtxt(meshPath+"/coordinates.dat", dtype=float)
  x_fom, y_fom = fomCoords[:,1], fomCoords[:,2]
  x_fom = np.reshape(x_fom, (ny,nx))
  y_fom = np.reshape(y_fom, (ny,nx))
  fig = plt.figure(1)
  fomS = np.reshape(yn, (nx*ny, 4))
  rho = fomS[:,0]
  u   = fomS[:,1]/rho
  v   = fomS[:,2]/rho
  p   = computePressure(rho, u, v, fomS[:,3])
  rho1 = np.reshape(rho, (nx,ny))
  p1 = np.reshape(p, (nx,ny))
  mycmap = cm.get_cmap('gist_rainbow_r', 20)
  fig = plt.figure(1)
  ax = plt.gca()
  h = plt.contourf(x_fom, y_fom, p1, cmap=mycmap)
  ax.set_aspect(aspect=1.)
  plt.colorbar()
  ax.set_xlim([-0.5, .5])
  ax.set_ylim([-0.5, .5])
  plt.show()


def test_eval_vel():
  meshPath = str(file_path)
  meshObj  = pda.loadCellCenterUniformMesh(meshPath)
  appObj   = pda.createProblem(meshObj,
                               pda.Euler2d.SedovFull,
                               pda.ReconstructionType.fifthOrderWeno)
  yn = appObj.initialCondition()

  dt = 0.0001
  Nsteps = int(0.05/dt)
  timec = 0.0
  v = appObj.createVelocity()
  y1 = yn.copy()
  y2 = yn.copy()
  start = time.time()
  for step in range(1, Nsteps+1):
    if step % 100 == 0: print("step = ", step)
    appObj.velocity(yn, timec, v)
    y1[:] = yn + dt * v

    appObj.velocity(y1, timec, v)
    y2[:] = (3./4.)*yn + 0.25*y1 + 0.25*dt*v

    appObj.velocity(y2, timec, v)
    yn[:] = (1./3.)*(yn + 2.*y2 + 2.*dt*v)

    timec = step * dt

  end = time.time()
  print(end - start)

  # makePlot(meshPath, yn)

  nx=50
  ny=50
  fomS = np.reshape(yn, (nx*ny, 4))
  rho = fomS[:,0]
  u   = fomS[:,1]/rho
  v   = fomS[:,2]/rho
  p   = computePressure(rho, u, v, fomS[:,3])
  gold = np.loadtxt(str(file_path)+"/p_gold.txt")
  error = LA.norm(p-gold)
  print(error)
  assert( error < 1e-13 )

# ---------------------------
if __name__ == '__main__':
# ---------------------------
  test_eval_vel()
