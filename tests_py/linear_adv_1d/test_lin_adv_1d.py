
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
from numpy import linalg as LA
from pressiodemoapps.enums import *
from pressiodemoapps.mesh import *
from pressiodemoapps.advection import *

def test_create_vel():
  meshPath = str(file_path)
  meshO    = loadCellCenterUniformMesh(meshPath)
  appObj   = createPeriodicLinearAdvection1d(meshO, reconstructWith.fifthOrderWeno)
  v = appObj.createVelocity()
  print(v.shape)

def analytical(x, t):
  return np.sin(np.pi * (x-t) )

def test_eval_vel():
  # create mesh obj
  meshPath = str(file_path)
  meshObj  = CellCenteredUniformMesh(meshPath)
  x = meshObj.viewX()

  # create app object
  appObj   = createPeriodicLinearAdvection1d(meshObj, reconstructWith.fifthOrderWeno)

  y0 = analytical(x, 0.)
  yn = y0.copy()
  v = appObj.createVelocity()

  dt = 0.001
  Nsteps = int(2./dt)
  time = 0.0
  for step in range(1, Nsteps+1):
    if step % 100 == 0:
      print("step = ", step)

    appObj.velocity(yn, time, v)
    k1 = dt * v

    y1 = yn+0.5*k1
    appObj.velocity(y1, time+0.5*dt, v)
    k2 = dt * v

    y2 = yn+0.5*k2
    appObj.velocity(y2, time+0.5*dt, v)
    k3 = dt * v

    y3 = yn+k3
    appObj.velocity(y3, time, v)
    k4 = dt * v

    yn = yn + (k1+2.*k2+2.*k3+k4)/6.
    time = step * dt

  anY = analytical(x, time)
  error = LA.norm(yn-anY, np.inf)
  print(error)
  assert( np.abs(error - 2.828379264019354e-06) < 1e-13 )
  print(yn)

if __name__ == '__main__':
  test_create_vel()
  test_eval_vel()
