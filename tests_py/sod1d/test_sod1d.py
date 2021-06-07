
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
from numpy import linalg as LA
from pressiodemoapps.enums import *
from pressiodemoapps.mesh import *
from pressiodemoapps.euler import *

def test_create_vel():
  meshPath = str(file_path)
  meshO    = loadCellCenterUniformMesh(meshPath)
  appObj   = createEuler1dProblem(meshO, reconstructWith.fifthOrderWeno, euler1d.sod)
  v = appObj.createVelocity()
  print(v.shape)
  assert(v.shape[0] == 300)

def test_eval_vel():
  meshPath = str(file_path)
  meshObj  = loadCellCenterUniformMesh(meshPath)
  x = meshObj.viewX()

  appObj = createEuler1dProblem(meshObj, reconstructWith.fifthOrderWeno, euler1d.sod)

  yn = appObj.initialCondition()
  v = appObj.createVelocity()

  dt = 0.001
  Nsteps = 100
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

  gold = np.loadtxt(str(file_path)+"/gold.txt")
  error = LA.norm(yn-gold)
  print(error)
  assert( error < 1e-13 )
  #import matplotlib.pyplot as plt
  #plt.plot(x, yn[0:-1:3])
  #plt.show()

if __name__ == '__main__':
  test_create_vel()
  test_eval_vel()
