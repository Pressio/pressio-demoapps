
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
from numpy import linalg as LA
import pressiodemoapps as pda

def test_create_vel():
  meshPath = str(file_path)
  meshO    = pda.loadCellCenterUniformMesh(meshPath)
  appObj   = pda.createProblem(meshO, pda.Euler1d.Sod, pda.ReconstructionType.fifthOrderWeno)
  v = appObj.createVelocity()
  print(v.shape)
  assert(v.shape[0] == 300)

def test_eval_vel():
  meshPath = str(file_path)
  meshObj  = pda.loadCellCenterUniformMesh(meshPath)
  appObj   = pda.createProblem(meshObj, pda.Euler1d.Sod, pda.ReconstructionType.fifthOrderWeno)

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
  #x = meshObj.viewX()
  #plt.plot(x, yn[0:-1:3])
  #plt.show()

if __name__ == '__main__':
  test_create_vel()
  test_eval_vel()
