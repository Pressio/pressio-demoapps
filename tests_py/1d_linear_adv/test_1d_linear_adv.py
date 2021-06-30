
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
from numpy import linalg as LA
import pressiodemoapps as pda

def test_weno3():
  def analytical(x, t):
    pixmt = np.pi * (x-t)
    return np.sin(pixmt)**4

  meshPath = str(file_path)+"/mesh_weno3"
  meshO    = pda.loadCellCenterUniformMesh(meshPath)
  probId   = pda.Advection1d.PeriodicLinear
  appObj   = pda.createProblem(meshO, probId, pda.InviscidFluxReconstruction.Weno3)

  x = meshO.viewX()
  yn = analytical(x, 0.)
  dt = 0.0001
  finalTime = 2.
  Nsteps = int(finalTime/dt)
  pda.advanceRK4(appObj, yn, dt, Nsteps)
  anY = analytical(x, finalTime)
  error = LA.norm(yn-anY, np.inf)
  print(error)
  assert( np.abs(error - 0.00374869100239660) < 1e-13 )
  print(yn)

def test_weno5():
  def analytical(x, t):
    return np.sin(np.pi * (x-t) )

  meshPath = str(file_path)+"/mesh_weno5"
  meshO    = pda.loadCellCenterUniformMesh(meshPath)
  probId   = pda.Advection1d.PeriodicLinear
  appObj   = pda.createProblem(meshO, probId, pda.InviscidFluxReconstruction.Weno5)

  x = meshO.viewX()
  yn = analytical(x, 0.)
  dt = 0.001
  finalTime = 2.
  Nsteps = int(finalTime/dt)
  pda.advanceRK4(appObj, yn, dt, Nsteps)
  anY = analytical(x, finalTime)
  error = LA.norm(yn-anY, np.inf)
  print(error)
  assert( np.abs(error - 2.828379264019354e-06) < 1e-13 )
  print(yn)

#----------------------------
if __name__ == '__main__':
#----------------------------
  test_weno3()
  test_weno5()
