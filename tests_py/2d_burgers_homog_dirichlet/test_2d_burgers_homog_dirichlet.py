
import pathlib, sys, math
file_path = pathlib.Path(__file__).parent.absolute()

import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
import pressiodemoapps as pda

def test_default_coefficients():
  meshPath = str(file_path)+"/mesh"
  meshO    = pda.load_cellcentered_uniform_mesh(meshPath)
  probId   = pda.AdvectionDiffusion2d.BurgersDirichlet
  invScheme  = pda.InviscidFluxReconstruction.Weno5
  viscScheme = pda.ViscousFluxReconstruction.FirstOrder
  appObj   = pda.create_problem(meshO, probId, invScheme, viscScheme)

  yn = appObj.initialCondition()
  dt = 0.005
  Nsteps = 1200
  pda.advanceRK4(appObj, yn, dt, Nsteps)
  goldD = np.loadtxt(str(file_path)+"/gold.txt")
  assert(np.allclose(yn.shape, goldD.shape))
  assert(np.isnan(yn).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(yn, goldD,rtol=1e-9, atol=1e-11))

#----------------------------
if __name__ == '__main__':
#----------------------------
  test_default_coefficients()
