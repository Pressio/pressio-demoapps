
import pathlib, sys, math
file_path = pathlib.Path(__file__).parent.absolute()

import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
import pressiodemoapps as pda


def test_default_coefficients():
  meshPath = str(file_path)+"/mesh"
  meshO    = pda.load_cellcentered_uniform_mesh(meshPath)
  probId   = pda.AdvectionDiffusionReaction2d.ProblemA
  inviscidScheme = pda.InviscidFluxReconstruction.Weno3
  viscousScheme = pda.ViscousFluxReconstruction.FirstOrder
  appObj   = pda.create_problem(meshO, probId, inviscidScheme, viscousScheme)

  yn = appObj.initialCondition()
  dt = 0.005
  Nsteps = 1000
  pda.advanceSSP3(appObj, yn, dt, Nsteps)
  goldD = np.loadtxt(str(file_path)+"/gold1.txt")
  assert(np.allclose(yn.shape, goldD.shape))
  assert(np.isnan(yn).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(yn, goldD,rtol=1e-9, atol=1e-11))

def test_custom_coefficients():
  meshPath = str(file_path)+"/mesh"
  meshO    = pda.load_cellcentered_uniform_mesh(meshPath)
  probId   = pda.AdvectionDiffusionReaction2d.ProblemA
  inviscidScheme = pda.InviscidFluxReconstruction.Weno3
  viscousScheme = pda.ViscousFluxReconstruction.FirstOrder
  ux = 0.5*np.cos(np.pi/3)
  uy = 0.5*np.sin(np.pi/3)
  diff = 0.0005
  sigma = 2.0
  appObj = pda.create_adv_diff_reac_2d_problem_A(meshO, probId, inviscidScheme, viscousScheme, ux, uy, diff, sigma)

  yn = appObj.initialCondition()
  dt = 0.005
  Nsteps = 1000
  pda.advanceSSP3(appObj, yn, dt, Nsteps)
  goldD = np.loadtxt(str(file_path)+"/gold2.txt")
  assert(np.allclose(yn.shape, goldD.shape))
  assert(np.isnan(yn).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(yn, goldD,rtol=1e-9, atol=1e-11))

#----------------------------
if __name__ == '__main__':
#----------------------------
  test_default_coefficients()
  test_custom_coefficients()
