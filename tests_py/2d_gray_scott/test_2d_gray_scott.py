
import pathlib, sys, math
file_path = pathlib.Path(__file__).parent.absolute()

import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
import pressiodemoapps as pda

def make_plot(state):
  nx, ny = 50, 50
  fomTotDofs = nx*nx*2

  fomCoords = np.loadtxt('./mesh/coordinates.dat', dtype=float)
  x_fom, y_fom = fomCoords[:,1], fomCoords[:,2]
  x_fom = np.reshape(x_fom, (ny,nx))
  y_fom = np.reshape(y_fom, (ny,nx))
  fig = plt.figure(1)
  A = np.reshape(state[0::2], (ny,nx))
  B = np.reshape(state[1::2], (ny,nx))
  print(np.min(A), np.max(A))
  print(np.min(B), np.max(B))
  plt.clf()
  plt.subplot(121)
  plt.contourf(x_fom, y_fom, A, 25)
  ax = plt.gca()
  ax.set_aspect(aspect=1.)
  ax.set_xlabel("x", fontsize=12)
  ax.set_ylabel("y", fontsize=12)
  ax.set_title("A", fontsize=15)

  plt.subplot(122)
  plt.contourf(x_fom, y_fom, B, 25)
  ax = plt.gca()
  ax.set_aspect(aspect=1.)
  ax.set_xlabel("x", fontsize=12)
  ax.set_title("B", fontsize=15)

  # plt.colorbar()
  # plt.pause(0.001)
  plt.show()

probId   = pda.DiffusionReaction2d.GrayScott

def test_default_coefficients():
  meshPath = str(file_path)+"/mesh"
  meshO    = pda.load_cellcentered_uniform_mesh(meshPath)
  scheme = pda.ViscousFluxReconstruction.FirstOrder
  appObj   = pda.create_problem(meshO, probId, scheme)

  yn = appObj.initialCondition()
  dt = 0.8
  Nsteps = int(4000./dt)
  pda.advanceRK4(appObj, yn, dt, Nsteps)
  goldD = np.loadtxt(str(file_path)+"/gold1.txt")
  assert(np.allclose(yn.shape, goldD.shape))
  assert(np.isnan(yn).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(yn, goldD,rtol=1e-9, atol=1e-11))
  # make_plot(yn)

def test_custom_coefficients():
  meshPath = str(file_path)+"/mesh"
  meshO    = pda.load_cellcentered_uniform_mesh(meshPath)
  scheme = pda.ViscousFluxReconstruction.FirstOrder
  appObj   = pda.create_problem(meshO, probId, scheme, 0.0005, 0.0001, 0.045, 0.065)

  yn = appObj.initialCondition()
  dt = 0.8
  Nsteps = int(4000./dt)
  pda.advanceRK4(appObj, yn, dt, Nsteps)
  # np.savetxt("field.txt", yn)
  goldD = np.loadtxt(str(file_path)+"/gold2.txt")
  assert(np.allclose(yn.shape, goldD.shape))
  assert(np.isnan(yn).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(yn, goldD,rtol=1e-9, atol=1e-11))
  # make_plot(yn)

#----------------------------
if __name__ == '__main__':
#----------------------------
  test_default_coefficients()
  test_custom_coefficients()
