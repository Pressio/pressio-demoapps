
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
import time
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import cm
import pressiodemoapps as pda

nx=65
ny=65

def makePlot(meshPath, yn):
  fomCoords = np.loadtxt(meshPath+"/coordinates.dat", dtype=float)
  x_fom, y_fom = fomCoords[:,1], fomCoords[:,2]
  x_fom = np.reshape(x_fom, (ny,nx))
  y_fom = np.reshape(y_fom, (ny,nx))
  fig = plt.figure(1)
  fomS = np.reshape(yn, (nx*ny, 3))
  h = fomS[:,0]
  h1 = np.reshape(h, (ny,nx))

  fig = plt.figure(1)
  ax = plt.gca()
  h = plt.contourf(x_fom, y_fom, h1)
  ax.set_aspect(aspect=1.)
  plt.colorbar()
  plt.show()

def test_run():
  meshPath = str(file_path)
  meshObj  = pda.load_cellcentered_uniform_mesh(meshPath)
  appObj   = pda.create_problem(meshObj,
                               pda.Swe2d.SlipWall,
                               pda.InviscidFluxReconstruction.Weno5)
  yn = appObj.initialCondition()
  dt = 0.005
  Nsteps = int(7.5/dt)
  pda.advanceRK4(appObj, yn, dt, Nsteps, showProgress=True)
  # makePlot(meshPath, yn)

  fomS = np.reshape(yn, (nx*ny, 3))
  h = fomS[:,0]
  # np.savetxt("h.txt", h)
  goldD = np.loadtxt(str(file_path)+"/h_gold.txt")
  assert(np.allclose(h.shape, goldD.shape))
  assert(np.isnan(h).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(h, goldD,rtol=1e-10, atol=1e-12))

# ---------------------------
if __name__ == '__main__':
# ---------------------------
  test_run()
