
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
import time
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import cm
import pressiodemoapps as pda

nx=80
ny=40

def makePlot(meshPath, yn):
  fomCoords = np.loadtxt(meshPath+"/coordinates.dat", dtype=float)
  x_fom, y_fom = fomCoords[:,1], fomCoords[:,2]
  x_fom = np.reshape(x_fom, (ny,nx))
  y_fom = np.reshape(y_fom, (ny,nx))
  fig = plt.figure(1)
  fomS = np.reshape(yn, (nx*ny, 4))
  rho = fomS[:,0]
  rho1 = np.reshape(rho, (ny,nx))

  fig = plt.figure(1)
  ax = plt.gca()
  h = plt.contourf(x_fom, y_fom, rho1)
  ax.set_aspect(aspect=1.)
  plt.colorbar()
  plt.show()

def test_run():
  meshPath = str(file_path)
  meshObj  = pda.loadCellCenterUniformMesh(meshPath)
  appObj   = pda.createProblem(meshObj,
                               pda.Euler2d.NormalShock,
                               pda.InviscidFluxReconstruction.Weno5)
  yn = appObj.initialCondition()
  dt = 0.0005
  Nsteps = int(0.1/dt)
  pda.advanceSSP3(appObj, yn, dt, Nsteps, showProgress=True)
  # makePlot(meshPath, yn)

  fomS = np.reshape(yn, (nx*ny, 4))
  rho = fomS[:,0]
  goldD = np.loadtxt(str(file_path)+"/rho_gold.txt")
  assert(np.allclose(rho.shape, goldD.shape))
  assert(np.isnan(rho).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(rho, goldD,rtol=1e-10, atol=1e-12))

# ---------------------------
if __name__ == '__main__':
# ---------------------------
  test_run()
