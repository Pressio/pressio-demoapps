
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
import time
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import cm
import pressiodemoapps as pda

def computePressure(rho, u, v, w, E):
  gamma_ = (5.+2.)/5.
  vel = u**2 + v**2 + w**2
  return (gamma_ - 1.) * (E - rho*vel*0.5)

def test_run():
  meshPath = str(file_path)
  meshObj  = pda.load_cellcentered_uniform_mesh(meshPath)
  appObj   = pda.create_problem(meshObj,
                                pda.Euler3d.SedovSymmetry,
                                pda.InviscidFluxReconstruction.Weno3)

  dt = 0.0001
  Nsteps = int(0.025/dt)
  yn = appObj.initialCondition()
  pda.advanceSSP3(appObj, yn, dt, Nsteps, showProgress=True)

  gold = np.loadtxt(str(file_path)+"/gold_state.txt")
  assert(np.allclose(yn.shape, gold.shape))
  assert(np.isnan(yn).all() == False)
  assert(np.isnan(gold).all() == False)
  assert(np.allclose(yn, gold,rtol=1e-8, atol=1e-10))

  error = LA.norm(yn-gold)
  print("error = ", error)

  # fig = plt.figure()
  # ax = plt.axes(projection ="3d")
  # coords = np.loadtxt('coordinates.dat', dtype=float)
  # x, y, z = coords[:,1], coords[:,2], coords[:,3]
  # rho = yn[0::5]
  # u   = yn[1::5]/rho
  # v   = yn[2::5]/rho
  # w   = yn[3::5]/rho
  # p   = computePressure(rho, u, v, w, yn[4::5])
  # print(np.min(p), np.max(p))
  # ax.scatter3D(x, y, z, c=p, marker='o', s=5)
  # ax.set_xlim([0., 0.4])
  # ax.set_ylim([0., 0.4])
  # ax.set_zlim([0., 0.4])
  # plt.show()

# ---------------------------
if __name__ == '__main__':
# ---------------------------
  test_run()
