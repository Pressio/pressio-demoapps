
import pathlib, sys, numpy as np
import matplotlib.pyplot as plt
import pressiodemoapps as pda
file_path = pathlib.Path(__file__).parent.absolute()

if __name__ == '__main__':
  meshPath = str(file_path) + "/mesh"
  meshObj  = pda.load_cellcentered_uniform_mesh(meshPath)
  problem  = pda.Euler1d.Sod
  scheme   = pda.InviscidFluxReconstruction.Weno3
  appObj   = pda.create_problem(meshObj, problem, scheme)

  yn = appObj.initialCondition()
  dt = 0.001
  Nsteps = 200
  pda.advanceRK4(appObj, yn, dt, Nsteps)

  x = meshObj.viewX()
  # plot only density
  plt.plot(x, yn[0:-1:3])
  plt.show()
