
import pathlib, sys, numpy as np
import matplotlib.pyplot as plt
import pressiodemoapps as pda
file_path = pathlib.Path(__file__).parent.absolute()

if __name__ == '__main__':
  meshPath = str(file_path)
  meshObj  = pda.load_cellcentered_uniform_mesh(meshPath)
  problem  = pda.Euler1d.ShuOsher
  scheme   = pda.InviscidFluxReconstruction.Weno5
  appObj   = pda.create_problem(meshObj, problem, scheme)

  yn = appObj.initialCondition()
  dt = 0.001
  Nsteps = int(1.8/dt)
  pda.advanceSSP3(appObj, yn, dt, Nsteps)

  goldD = np.loadtxt(str(file_path)+"/gold.txt")
  assert(np.allclose(yn.shape, goldD.shape))
  assert(np.isnan(yn).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(yn, goldD,rtol=1e-9, atol=1e-11))

  x = meshObj.viewX()
  plt.plot(x, yn[0:-1:3], '-k')
  plt.plot(x, goldD[0:-1:3], 'or')
  plt.show()
