
import numpy as np
import sys, os

gamma = (5.+2.)/5.

def computePressure(rho, u, v, E):
  vel = u**2 + v**2
  return (gamma - 1.) * (E - rho*vel*0.5)

if __name__== "__main__":
  # st = int(sys.argv[1])

  nx=25
  ny=25
  fomTotDofs = nx*ny*4

  D = np.fromfile("sedov2d_solution.bin")
  nt = int(np.size(D)/fomTotDofs)
  D = np.reshape(D, (nt, fomTotDofs))
  D = D[-1, :]
  D = np.reshape(D, (nx*ny, 4))
  rho = D[:,0]
  u   = D[:,1]/rho
  v   = D[:,2]/rho
  p   = computePressure(rho, u, v, D[:,3])

  np.savetxt("p.txt", p)
  # import matplotlib.pyplot as plt
  # coords = np.loadtxt('coordinates.dat', dtype=float)
  # x_fom, y_fom = coords[:,1], coords[:,2]
  # x_fom = np.reshape(x_fom, (ny,nx))
  # y_fom = np.reshape(y_fom, (ny,nx))
  # rho1 = np.reshape(p, (nx,ny))
  # fig = plt.figure(1)
  # ax = plt.gca()
  # h = plt.contourf(x_fom, y_fom, rho1)
  # ax.set_aspect(aspect=1.)
  # plt.colorbar()
  # # ax.set_xlim([-1., 1.])
  # # ax.set_ylim([-1., 1.])
  # plt.show()

  goldD = np.loadtxt("p_gold.txt")
  assert(np.allclose(p.shape, goldD.shape))
  assert(np.isnan(p).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(p, goldD,rtol=1e-10, atol=1e-12))
