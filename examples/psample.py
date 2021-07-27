
import scipy.linalg
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
import time, os
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import cm
import pressiodemoapps as pda

def computePressure(rho, u, v, E):
  gamma_ = (5.+2.)/5.
  vel = u**2 + v**2
  return (gamma_ - 1.) * (E - rho*vel*0.5)

def splitPrimitive(S):
  rho = S[0::4, :]
  u = S[1::4, :]
  v = S[2::4, :]
  E = S[3::4, :]
  p = computePressure(rho, u, v, E)
  return [rho, u, v, p]

if __name__ == '__main__':
  dataPath = str(file_path)
  S = np.loadtxt(dataPath+"/state_snapshots.txt")
  # transpose to make each column a sample
  S = S.transpose()
  nCells   = S.shape[0]/4
  nSamples = S.shape[1]
  print("Snapshot shape = ", S.shape[0], S.shape[1])

  # as many mesh points as number of samples
  sampleMeshSize = int(nSamples)
  print("sampleMeshSize = ", sampleMeshSize)

  rho,u,v,pressure = splitPrimitive(S)
  # use density to run things
  U,_,_ = np.linalg.svd(rho, full_matrices=False)
  _,_,P = scipy.linalg.qr(U[:,0:sampleMeshSize].transpose(), pivoting=True )
  sample_mesh = P[0:sampleMeshSize]
  sample_mesh = np.sort(sample_mesh)
  print(sample_mesh)

  # check if subdir exists, if not, make it
  whereToSave = dataPath+"/samplemesh"
  if not os.path.exists(whereToSave):
    os.system('mkdir -p ' + whereToSave)
  np.savetxt(dataPath+'/samplemesh/sample_mesh_gids.dat', sample_mesh, fmt='%8i')
