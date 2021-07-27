
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
import time
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import cm
import pressiodemoapps as pda

if __name__ == '__main__':
  dataPath = str(file_path)
  coords = np.loadtxt(dataPath+"/coordinates.dat", dtype=float)
  x, y = coords[:,1], coords[:,2]

  # load data
  snaps = np.loadtxt(dataPath+"/state_snapshots.txt")

  # plot
  fig = plt.figure(1)
  ax = plt.gca()
  mycmap = cm.get_cmap('gist_rainbow_r', 20)
  plt.scatter(x,y, c=snaps[-1,0::4], s=10, marker='s', cmap=mycmap)
  ax.set_aspect(aspect=1.)
  plt.show()
