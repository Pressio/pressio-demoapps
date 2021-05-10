
import matplotlib.pyplot as plt
import sys, os, time
import numpy as np
from numpy import linspace, meshgrid
from matplotlib import cm
import collections
from argparse import ArgumentParser
import random
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

# d is a dictionary with the full mesh graph
# returns the graph in sparse matrix format
def convertGraphDicToSparseMatrix(d):
  # vectorize all entries of the graph so that each entry in the new
  # arrays contains (row_index, col_index, 1) refereing to one entry in the matrix

  # the final +1 is to account for the key itself which represents the diagonal
  row_ind = [k for k, v in d.items() for _ in range(len(v)+1)]
  col_ind = []
  for k, v in d.items():
    col_ind.append(k)
    for i in v: col_ind.append(i)

  #col_ind = [int(i) for ids in d.values() for i in ids]
  val = np.ones(len(row_ind)) # can just put ones, since this is a graph

  # from the graph, create a sparse matrix
  spM = sp.csr_matrix(( val, (row_ind, col_ind)) )
  return spM


# spMat is a sparse matrix containing the graph
def reverseCuthillMckee(spMat, symmetric):
  rc = reverse_cuthill_mckee(spMat, symmetric_mode=symmetric)
  I,J=np.ix_(rc,rc)
  X2 = spMat[I,J]
  return [rc, X2]


def plotLabels(x, y, dx, dy, gids, ax, m='o', col='b', ms='0', fontSz=7):
  ax.plot(x, y, m, markerfacecolor=col, markersize=ms)
  cells = []
  for i in range(0, len(x)):
    if fontSz!=0:
      ax.text(x[i], y[i], str(np.int64(gids[i])),
              verticalalignment='center',
              horizontalalignment='center',
              fontsize=fontSz)
    rect = Rectangle((x[i]-dx*0.5, y[i]-dy*0.5), dx, dy)
    cells.append(rect)

  pc = PatchCollection(cells, facecolor='none', edgecolor='k', linewidths=(0.15))
  ax.add_collection(pc)
  ax.set_aspect(aspect=1)


def printDicPretty(d):
  for key, value in d.items():
    print(str(key), value)


def createRandomListTargetResidualGIDs(nx, ny, numCells, targetPct):
  # how many elements of the full mesh the targetPct corresponds to
  # need to get the floor of this to get an integer
  targetSMSize = targetPct * 1e-2 * numCells
  targetSMSize = np.int64(np.floor(targetSMSize))
  rGIDs = random.sample(range(numCells), targetSMSize)
  return rGIDs
