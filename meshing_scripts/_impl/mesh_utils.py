
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
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

# #---------------------------------------------------------------------
# def createRandomListTargetResidualGIDs(nx, ny, numCells, targetPct):
#   # how many elements of the full mesh the targetPct corresponds to
#   # need to get the floor of this to get an integer
#   targetSMSize = targetPct * 1e-2 * numCells
#   targetSMSize = np.int64(np.floor(targetSMSize))
#   rGIDs = random.sample(range(numCells), targetSMSize)
#   return rGIDs

#---------------------------------------------------------------------
def pretty_figure(ax, leg=None):
  mycolor = 'w'

  # # makes all axes and text whie
  # for l in ['bottom', 'left', 'right', 'top']:
  #   ax.spines[l].set_color(mycolor)

  # ax.xaxis.label.set_color(mycolor);
  # ax.tick_params(axis='x', colors=mycolor)
  # ax.yaxis.label.set_color(mycolor);
  # ax.tick_params(axis='y', colors=mycolor)
  # if leg!=None:
  #   for text in leg.get_texts(): text.set_color(mycolor)

#---------------------------------------------------------------------
def convertGraphDicToSparseMatrix(d):
  # vectorize all entries of the graph so that each entry in the new
  # arrays contains (row_index, col_index, 1)
  # refereing to one entry in the matrix

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


#---------------------------------------------------------------------
def reverseCuthillMckee(spMat, symmetric):
  rc = reverse_cuthill_mckee(spMat, symmetric_mode=symmetric)
  I,J=np.ix_(rc,rc)
  X2 = spMat[I,J]
  return [rc, X2]

#---------------------------------------------------------------------
def plotCells1d(x, dx, gids, ax, darkMode, facecol='w', fontSz=7, alpha=1.):
  cells = []
  for i in range(0, len(x)):
    if fontSz!=0:
      textcolor='w' if darkMode else 'k'
      ax.text(x[i], 0, str(np.int64(gids[i])), verticalalignment='center',
              horizontalalignment='center', fontsize=fontSz, color=textcolor)
    rect = Rectangle((x[i]-dx*0.5, -0.05), dx, 0.1)
    cells.append(rect)

  edgecol = 'k' if darkMode==0 else 'grey'
  pc = PatchCollection(cells, facecolor=facecol, edgecolor=edgecol, linewidths=0.5, alpha=alpha)
  ax.add_collection(pc)
  ax.set_aspect(aspect=1)
  if darkMode==1: pretty_figure(ax)

#---------------------------------------------------------------------
def plotCells2d(x, y, dx, dy, gids, ax, darkMode, facecol='w', fontSz=7, alpha=1.):
  cells = []
  for i in range(0, len(x)):
    if fontSz!=0:
      textcolor='w' if darkMode else 'k'
      ax.text(x[i], y[i], str(np.int64(gids[i])), verticalalignment='center',
              horizontalalignment='center', fontsize=fontSz, color=textcolor)
    rect = Rectangle((x[i]-dx*0.5, y[i]-dy*0.5), dx, dy)
    cells.append(rect)

  edgecol = 'k' if darkMode==0 else 'grey'
  pc = PatchCollection(cells, facecolor=facecol, edgecolor=edgecol, linewidths=0.5, alpha=alpha)
  ax.add_collection(pc)
  ax.set_aspect(aspect=1)
  if darkMode==1: pretty_figure(ax)

#---------------------------------------------------------------------
def cuboid_data2(o, size=(1,1,1)):
  X = [[[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]],
       [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]],
       [[1, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1]],
       [[0, 0, 1], [0, 0, 0], [0, 1, 0], [0, 1, 1]],
       [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]],
       [[0, 1, 1], [0, 0, 1], [1, 0, 1], [1, 1, 1]]]
  X = np.array(X).astype(float)
  for i in range(3):
    X[:,:,i] *= size[i]
  X += np.array(o)
  return X

#---------------------------------------------------------------------
def plotCubeAt2(positions, sizes=None, colors=None, **kwargs):
  if not isinstance(colors,(list,np.ndarray)): colors=["C0"]*len(positions)
  if not isinstance(sizes,(list,np.ndarray)): sizes=[(1,1,1)]*len(positions)
  g = []
  for p,s,c in zip(positions,sizes,colors):
    g.append( cuboid_data2(p, size=s) )
  #return Poly3DCollection(np.concatenate(g),
  #                        facecolors=np.repeat(colors,6), **kwargs)
  #return Line3DCollection(np.concatenate(g), colors='k', linewidths=0.1, linestyles=':')
  return Poly3DCollection(np.concatenate(g), **kwargs)

#---------------------------------------------------------------------
def plotCells3d(x, y, z, dx, dy, dz, ax, darkMode, facecol, fontSz=7, alpha=0.1):
  edgecol = 'k' if darkMode==0 else 'grey'
  for i in range(0, len(x)):
    cellOrigin = [(x[i]-dx*0.5, y[i]-dy*0.5, z[i]-dz*0.5)]
    cellSize = [(dx, dy, dz)]
    pc = plotCubeAt2(cellOrigin, cellSize, edgecolor='k', alpha=alpha)
    pc.set_facecolor(facecol)
    ax.add_collection3d(pc)

  if darkMode==1: pretty_figure(ax)
