
import matplotlib.pyplot as plt
import sys, os, time
import numpy as np
from matplotlib import cm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

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
  if not isinstance(colors,(list,np.ndarray)):
    colors=["C0"]*len(positions)

  if not isinstance(sizes,(list,np.ndarray)):
    sizes=[(1,1,1)]*len(positions)

  g = []
  for p,s,c in zip(positions,sizes,colors):
    g.append( cuboid_data2(p, size=s) )
  #return Poly3DCollection(np.concatenate(g),
  #                        facecolors=np.repeat(colors,6), **kwargs)
  #return Line3DCollection(np.concatenate(g), colors='k', \
  #  linewidths=0.1, linestyles=':')
  return Poly3DCollection(np.concatenate(g), **kwargs)

#---------------------------------------------------------------------
def pretty_figure(ax, leg=None):
  mycolor = 'w'

  # makes all axes and text whie
  for l in ['bottom', 'left', 'right', 'top']:
    ax.spines[l].set_color(mycolor)
  ax.xaxis.label.set_color(mycolor);
  ax.tick_params(axis='x', colors=mycolor)
  ax.yaxis.label.set_color(mycolor);
  ax.tick_params(axis='y', colors=mycolor)
  if leg!=None:
    for text in leg.get_texts(): text.set_color(mycolor)

#---------------------------------------------------------------------
def remove_ticks(ax, dim):
  if dim >=1:
    ax.set_xticks([])
  if dim>=2:
    ax.set_yticks([])
  if dim==3:
    ax.set_zticks([])

#---------------------------------------------------------------------
def plotCells1d(x, dx, gids, ax, darkMode, \
                facecol='w', fontSz=7, alpha=1.):
  cells = []
  for i in gids:
    ptX = x[i]
    if fontSz!=0:
      textcolor='w' if darkMode else 'k'
      ax.text(ptX, 0, str(np.int64(gids[i])), \
              verticalalignment='center',
              horizontalalignment='center', \
              fontsize=fontSz, color=textcolor)
    rect = Rectangle((ptX-dx*0.5, -0.05), dx, 0.1)
    cells.append(rect)

  edgecol = 'k' if darkMode==0 else 'grey'
  pc = PatchCollection(cells, facecolor=facecol, \
                       edgecolor=edgecol, \
                       linewidths=0.5, alpha=alpha)
  ax.add_collection(pc)
  ax.set_aspect(aspect=1)
  if darkMode==1: 
    pretty_figure(ax)
    remove_ticks(ax,1)

#---------------------------------------------------------------------
def plotCells2d(x, y, dx, dy, gids, ax, darkMode, facecol='w', fontSz=7, alpha=1.):
  cells = []
  for i in gids:
    ptX = x[i]
    ptY = y[i]
    if fontSz!=0:
      textcolor='w' if darkMode else 'k'
      ax.text(ptX, ptY, str(np.int64(i)), verticalalignment='center',
              horizontalalignment='center', fontsize=fontSz, color=textcolor)
    rect = Rectangle((ptX-dx*0.5, ptY-dy*0.5), dx, dy)
    cells.append(rect)

  edgecol = 'k' if darkMode==0 else 'grey'
  pc = PatchCollection(cells, facecolor=facecol, \
                       edgecolor=edgecol, linewidths=0.5, alpha=alpha)
  ax.add_collection(pc)
  ax.set_aspect(aspect=1)
  if darkMode==1: 
    pretty_figure(ax)
    remove_ticks(ax,2)

#---------------------------------------------------------------------
def plotCells3d(x, y, z, dx, dy, dz, gids, ax, \
                darkMode, facecol, fontSz=7, alpha=0.1):

  edgecol = 'k' if darkMode==0 else 'grey'

  for i in gids:
    ptX = x[i]
    ptY = y[i]
    ptZ = z[i]
    cellOrigin = [(ptX-dx*0.5, ptY-dy*0.5, ptZ-dz*0.5)]
    cellSize = [(dx, dy, dz)]
    pc = plotCubeAt2(cellOrigin, cellSize, edgecolor='k', alpha=alpha)
    pc.set_facecolor(facecol)
    ax.add_collection3d(pc)

  if darkMode==1: 
    pretty_figure(ax)
    remove_ticks(ax,3)

