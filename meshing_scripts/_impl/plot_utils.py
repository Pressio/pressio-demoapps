
import matplotlib.pyplot as plt
import sys, os, time
import numpy as np
from matplotlib import cm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

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
def plotCells1d(x, dx, gids, ax, darkMode, \
                facecol='w', fontSz=7, alpha=1.):
  cells = []
  for i in range(0, len(x)):
    if fontSz!=0:
      textcolor='w' if darkMode else 'k'
      ax.text(x[i], 0, str(np.int64(gids[i])), \
              verticalalignment='center',
              horizontalalignment='center', \
              fontsize=fontSz, color=textcolor)
    rect = Rectangle((x[i]-dx*0.5, -0.05), dx, 0.1)
    cells.append(rect)

  edgecol = 'k' if darkMode==0 else 'grey'
  pc = PatchCollection(cells, facecolor=facecol, \
                       edgecolor=edgecol, \
                       linewidths=0.5, alpha=alpha)
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
  pc = PatchCollection(cells, facecolor=facecol, \
                       edgecolor=edgecol, linewidths=0.5, alpha=alpha)
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
def plotCells3d(x, y, z, dx, dy, dz, ax, \
                darkMode, facecol, fontSz=7, alpha=0.1):

  edgecol = 'k' if darkMode==0 else 'grey'
  for i in range(0, len(x)):
    cellOrigin = [(x[i]-dx*0.5, y[i]-dy*0.5, z[i]-dz*0.5)]
    cellSize = [(dx, dy, dz)]
    pc = plotCubeAt2(cellOrigin, cellSize, edgecolor='k', alpha=alpha)
    pc.set_facecolor(facecol)
    ax.add_collection3d(pc)

  if darkMode==1: pretty_figure(ax)

# -----------------------------------------------------------
def plotReducedMesh(axSM, dim, plotFontSize, \
                    x, y, z, dx, dy, dz, \
                    gids_sm, fm_to_sm_map, \
                    stencilMeshGIDs, sampleMeshGIDs, darkMode):

  xStencilMesh    = x[stencilMeshGIDs]
  xSampleMeshOnly = x[sampleMeshGIDs]
  stencilCellsColor  = 'w' if not darkMode else 'none'
  residualCellsColor = 'cyan'

  if dim==1:
    plotCells1d(xStencilMesh, dx, gids_sm, axSM, darkMode,
                facecol=stencilCellsColor, fontSz=plotFontSize)
    plotCells1d(xSampleMeshOnly, dx, gids_sm, axSM, darkMode,
                fontSz=0, facecol=residualCellsColor, alpha=0.25)
    axSM.set_xlim(np.min(x)-dx*0.5, np.max(x)+dx*0.5)
    axSM.set_ylim(-0.02, 0.02)
    axSM.set_yticks([])

  if dim==2:
    yStencilMesh    = y[stencilMeshGIDs]
    ySampleMeshOnly = y[sampleMeshGIDs]
    plotCells2d(xStencilMesh, yStencilMesh, dx, dy, gids_sm, axSM,
                darkMode, facecol=stencilCellsColor, fontSz=plotFontSize)
    plotCells2d(xSampleMeshOnly, ySampleMeshOnly, dx, dy, gids_sm, axSM,
                darkMode, facecol=residualCellsColor, fontSz=0, alpha=0.5)
    axSM.set_aspect(1.0)
    axSM.set_xlim(np.min(x)-dx*0.5, np.max(x)+dx*0.5)
    axSM.set_ylim(np.min(y)-dy*0.5, np.max(y)+dy*0.5)

  if dim==3:
    yStencilMesh    = y[stencilMeshGIDs]
    ySampleMeshOnly = y[sampleMeshGIDs]
    zStencilMesh    = z[stencilMeshGIDs]
    zSampleMeshOnly = z[sampleMeshGIDs]
    plotCells3d(xStencilMesh, yStencilMesh, zStencilMesh, dx, dy, dz, \
                axSM, darkMode, facecol='w', fontSz=plotFontSize)
    plotCells3d(xSampleMeshOnly, ySampleMeshOnly, zSampleMeshOnly, dx, dy, dz, \
                axSM, darkMode, facecol='cyan', fontSz=0, alpha=0.25)
