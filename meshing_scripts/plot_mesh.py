
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import sys, os, time, glob, collections
from argparse import ArgumentParser
import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

from _impl.read_mesh_files import *
from _impl.plot_utils import *

#====================================================
def initFigures(dim):
  fig = plt.figure(0)
  if dim<=2:
    ax = fig.gca()
  elif dim==3:
    ax = fig.gca(projection='3d')
  return fig, ax

#====================================================
def printDicPretty(d):
  for key, value in d.items():
    print(str(key), value)

#====================================================
def str2bool(v):
  if isinstance(v, bool):
    return v
  if v.lower() in ('yes', 'true', 't', 'y', '1'):
    return True
  elif v.lower() in ('no', 'false', 'f', 'n', '0'):
    return False
  else:
    raise argparse.ArgumentTypeError('Boolean value expected.')

#====================================================
def main(workDir, debug, plotMode, plotFontSize, darkMode):
  dim, dx,dy,dz, samplesMeshSize, stencilMeshSize, domainBounds, stencilSize = readMeshInfo(workDir)
  stencilMeshIndices, x,y,z  = readMeshCoordinates(workDir, dim)
  G,gids = readMeshConnec(workDir)

  fig, ax = initFigures(dim)
  if dim==1:
    cellsColor  = 'w' if not darkMode else 'none'
    plotCells1d(x, dx, \
                stencilMeshIndices, ax, darkMode, \
                facecol=cellsColor, fontSz=plotFontSize)
    ax.set_xlim(domainBounds[0]-dx*0.5, domainBounds[1]+dx*0.5)
    ax.set_ylim(-0.02, 0.02)
    ax.set_yticks([])

  if dim==2:
    cellsColor  = 'w' if not darkMode else 'none'
    plotCells2d(x, y, dx, dy, \
                stencilMeshIndices, ax, darkMode, \
                facecol=cellsColor, fontSz=plotFontSize)
    ax.set_aspect(1.0)
    ax.set_xlim(np.min(x)-dx*0.5, np.max(x)+dx*0.5)
    ax.set_ylim(np.min(y)-dy*0.5, np.max(y)+dy*0.5)

  if dim==3:
    plotCells3d(x, y, z, dx, dy, dz, ax, darkMode, \
                facecol='w', fontSz=plotFontSize)

  if debug:
    print("natural order full mesh connectivity")
    printDicPretty(G)
    print("\n")

  # gids_sm = list(sm_to_fm_map.keys())
  # # if plotting != "none":
  # #   plotReducedMesh(axSM, dim, plotFontSize, \
  # #                   x, y, z, dx, dy, dz, \
  # #                   gids_sm, fm_to_sm_map,
  # #                   stencilMeshGIDs, sampleMeshGIDs,\
  # #                   darkMode)

  if plotMode == "show":
    plt.show()
  elif plotMode =="print":
    if darkMode==1:
      fig.savefig(workDir+'/mesh.png'  ,dpi=250,format="png", transparent=True)
    else:
      fig.savefig(workDir+'/mesh.png'  ,dpi=250,format="png")


###############################
if __name__== "__main__":
###############################
  parser = ArgumentParser()

  parser.add_argument(
    "--workDir", "--workdir", "--wdir",
    dest="wdir",
    help="Full path to where to mesh files are located.")

  parser.add_argument("-p",
                      nargs="*",
                      dest="plottingInfo",
                      default="none",
                      help="What type of plotting you want:\n"+
                           "use <show> for showing plots,\n"+
                           "use <print> for printing only")

  parser.add_argument(
    "-d", "-debug", "--debug",
    type=bool, dest="debug", default=False)

  args = parser.parse_args()
  # ------------------------------------

  assert(args.wdir != None)

  plotMode = args.plottingInfo[0]
  plotFontSize = 0
  darkMode = 0
  if len(args.plottingInfo) == 2:
    plotFontSize = int(args.plottingInfo[1])
  if len(args.plottingInfo) == 3:
    plotFontSize = int(args.plottingInfo[1])
    darkMode = int(args.plottingInfo[2])

  main(args.wdir, args.debug, plotMode, plotFontSize, darkMode)
