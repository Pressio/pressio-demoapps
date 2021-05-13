
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

from natural_order_mesh import NatOrdMeshRow
from mesh_utils import *

def main(workDir, debug, Nx, Ny, \
         plotting, orderingType, \
         xL, xR, yL, yR, \
         stencilSize, plotFontSize, enablePeriodicBc):

  L = [xR-xL, yR-yL]
  numCells = Nx*Ny
  dx,dy = L[0]/Nx, L[1]/Ny
  print ("Bounds for x = ", xL, xR, " with dx = ", dx)
  print ("Bounds for y = ", yL, yR, " with dy = ", dy)

  if (plotting != "none"):
    figFM = plt.figure(0)
    axFM = figFM.gca()

  # ---------------------------------------------
  if orderingType in ["naturalRow", "naturalRowRcm"]:
    meshObj = NatOrdMeshRow(Nx, Ny, dx, dy, xL,xR,yL,yR, stencilSize)
  else:
    print("invalid orderingType = ", orderingType)
    sys.exit(1)


  [x, y] = meshObj.getXY()
  gids = meshObj.getGIDs()
  G = meshObj.getGraph()
  spMat = meshObj.getSparseMatrixRepr()

  if debug:
    print("full mesh connectivity")
    printDicPretty(G)
    print("\n")

  # -----------------------------------------------------
  # if needed, apply reverse cuthill mckee (RCM)
  # -----------------------------------------------------
  if (orderingType == "rcm"):
    # the starting matrix will always be symmetric because it is full mesh
    # with a symmetric stencil
    [rcmPermInd, spMatRCM] = reverseCuthillMckee(spMat, symmetric=True)

    # hash table to map permutation indices to old indices
    # key = the old cell index
    # value = the new index of the cell after RCM
    ht = {}
    for i in range(len(rcmPermInd)):
      ht[rcmPermInd[i]] = i

    # use the permuation indices to find the new indexing of the cells
    newFullMeshGraph = {}
    for key, value in G.items():
      newFullMeshGraph[ht[key]] = [ht[i] for i in value]
    # sort by key
    newFullMeshGraph = collections.OrderedDict(sorted(newFullMeshGraph.items()))
    print("full mesh connectivity after RCM")
    #printDicPretty(newFullMeshGraph)
    print("\n")

    # this is ugly but works, replace globals
    G = newFullMeshGraph
    spMat = spMatRCM
    #with the permuation indices, update the gids labeling
    x, y = x[rcmPermInd], y[rcmPermInd]
    # if plotting != "none":
    #  plotLabels(x,y,dx, dy, gids,ax10)
    #  ax11.spy(spMatRCM)
    #  ax10.set_aspect(1.0)
    #  ax10.set_xlim(xL,xR)
    #  ax10.set_ylim(yL,yR)

  # -----------------------------------------------------
  meshSize = len(G)
  print ("meshSize = ", meshSize)

  # -----------------------------------------------------
  # info file
  f = open(workDir + "/info.dat","w+")
  f.write("xMin %.14f\n" % xL)
  f.write("xMax %.14f\n" % xR)
  f.write("yMin %.14f\n" % yL)
  f.write("yMax %.14f\n" % yR)
  f.write("dx %.14f\n" % dx)
  f.write("dy %.14f\n" % dy)
  f.write("sampleMeshSize %8d\n"  % meshSize)
  f.write("stencilMeshSize %8d\n" % meshSize)
  f.write("stencilSize %2d\n" % stencilSize)
  f.write("nx %8d\n" % Nx)
  f.write("ny %8d\n" % Ny)
  f.close()

  # -----------------------------------------------------
  # connectivity file
  f = open(workDir + "/connectivity.dat","w+")
  for k in sorted(G.keys()):
    f.write("%8d " % k)
    for i in G[k]:
      f.write("%8d " % i)
    f.write("\n")
  f.close()

  # -----------------------------------------------------
  # print coordinates file
  f = open(workDir + "/coordinates.dat","w+")
  for k in sorted(G.keys()):
    f.write("%8d " % k)
    f.write("%.14f " % x[k])
    f.write("%.14f " % y[k])
    f.write("\n")
  f.close()

  if plotting != "none":
    plotLabels(x, y, dx, dy, gids, axFM, fontSz=plotFontSize)
    axFM.set_aspect(1.0)
    axFM.set_xlim(xL,xR)
    axFM.set_ylim(yL,yR)

  if plotting == "show":
    plt.show()
  elif plotting =="print":
    plt.savefig(workDir+'/full.pdf', dpi=150, format='pdf')


###############################
if __name__== "__main__":
###############################
  # generates a mesh graph for 2d (x,y) cell-centered grid
  # where x is horizontal, y is vertical

  parser = ArgumentParser()
  parser.add_argument(
    "--outDir", dest="outDir",
    help="Full path to where to store all the mesh output files.")

  parser.add_argument(
    "-n", "--numCells",
    nargs="*", type=int, dest="numCells",
    help="Num of cells along x and y. \
If you only pass one value, I use the same for both axis.")

  parser.add_argument(
    "-b", "--bounds",
    nargs="*", dest="bounds", type=float,
    help="Domain bounds along x and y. \
If you only pass the bounds for x, then I use the same for y.")

  parser.add_argument(
    "-o", "--ordering",
    dest="orderingType",
    default="naturalRow",
    help="What cell ordering you want: use <naturalRow> for \
row natural ordering of the mesh cells, use <naturalRowRcm> for \
applying reverse Cuthill-Mckee")

  parser.add_argument(
    "-p", "--plotting",
    nargs="*",
    dest="plottingInfo",
    default="none",
    help="What type of plotting you want: use <show> for showing plots \
use <print> for printing only, use <none> for no plots")

  parser.add_argument(
    "-s", "--stencilSize",
    type=int, dest="stencilSize", default=3,
    help="Stencil size for connectivity: choices are 3,5,7.")

  parser.add_argument(
    "--periodic",
    type=bool, dest="periodic", default=True,
    help="True/False to enable/disable periodic BC so that mesh connectivity\
has that info embedded in.")

  parser.add_argument(
    "-d", "--debug",
    type=bool, dest="debug", default=False,
    help="True/False to eanble debug mode.")

  args = parser.parse_args()

  # -------------------------------------------------------
  # check that args make sense
  assert(args.periodic == True)
  assert(len(args.numCells) == 1 or len(args.numCells) == 2)
  assert(len(args.bounds) == 4 or len(args.bounds) == 2)
  assert( args.stencilSize in [3, 5, 7] )
  assert( args.orderingType in ["naturalRow", "naturalRowRcm"] )

  nx = int(args.numCells[0])
  ny = nx
  if len(args.numCells) == 2:
    ny = int(args.numCells[1])

  xL = float(args.bounds[0])
  xR = float(args.bounds[1])
  yL,yR = xL, xR
  if len(args.bounds) == 4:
      yL = float(args.bounds[2])
      yR = float(args.bounds[3])

  if (xR <= xL):
    print(" xRight <= xLeft: the right bound must be larger than left bound")
    sys.exit(1)
  if (yR <= yL):
    print(" yRight <= yLeft: the right bound must be larger than left bound")
    sys.exit(1)

  # check if working dir exists, if not, make it
  if not os.path.exists(args.outDir):
    os.system('mkdir -p ' + args.outDir)

  plotFontSize = 0
  if len(args.plottingInfo) == 2:
    plotFontSize = int(args.plottingInfo[1])

  main(args.outDir, args.debug,
       nx, ny,
       args.plottingInfo[0],  args.orderingType,
       float(xL), float(xR), float(yL), float(yR),
       args.stencilSize, plotFontSize,
       args.periodic)
