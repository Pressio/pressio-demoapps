
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

from one_dim_mesh import OneDimMesh
from natural_order_mesh import NatOrdMeshRow
from mesh_utils import *
from remove_cells_step import *

# ------------------------------------------------
def checkDomainBoundsCoordinates(xValues, yValues):
  # x coordiantes must be in increasing order
  for i in range(1, len(xValues)):
    assert(xValues[i] > xValues[i-1])

  # y coordiantes must be in increasing order
  for i in range(1, len(yValues)):
    assert(yValues[i] > yValues[i-1])

# ------------------------------------------------
def str2bool(v):
  if isinstance(v, bool):
    return v
  if v.lower() in ('yes', 'true', 't', 'y', '1'):
    return True
  elif v.lower() in ('no', 'false', 'f', 'n', '0'):
    return False
  else:
    raise argparse.ArgumentTypeError('Boolean value expected.')

# ------------------------------------------------
def main(workDir, debug, Nx, Ny, \
         plotting, orderingType, \
         xBd, yBd, stencilSize, \
         plotFontSize, enablePeriodicBc, dim):

  # figure out if domain is a plain rectangle or has a step in it
  plainDomain = True
  print(len(xBd))
  print(xBd)
  if dim!=1 and len(xBd) > 4:
    plainDomain = False

  if plainDomain:
    numCells = Nx*Ny
    L = [xBd[1]-xBd[0], yBd[1]-yBd[0]]
    dx,dy = L[0]/Nx, L[1]/Ny
    print ("Bounds for x = ", xBd[0], xBd[1], " with dx = ", dx)
    print ("Bounds for y = ", yBd[0], yBd[1], " with dy = ", dy)

  if (plotting != "none"):
    figFM = plt.figure(0)
    axFM = figFM.gca()

  # ---------------------------------------------
  x = None
  y = None
  G = None
  gids = None
  if dim == 1:
    meshObj = OneDimMesh(Nx, dx, xBd[0], xBd[1], stencilSize, enablePeriodicBc)
    # get mesh coordinates, gids and graph
    [x, y] = meshObj.getCoordinates()
    gids = meshObj.getGIDs()
    G = meshObj.getGraph()

  elif dim==2 and plainDomain:
    meshObj = NatOrdMeshRow(Nx,Ny, dx,dy,\
                            xBd[0], xBd[1], yBd[0], yBd[1], \
                            stencilSize, enablePeriodicBc)

    # get mesh coordinates, gids and graph
    [x, y] = meshObj.getCoordinates()
    gids = meshObj.getGIDs()
    G = meshObj.getGraph()

  elif dim==2 and not plainDomain:
    minX, maxX = min(xBd), max(xBd)
    minY, maxY = min(yBd), max(yBd)
    L = [maxX-minX, maxY-minY]
    dx,dy = L[0]/Nx, L[1]/Ny

    meshObj = NatOrdMeshRow(Nx, Ny, dx, dy,\
                            minX, maxX, minY, maxY, \
                            stencilSize, False)

    x,y,G,gids = removeStepCells(meshObj, xBd, yBd)

  if debug:
    print("full mesh connectivity")
    printDicPretty(G)
    print("\n")

  # -----------------------------------------------------
  # if needed, apply reverse cuthill mckee (RCM)
  # -----------------------------------------------------
  if (orderingType == "rcm"):
    spMat = meshObj.getSparseMatrixRepr()

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
  if dim==1:
    f.write("dim %1d\n" % 1)
    f.write("xMin %.14f\n" % xBd[0])
    f.write("xMax %.14f\n" % xBd[1])
    f.write("dx %.14f\n" % dx)
    f.write("sampleMeshSize %8d\n"  % meshSize)
    f.write("stencilMeshSize %8d\n" % meshSize)
    f.write("stencilSize %2d\n" % stencilSize)
    f.write("nx %8d\n" % Nx)
    f.close()

  else:
    f.write("dim %1d\n" % 2)
    f.write("xMin %.14f\n" % min(xBd))
    f.write("xMax %.14f\n" % max(xBd))
    f.write("yMin %.14f\n" % min(yBd))
    f.write("yMax %.14f\n" % max(yBd))
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
    axFM.set_xlim(min(xBd), max(xBd))
    if dim!=1:
      axFM.set_ylim(min(yBd), max(yBd))

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
If you only pass one value, I assume 1d.\
If you pass two values, I assume 2d.")

  parser.add_argument(
    "-b", "--bounds",
    nargs="*", dest="bounds", type=float,
    help="Domain bounds along x and y. \
First, you pass all values of x, then pass all values of y.")

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
    dest="periodic", type=str2bool, default=True,
    help="True/False to enable/disable periodic BC so that mesh connectivity\
has that info embedded in.")

  parser.add_argument(
    "-d", "--debug",
    dest="debug", type=str2bool, default=False,
    help="True/False to eanble debug mode.")

  args = parser.parse_args()

  # -------------------------------------------------------
  # check that args make sense
  assert( len(args.numCells) == 1 or len(args.numCells) == 2 )
  assert( args.stencilSize in [3, 5, 7] )
  assert( args.orderingType in ["naturalRow", "naturalRowRcm"] )

  # check if working dir exists, if not, make it
  if not os.path.exists(args.outDir):
    os.system('mkdir -p ' + args.outDir)

  nx = int(args.numCells[0])
  ny = 1
  if len(args.numCells) == 2:
    ny = int(args.numCells[1])
  print(nx, ny)

  dim = -1
  if ny==1:
    if len(args.bounds)!=2:
      print("For 1d, you need to pass ny=1 and the domain's bounds")
      sys.exit()
    else:
      xCoords = [args.bounds[0], args.bounds[1]]
      yCoords = [0.0, 0.0]
      dim = 1

  elif nx!=1 and ny!=1 and len(args.bounds) == 4:
    xCoords = [args.bounds[0], args.bounds[1]]
    yCoords = [args.bounds[2], args.bounds[3]]
    dim = 2

  elif nx!=1 and ny!=1 and len(args.bounds) == 12:
    xCoords = []
    yCoords = []
    print(args.bounds)
    for i,v in enumerate(args.bounds):
      if i % 2 == 0:
        xCoords.append(v)
      else:
        yCoords.append(v)
    dim = 2

  # other things
  plotFontSize = 0
  if len(args.plottingInfo) == 2:
    plotFontSize = int(args.plottingInfo[1])

  main(args.outDir, args.debug,
       nx, ny, args.plottingInfo[0],
       args.orderingType,
       xCoords, yCoords,
       args.stencilSize, plotFontSize,
       args.periodic, dim)



  # coords = args.bounds
  # assert(numCoordsIn % 2 == 0)
  # xPoints = args.bounds[0:int(numCoordsIn/2)] if ny!=1 else args.bounds
  # if nx != 1:
  #   yPoints = args.bounds[int(numCoordsIn/2):]
  # checkDomainBoundsCoordinates(xPoints, yPoints)
  # print(xPoints)
