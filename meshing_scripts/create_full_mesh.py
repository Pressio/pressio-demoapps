
import sys, os, time, collections
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

from _impl.one_dim_mesh import OneDimMesh
from _impl.natural_order_mesh_2d import NatOrdMeshRow2d
from _impl.natural_order_mesh_3d import NatOrdMesh3d
from _impl.remove_cells_step import *

#------------------------------------------------
def printDicPretty(d):
  for key, value in d.items():
    print(str(key), value)

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
def main_impl(workDir, debug, dimensionality, \
              Nx, Ny, Nz, orderingType, \
              xBd, yBd, zBd, stencilSize, \
              periodicX, periodicY, periodicZ):

  # figure out if domain is a plain rectangle or has a step in it
  plainDomain = True
  if dimensionality==2 and len(xBd) > 4:
    plainDomain = False
  if dimensionality==3:
    plainDomain = True

  if plainDomain:
    numCells = Nx*Ny*Nz
    L = [xBd[1]-xBd[0], yBd[1]-yBd[0], zBd[1]-zBd[0]]
    dx,dy,dz = L[0]/Nx, L[1]/Ny, L[2]/Nz
    print ("Bounds for x = ", xBd[0], xBd[1], " with dx = ", dx)
    if dimensionality>1:
      print ("Bounds for y = ", yBd[0], yBd[1], " with dy = ", dy)
    if dimensionality==3:
      print ("Bounds for z = ", zBd[0], zBd[1], " with dz = ", dz)

  # ---------------------------------------------
  x = None
  y = None
  z = None
  G = None
  gids = None
  if dimensionality == 1:
    meshObj = OneDimMesh(Nx, dx, xBd[0], xBd[1], stencilSize, periodicX)
    [x, y] = meshObj.getCoordinates()
    gids = meshObj.getGIDs()
    G = meshObj.getGraph()

  elif dimensionality==2 and plainDomain:
    meshObj = NatOrdMeshRow2d(Nx,Ny, dx,dy,\
                              xBd[0], xBd[1], \
                              yBd[0], yBd[1], \
                              stencilSize, periodicX, periodicY)

    # get mesh coordinates, gids and graph
    [x, y] = meshObj.getCoordinates()
    gids = meshObj.getGIDs()
    G = meshObj.getGraph()

  elif dimensionality==2 and not plainDomain:
    minX, maxX = min(xBd), max(xBd)
    minY, maxY = min(yBd), max(yBd)
    L = [maxX-minX, maxY-minY]
    dx,dy = L[0]/Nx, L[1]/Ny

    meshObj = NatOrdMeshRow2d(Nx, Ny, dx, dy,\
                              minX, maxX, minY, maxY, \
                              stencilSize, False)

    x,y,G,gids = removeStepCells(meshObj, xBd, yBd)

  elif dimensionality==3 and plainDomain:
    minX, maxX = min(xBd), max(xBd)
    minY, maxY = min(yBd), max(yBd)
    minZ, maxZ = min(zBd), max(zBd)
    L = [maxX-minX, maxY-minY, maxZ-minZ]
    dx,dy,dz = L[0]/Nx, L[1]/Ny, L[2]/Nz

    meshObj = NatOrdMesh3d(Nx, Ny, Nz, dx, dy, dz,\
                           xBd[0], xBd[1], \
                           yBd[0], yBd[1], \
                           zBd[0], zBd[1], \
                           stencilSize, \
                           periodicX, periodicY, periodicZ)

    # get mesh coordinates, gids and graph
    [x, y, z] = meshObj.getCoordinates()
    gids = meshObj.getGIDs()
    G = meshObj.getGraph()

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

  # -----------------------------------------------------
  meshSize = len(G)
  print ("meshSize = ", meshSize)

  # -----------------------------------------------------
  # info file
  f = open(workDir + "/info.dat","w+")
  if dimensionality==1:
    f.write("dim %1d\n" % 1)
    f.write("xMin %.14f\n" % xBd[0])
    f.write("xMax %.14f\n" % xBd[1])
    f.write("dx %.14f\n" % dx)
    f.write("sampleMeshSize %8d\n"  % meshSize)
    f.write("stencilMeshSize %8d\n" % meshSize)
    f.write("stencilSize %2d\n" % stencilSize)
    f.write("nx %8d\n" % Nx)
    f.close()

  elif dimensionality==2:
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

  elif dimensionality==3:
    f.write("dim %1d\n" % 3)
    f.write("xMin %.14f\n" % min(xBd))
    f.write("xMax %.14f\n" % max(xBd))
    f.write("yMin %.14f\n" % min(yBd))
    f.write("yMax %.14f\n" % max(yBd))
    f.write("zMin %.14f\n" % min(zBd))
    f.write("zMax %.14f\n" % max(zBd))
    f.write("dx %.14f\n" % dx)
    f.write("dy %.14f\n" % dy)
    f.write("dz %.14f\n" % dz)
    f.write("sampleMeshSize %8d\n"  % meshSize)
    f.write("stencilMeshSize %8d\n" % meshSize)
    f.write("stencilSize %2d\n" % stencilSize)
    f.write("nx %8d\n" % Nx)
    f.write("ny %8d\n" % Ny)
    f.write("nz %8d\n" % Nz)
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
    if dimensionality==3:
      f.write("%.14f " % z[k])
    f.write("\n")
  f.close()

###############################
if __name__== "__main__":
###############################
  # the RawTextHelpFormatter is needed so that help is formatted properly
  parser = ArgumentParser(description='help', formatter_class=RawTextHelpFormatter)

  parser.add_argument(
    "-d", "--debug",
    dest="debug",
    type=str2bool,
    default=False,
    choices=[False,True],
    help="Enable/disable debug mode.")

  parser.add_argument(
    "--outDir", "--outdir",
    dest="outDir",
    help="Full path to output directory where all the mesh files will be generated.")

  parser.add_argument(
    "-n", "--numCells",
    nargs="*",
    type=int,
    dest="numCells",
    help="Number of cells to use along each axis. This determines the dimensionality.\n"+
    "If you pass one value ,   I assume a 1d domain.\n"+
    "If you pass two values,   I assume a 2d domain.\n"+
    "If you pass three values, I assume a 3d domain.")

  parser.add_argument(
    "-b", "--bounds",
    nargs="*",
    type=float,
    dest="bounds",
    help="Domain physical bounds along each axis. Must be a list of pairs.\n"+
    "First, you pass bounds for x axis, then, if needed, for y, and z.\n"+
    "For example: \n"+
    "  --bounds 0.0 1.0           implies that x is in (0,1)\n"+
    "  --bounds 0.0 1.0 -1.0 1.0: implies x in (0,1) and y in (-1,1).\n"+
    "NOTE: the number of pairs passed to --bounds must be consistent \n"+
    "      with how many args you specify for --numCells. ")

  parser.add_argument(
    "-o", "--ordering",
    dest="orderingType",
    default="naturalRow",
    choices=["naturalRow"],
    help="Desired ordering to use for enumrating cells.\n"+
    "naturalRow: row natural ordering of the mesh cells (default)\n")
    #use <naturalRowRcm> to apply the reverse Cuthill-Mckee")

  parser.add_argument(
    "-s", "--stencilSize", "--stencilsize",
    type=int, dest="stencilSize", default=3,
    choices=[3,5,7],
    help="Stencil size to use for assembling the connectivity")

  parser.add_argument(
    "--periodic",
    nargs="*",
    type=str,
    dest="periodic",
    choices=["x", "y", "z", "x y", "x y z", "x z", "y z"],
    help="Specify which axes to assume periodic.\n"+
    "Choices are either 'x, y, z'. Default is nonperiodic.")

  args = parser.parse_args()

  # -------------------------------------------------------
  # check args
  assert( len(args.numCells) in [1,2,3] )
  assert( args.stencilSize   in [3,5,7] )
  assert( args.orderingType  in ["naturalRow", "naturalRowRcm"] )

  # check if working dir exists, if not, make it
  if not os.path.exists(args.outDir):
    os.system('mkdir -p ' + args.outDir)

  nx = int(args.numCells[0])
  ny = 1
  nz = 1
  if len(args.numCells) == 2:
    ny = int(args.numCells[1])
  if len(args.numCells) == 3:
    ny = int(args.numCells[1])
    nz = int(args.numCells[2])

  dim = -1
  if ny==1 and nz==1:
    if len(args.bounds)!=2:
      print("For 1d, you need to pass ny=1 and the domain's bounds")
      sys.exit()
    else:
      xCoords = [args.bounds[0], args.bounds[1]]
      yCoords = [0.0, 0.0]
      zCoords = [0.0, 0.0]
      dim = 1

  elif nx!=1 and ny!=1 and nz==1 and len(args.bounds) == 4:
    xCoords = [args.bounds[0], args.bounds[1]]
    yCoords = [args.bounds[2], args.bounds[3]]
    zCoords = [0.0, 0.0]
    dim = 2

  elif nx!=1 and ny!=1 and nz==1 and len(args.bounds) == 12:
    xCoords = []
    yCoords = []
    zCoords = [0.0, 0.0]
    print(args.bounds)
    for i,v in enumerate(args.bounds):
      if i % 2 == 0:
        xCoords.append(v)
      else:
        yCoords.append(v)
    dim = 2

  elif nx!=1 and ny!=1 and nz!=1:
    xCoords = [args.bounds[0], args.bounds[1]]
    yCoords = [args.bounds[2], args.bounds[3]]
    zCoords = [args.bounds[4], args.bounds[5]]
    dim = 3

  # handle periodic
  periodicX = False
  periodicY = False
  periodicZ = False
  ap = args.periodic

  if isinstance(ap, list):
    xIsInPeriodicArg = "x" in ap if isinstance(ap, list) else ap=="x"
    yIsInPeriodicArg = "y" in ap if isinstance(ap, list) else ap=="y"
    zIsInPeriodicArg = "z" in ap if isinstance(ap, list) else ap=="z"
  else:
    xIsInPeriodicArg = ap=="x"
    yIsInPeriodicArg = ap=="y"
    zIsInPeriodicArg = ap=="z"

  # for 1d, we cannot have y or z in periodic arg
  if dim==1 and yIsInPeriodicArg:
    sys.exit("For 1d problem you cannot specify --periodic {}".format("y"))
  if dim==1 and zIsInPeriodicArg:
    sys.exit("For 1d problem you cannot specify --periodic {}".format("z"))

  # for 2d, we cannot have z in periodic arg
  if dim==2 and zIsInPeriodicArg:
    sys.exit("For 2d problem you cannot specify --periodic {}".format("z"))

  periodicX = xIsInPeriodicArg
  periodicY = yIsInPeriodicArg
  periodicZ = zIsInPeriodicArg

  main_impl(args.outDir, args.debug, dim, \
            nx, ny, nz, args.orderingType, \
            xCoords, yCoords, zCoords, \
            args.stencilSize, \
            periodicX, periodicY, periodicZ)
