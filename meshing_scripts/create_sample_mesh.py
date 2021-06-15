
import matplotlib.pyplot as plt
import sys, os, time, glob
import numpy as np
from numpy import linspace, meshgrid
from matplotlib import cm
import collections
from argparse import ArgumentParser
import random
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

from mesh_utils import *

#=======================================================
def str2bool(v):
  if isinstance(v, bool):
    return v
  if v.lower() in ('yes', 'true', 't', 'y', '1'):
    return True
  elif v.lower() in ('no', 'false', 'f', 'n', '0'):
    return False
  else:
    raise argparse.ArgumentTypeError('Boolean value expected.')

#====================================================================
def readFullMeshConnec(fullMeshDir):
  G = {}
  gids = []
  inFile = fullMeshDir+"/connectivity.dat"
  with open(inFile) as fp:
    line = fp.readline()
    cnt = 1
    while line:
      lineList = line.split()
      numNeigh = len(lineList[1:])
      connections = [np.int64(i) for i in lineList[1:]]
      pt = np.int64(lineList[0])
      gids.append(pt)
      G[pt] = np.copy(connections)
      line = fp.readline()
      cnt += 1

  gids = np.array(gids)
  return [G,gids]

#====================================================================
def readFullMeshInfo(fullMeshDir):
  dim=0
  bounds=[None]*6
  dx,dy,dz = 0., 0., 0.
  sampleMeshSize = 0
  stencilMeshSize = 0
  stencilSize = 0
  infoFile = fullMeshDir+"/info.dat"

  with open(infoFile) as fp:
    line = fp.readline()
    cnt = 1
    while line:
      lineList = line.split()
      if "dim" in lineList:
        dim = int(lineList[1])

      elif "dx" in lineList:
        dx = float(lineList[1])

      elif "dy" in lineList:
        dy = float(lineList[1])

      elif "dz" in lineList:
        dz = float(lineList[1])

      elif "xMin" in lineList:
        bounds[0] = float(lineList[1])

      elif "xMax" in lineList:
        bounds[1] = float(lineList[1])

      elif "yMin" in lineList:
        bounds[2] = float(lineList[1])

      elif "yMax" in lineList:
        bounds[3] = float(lineList[1])

      elif "zMin" in lineList:
        bounds[4] = float(lineList[1])

      elif "zMax" in lineList:
        bounds[5] = float(lineList[1])

      elif "sampleMeshSize" in lineList:
        sampleMeshSize = np.int64(lineList[1])

      elif "stencilMeshSize" in lineList:
        stencilMeshSize = np.int64(lineList[1])

      elif "stencilSize" in lineList:
        stencilSize = np.int64(lineList[1])

      line = fp.readline()
      cnt += 1

  # in the full mesh, the sample and stencil size must be equal
  assert(sampleMeshSize == stencilMeshSize)
  return [dim, dx,dy,dz, int(sampleMeshSize), bounds, stencilSize]

#====================================================================
def readFullMeshCoordinates(fullMeshDir,dim):
  x,y,z = [], [], []
  inFile = fullMeshDir+"/coordinates.dat"
  with open(inFile) as fp:
    line = fp.readline()
    cnt = 1
    while line:
      lineList = line.split()
      x.append(float(lineList[1]))
      y.append(float(lineList[2]))
      if dim==3:
        z.append(float(lineList[3]))

      line = fp.readline()
      cnt += 1
  return np.array(x), np.array(y), np.array(z)


def initFigures(dim):
  figFM = plt.figure(0)
  figSM = plt.figure(1)
  if dim<=2:
    axFM = figFM.gca()
    axSM = figSM.gca()
  elif dim==3:
    axFM = figFM.gca(projection='3d')
    axSM = figSM.gca(projection='3d')
  return figFM, figSM, axFM, axSM

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


#=========================================================================
def main(workDir, debug, fullMeshDir, tilingDir, plotting, plotFontSize, darkMode):
#=========================================================================
  dim,dx,dy,dz,numCells,domainBounds,stencilSize = readFullMeshInfo(fullMeshDir)
  x,y,z = readFullMeshCoordinates(fullMeshDir, dim)
  G,gids = readFullMeshConnec(fullMeshDir)

  if (plotting != "none"):
    figFM, figSM, axFM, axSM = initFigures(dim)

  if plotting != "none" and dim==1:
    cellsColor  = 'w' if not darkMode else 'none'
    plotCells1d(x, dx, gids, axFM, darkMode, facecol=cellsColor, fontSz=plotFontSize)
    axFM.set_xlim(np.min(x)-dx*0.5, np.max(x)+dx*0.5)
    axFM.set_ylim(-0.02, 0.02)
    axFM.set_yticks([])

  if plotting != "none" and dim==2:
    cellsColor  = 'w' if not darkMode else 'none'
    plotCells2d(x, y, dx, dy, gids, axFM, darkMode, facecol=cellsColor, fontSz=plotFontSize)
    axFM.set_aspect(1.0)
    axFM.set_xlim(np.min(x)-dx*0.5, np.max(x)+dx*0.5)
    axFM.set_ylim(np.min(y)-dy*0.5, np.max(y)+dy*0.5)

  if plotting != "none" and dim==3:
    plotCells3d(x, y, z, dx, dy, dz, axFM, darkMode, facecol='w', fontSz=plotFontSize)

  if debug:
    print("natural order full mesh connectivity")
    printDicPretty(G)
    print("\n")

  # -----------------------------------------------------
  # read list of GIDs for the sample mesh
  # -----------------------------------------------------
  sampleMeshGIDs = np.loadtxt(workDir+"/sample_mesh_gids.dat", dtype=np.int64)
  sampleMeshGIDs = np.sort(sampleMeshGIDs)
  if debug:
    print("Sample mesh gids:")
    print(sampleMeshGIDs)

  # -----------------------------------------------------
  # store the connectivity for the selected points
  # -----------------------------------------------------
  smGraph0 = collections.OrderedDict()
  for rPt in sampleMeshGIDs:
    smGraph0[rPt] = G[rPt]
  print("\n")
  if debug:
    print("sample mesh graph0 (IDs wrt full mesh)")
    for k,v in smGraph0.items(): print(k, v)

  stencilMeshGIDs = []
  # loop over target cells wherw we want residual
  for k,v in smGraph0.items():
    # append the GID of this cell
    stencilMeshGIDs.append(k)
    # append GID of stencils/neighborin cells
    for j in v:
      if j != -1:
        stencilMeshGIDs.append(np.int64(j))

  # remove duplicates and sort
  stencilMeshGIDs = list(dict.fromkeys(stencilMeshGIDs))
  stencilMeshGIDs.sort()
  if debug:
    print("\n")
    print("stencilMeshGIDs (IDs wrt full mesh)")
    print(stencilMeshGIDs)

  np.savetxt(workDir+"/stencil_mesh_gids.dat", stencilMeshGIDs, fmt="%8d")

  # -----------------------------------------------------
  # if we get here, the sample mesh graph contains the GIDs
  # of the sample mesh wrt the FULL mesh ordering becaus we
  # simply extracted a subset from the full mesh.
  # However, what we need is to reenumerate the sample mesh cells.
  # We numerate the sample mesh points with new indexing
  # and create a map of full-mesh gids to new gids
  # -----------------------------------------------------
  # fm_to_sm_map is such that:
  # - key   = GID wrt full mesh
  # - value = GID wrt sample mesh

  if tilingDir!=None:
    gidsfiles = glob.glob(tilingDir+"/cell_gids_p_*.txt")

    # sort based on the ID, so need to extract ID which is last of dir name
    def func(elem): return int(elem.split('_')[-1].split('.')[0])
    gidsfiles = sorted(gidsfiles, key=func)
    numPartitions = len(gidsfiles)

    fm_to_sm_map = collections.OrderedDict()
    i=0
    for p in range(numPartitions):
      print(p)
      pgids = np.loadtxt(gidsfiles[p], dtype=int)
      b = list(set(stencilMeshGIDs).intersection(pgids))
      for pt in np.sort(np.array(b)):
        fm_to_sm_map[pt] = i
        i+=1
  else:
    fm_to_sm_map = collections.OrderedDict()
    i = 0
    for pt in stencilMeshGIDs:
      fm_to_sm_map[pt] = i
      i+=1

  if debug:
    print("Done with fm_to_sm_map")
    for k,v in fm_to_sm_map.items(): print(k, v)

  if debug:
    print("doing now the sm -> fm gids mapping ")
  sm_to_fm_map = collections.OrderedDict()
  for k,v in fm_to_sm_map.items():
    sm_to_fm_map[v] = k
  if debug:
    print("Done with sm_to_fm_map")
    for k,v in sm_to_fm_map.items():
      print(k, v)

  # -----------------------------------------------------
  # Here we have a list of unique GIDs for the sample mesh.
  # Map the GIDs from the full to sample mesh indexing
  # so that we know what is what.
  # -----------------------------------------------------
  sampleMeshGraph = collections.OrderedDict()
  for rGidFM, v in smGraph0.items():
    smGID = fm_to_sm_map[rGidFM]
    smStencilGIDs = v
    for i in range(len(smStencilGIDs)):
      thisGID = smStencilGIDs[i]
      if thisGID != -1:
        smStencilGIDs[i] = fm_to_sm_map[thisGID]
    sampleMeshGraph[smGID] = smStencilGIDs

  if debug:
    print("\n")
    print("Done with sampleMeshGraph")
    print("sample mesh connectivity")
    printDicPretty(sampleMeshGraph)

  gids_sm = list(sm_to_fm_map.keys())
  if plotting != "none":
    plotReducedMesh(axSM, dim, plotFontSize, \
                    x, y, z, dx, dy, dz, \
                    gids_sm, fm_to_sm_map,
                    stencilMeshGIDs, sampleMeshGIDs,\
                    darkMode)

  # -----------------------------------------------------
  sampleMeshSize = len(sampleMeshGraph)
  print ("sampleMeshSize = ", sampleMeshSize,
         " which is = ", sampleMeshSize/numCells*100, " % of full mesh")
  stencilMeshSize = len(fm_to_sm_map)
  print ("stencilMeshSize = ", stencilMeshSize,
         " which is = ", stencilMeshSize/numCells*100, " % of full mesh")

  # -----------------------------------------------------
  # print info file
  f = open(workDir+"/info.dat","w+")
  if dim==1:
    f.write("dim %1d\n" % dim)
    f.write("xMin %.14f\n" % domainBounds[0])
    f.write("xMax %.14f\n" % domainBounds[1])
    f.write("dx %.14f\n" % dx)
    f.write("sampleMeshSize %8d\n"  % sampleMeshSize)
    f.write("stencilMeshSize %8d\n" % stencilMeshSize)
    f.write("stencilSize %2d\n" % stencilSize)
    f.close()

  elif dim==2:
    f.write("dim %1d\n" % dim)
    f.write("xMin %.14f\n" % domainBounds[0])
    f.write("xMax %.14f\n" % domainBounds[1])
    f.write("yMin %.14f\n" % domainBounds[2])
    f.write("yMax %.14f\n" % domainBounds[3])
    f.write("dx %.14f\n" % dx)
    f.write("dy %.14f\n" % dy)
    f.write("sampleMeshSize %8d\n"  % sampleMeshSize)
    f.write("stencilMeshSize %8d\n" % stencilMeshSize)
    f.write("stencilSize %2d\n" % stencilSize)
    f.close()

  elif dim==3:
    f.write("dim %1d\n" % dim)
    f.write("xMin %.14f\n" % domainBounds[0])
    f.write("xMax %.14f\n" % domainBounds[1])
    f.write("yMin %.14f\n" % domainBounds[2])
    f.write("yMax %.14f\n" % domainBounds[3])
    f.write("zMin %.14f\n" % domainBounds[4])
    f.write("zMax %.14f\n" % domainBounds[5])
    f.write("dx %.14f\n" % dx)
    f.write("dy %.14f\n" % dy)
    f.write("dz %.14f\n" % dz)
    f.write("sampleMeshSize %8d\n"  % sampleMeshSize)
    f.write("stencilMeshSize %8d\n" % stencilMeshSize)
    f.write("stencilSize %2d\n" % stencilSize)
    f.close()

  # -----------------------------------------------------
  # print connectivity file
  f = open(workDir+"/connectivity.dat","w+")
  for k in sorted(sampleMeshGraph.keys()):
    f.write("%8d " % k)
    for i in sampleMeshGraph[k]:
      f.write("%8d " % i)
    f.write("\n")
  f.close()

  # -----------------------------------------------------
  # print coordinate file
  f = open(workDir+"/coordinates.dat","w+")
  for k,v in sm_to_fm_map.items():
    f.write("%8d " % k)
    f.write("%.14f " % x[v])
    f.write("%.14f " % y[v])
    f.write("\n")
  f.close()

  # # -----------------------------------------------------
  # # print gids mapping: mapping between GIDs in sample mesh to full mesh
  # # I only need this when I use the sample mesh
  # # -----------------------------------------------------
  # if (samplingType=="random"):
  #   f1 = open("sm_to_fm_gid_mapping.dat","w+")
  #   for k,v in fm_to_sm_map.items():
  #     f1.write("%8d %8d\n" % (v, k))
  #   f1.close()

  if plotting == "show":
    plt.show()
  elif plotting =="print":
    if darkMode==1:
      figFM.savefig(workDir+'/full.png'  ,dpi=250,format="png", transparent=True)
      figSM.savefig(workDir+'/sample.png',dpi=250,format="png", transparent=True)
    else:
      figFM.savefig(workDir+'/full.png'  ,dpi=250,format="png")
      figSM.savefig(workDir+'/sample.png',dpi=250,format="png")


###############################
if __name__== "__main__":
###############################
  parser = ArgumentParser()

  parser.add_argument("-wdir", "--wdir", dest="wdir")

  parser.add_argument("-p", "-plotting", "--plotting",
                      nargs="*",
                      dest="plottingInfo",
                      default="none",
                      help="What type of plotting you want:\n"+
                           "use <show> for showing plots,\n"+
                           "use <print> for printing only, \n"+
                           "use <none> for no plots")

  parser.add_argument("-d", "-debug", "--debug",
                      type=bool, dest="debug", default=False)

  parser.add_argument("-fullMeshDir", "--fullMeshDir",
                      dest="fullMeshDir", default=None)

  parser.add_argument("-useTilingFrom", "--useTilingFrom",
                      dest="tilingDir", default=None)
  args = parser.parse_args()

  assert(args.fullMeshDir != None)

  plotFontSize = 0
  darkMode = 0
  if len(args.plottingInfo) >= 2:
    plotFontSize = int(args.plottingInfo[1])
  if len(args.plottingInfo) == 3:
    darkMode = int(args.plottingInfo[2])

  main(args.wdir, args.debug, \
       args.fullMeshDir, args.tilingDir, \
       args.plottingInfo[0], plotFontSize, darkMode)
