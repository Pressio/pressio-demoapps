
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
      #np.zeros(numNeigh, dtype=np.int64)
      #for i in range(numNeigh):
      #connections[i] = np.int64(lineList[1+i])
      pt = np.int64(lineList[0])
      gids.append(pt)
      G[pt] = np.copy(connections)

      line = fp.readline()
      cnt += 1

  gids = np.array(gids)
  return [G,gids]

#====================================================================
def readFullMeshInfo(fullMeshDir):
  bounds=[None]*4
  dx,dy = 0., 0.
  sampleMeshSize = 0
  stencilMeshSize = 0
  stencilSize = 0
  infoFile = fullMeshDir+"/info.dat"
  with open(infoFile) as fp:
    line = fp.readline()
    cnt = 1
    while line:
      lineList = line.split()
      if "dx" in lineList:
        dx = float(lineList[1])
      elif "dy" in lineList:
        dy = float(lineList[1])
      elif "xMin" in lineList:
        bounds[0] = float(lineList[1])
      elif "xMax" in lineList:
        bounds[1] = float(lineList[1])
      elif "yMin" in lineList:
        bounds[2] = float(lineList[1])
      elif "yMax" in lineList:
        bounds[3] = float(lineList[1])
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
  return [dx,dy, int(sampleMeshSize), bounds, stencilSize]

#====================================================================
def readFullMeshCoordinates(fullMeshDir):
  x,y = [], []
  inFile = fullMeshDir+"/coordinates.dat"
  with open(inFile) as fp:
    line = fp.readline()
    cnt = 1
    while line:
      lineList = line.split()
      x.append(float(lineList[1]))
      y.append(float(lineList[2]))

      line = fp.readline()
      cnt += 1
  return np.array(x), np.array(y)

#====================================================================
def main(workDir, debug, fullMeshDir, tilingDir, plotting, plotFontSize):
  dx,dy,numCells,domainBounds,stencilSize = readFullMeshInfo(fullMeshDir)
  x,y = readFullMeshCoordinates(fullMeshDir)
  G,gids = readFullMeshConnec(fullMeshDir)

  if (plotting != "none"):
    figFM = plt.figure(0)
    axFM = figFM.gca()
    figSM = plt.figure(1)
    axSM = figSM.gca()

  if plotting != "none":
    plotLabels(x, y, dx, dy, gids, axFM, fontSz=plotFontSize)
    axFM.set_aspect(1.0)
    axFM.set_xlim(np.min(x)-0.5,np.max(x)+0.5)
    axFM.set_ylim(np.min(y)-0.5,np.max(y)+0.5)

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

      # for pt in stencilMeshGIDs:
      #   if pt in pgids:
      #     fm_to_sm_map[pt] = i
      #     i+=1
  else:
    fm_to_sm_map = collections.OrderedDict()
    i = 0
    for pt in stencilMeshGIDs:
      fm_to_sm_map[pt] = i
      i+=1

  print("Done with fm_to_sm_map")
  if debug:
    for k,v in fm_to_sm_map.items(): print(k, v)


  print("doing now the sm -> fm gids mapping ")
  sm_to_fm_map = collections.OrderedDict()
  for k,v in fm_to_sm_map.items():
    sm_to_fm_map[v] = k
  print("Done with sm_to_fm_map")
  if debug:
    for k,v in sm_to_fm_map.items(): print(k, v)

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
          smStencilGIDs[i] = fm_to_sm_map[thisGID]
      sampleMeshGraph[smGID] = smStencilGIDs

  print("\n")
  print("Done with sampleMeshGraph")
  print("sample mesh connectivity")
  if debug:
    printDicPretty(sampleMeshGraph)

  gids_sm = list(sm_to_fm_map.keys())
  if plotting != "none":
    x2, y2 = x[ list(fm_to_sm_map.keys()) ], y[ list(fm_to_sm_map.keys()) ]
    plotLabels(x2, y2, dx, dy, gids_sm, axSM, 's', 'r', 0, fontSz=plotFontSize)
    axSM.set_aspect(1.0)
    axSM.set_xlim(np.min(x)-0.5,np.max(x)+0.5)
    axSM.set_ylim(np.min(y)-0.5,np.max(y)+0.5)

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
    figFM.savefig(workDir+'/full.pdf'  ,dpi=120,format="pdf")
    figSM.savefig(workDir+'/sample.pdf',dpi=120,format="pdf")


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
  if len(args.plottingInfo) == 2:
    plotFontSize = int(args.plottingInfo[1])

  main(args.wdir, args.debug, \
       args.fullMeshDir, args.tilingDir, \
       args.plottingInfo[0], plotFontSize)
