
import sys, os, time, glob, collections
import numpy as np

def readMeshConnec(fullMeshDir):
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
def readMeshInfo(fullMeshDir):
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

  return [dim, dx,dy,dz, int(sampleMeshSize), int(stencilMeshSize), bounds, stencilSize]

#====================================================================
def readMeshCoordinates(fullMeshDir,dim):
  indices, x,y,z = [], [], [], []
  inFile = fullMeshDir+"/coordinates.dat"
  with open(inFile) as fp:
    line = fp.readline()
    cnt = 1
    while line:
      lineList = line.split()
      indices.append(np.int32(lineList[0]))
      x.append(float(lineList[1]))
      y.append(float(lineList[2]))
      if dim==3:
        z.append(float(lineList[3]))

      line = fp.readline()
      cnt += 1
  return np.array(indices), np.array(x), np.array(y), np.array(z)
