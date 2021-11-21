
import sys, os, time, glob, collections
from argparse import ArgumentParser
import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

from _impl.read_mesh_files import *

#=========================================================================
def printDicPretty(d):
  for key, value in d.items():
    print(str(key), value)

#=========================================================================
def str2bool(v):
  if isinstance(v, bool):
    return v
  if v.lower() in ('yes', 'true', 't', 'y', '1'):
    return True
  elif v.lower() in ('no', 'false', 'f', 'n', '0'):
    return False
  else:
    raise argparse.ArgumentTypeError('Boolean value expected.')

#=========================================================================
def main(workDir, debug, fullMeshDir, tilingDir, smGIDsFile):
  dim,dx,dy,dz,samplesMeshSize, stencilMeshSize, domainBounds,stencilSize = readMeshInfo(fullMeshDir)
  _,x,y,z = readMeshCoordinates(fullMeshDir, dim)
  G,gids = readMeshConnec(fullMeshDir)

  if debug:
    print("natural order full mesh connectivity")
    printDicPretty(G)
    print("\n")

  # read list of GIDs for the sample mesh
  sampleMeshGIDs = np.loadtxt(smGIDsFile, dtype=np.int64)
  sampleMeshGIDs = np.sort(sampleMeshGIDs)
  if debug:
    print("Sample mesh gids:")
    print(sampleMeshGIDs)

  # store the connectivity for the selected points
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

  sampleMeshSize = len(sampleMeshGraph)
  print ("sampleMeshSize = ", sampleMeshSize,
         " which is = ", sampleMeshSize/samplesMeshSize*100, " % of full mesh")
  stencilMeshSize = len(fm_to_sm_map)
  print ("stencilMeshSize = ", stencilMeshSize,
         " which is = ", stencilMeshSize/samplesMeshSize*100, " % of full mesh")

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

  # print connectivity file
  f = open(workDir+"/connectivity.dat","w+")
  for k in sorted(sampleMeshGraph.keys()):
    f.write("%8d " % k)
    for i in sampleMeshGraph[k]:
      f.write("%8d " % i)
    f.write("\n")
  f.close()

  # print coordinate file
  f = open(workDir+"/coordinates.dat","w+")
  for k,v in sm_to_fm_map.items():
    f.write("%8d " % k)
    f.write("%.14f " % x[v])
    f.write("%.14f " % y[v])
    if dim==3:
      f.write("%.14f " % z[v])
    f.write("\n")
  f.close()

###############################
if __name__== "__main__":
###############################
  parser = ArgumentParser()

  parser.add_argument(
    "-d", "-debug", "--debug",
    type=bool, dest="debug", default=False)

  parser.add_argument(
    "--sampleMeshIndices",
    dest="sampleMeshGIDsFile", default=None,
    help="Full path to text file with indices of cells to use as sample mesh celles")

  parser.add_argument(
    "-o", "--outDir", "--outdir", "--workDir",
    dest="wdir",
    help="Full path to where to store all the mesh output files.")

  parser.add_argument(
    "--fullMeshDir", "--fullmeshdir", "--fullmeshDir", "--fullMeshdir",
    dest="fullMeshDir", default=None)

  parser.add_argument(
    "--useTilingFrom",
    dest="tilingDir", default=None)

  args = parser.parse_args()
  # ------------------------------------

  assert(args.fullMeshDir != None)

  # by default, we look for file with GIDs in workdir
  smGIDsFile = args.wdir + "/sample_mesh_gids.dat",
  # if file is specified by user, use that instead
  if (args.sampleMeshGIDsFile != None):
    smGIDsFile = args.sampleMeshGIDsFile

  # check if working dir exists, if not, make it
  if not os.path.exists(args.wdir):
    os.system('mkdir -p ' + args.wdir)

  print(args.sampleMeshGIDsFile)
  main(args.wdir, args.debug, args.fullMeshDir, args.tilingDir, smGIDsFile)
