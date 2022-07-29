
#==============================================================
# imports
#==============================================================

import os, random, sys, subprocess
import numpy as np

try:
  import pressiotools.linalg as ptla
  from pressiotools.samplemesh.withLeverageScores import computeNodes
except ImportError:
  raise ImportError("Unable to import classes from pressiotools")

from generate_sample_mesh_files import \
  _generate_sample_mesh_files
from extract_from_mesh_directory \
  import find_num_cells_from_info_file
from extract_from_run_directory \
  import find_mesh_from_run_directory
from load_snapshots import \
  load_rhs_snapshot_matrix_if_present, \
  load_state_snapshot_matrix

#==============================================================
# functions
#==============================================================

def _make_random_sample_mesh_impl(dryRun, workDir, fraction, \
                                  outDir, dimensionality):
  # given the workdir, find the full mesh used
  fullMeshFullPaths = [workDir+'/'+d for d in os.listdir(workDir) \
                       if "full_mesh" in d]
  # ensure that only one mesh is found
  if len(fullMeshFullPaths) != 1:
    sys.exit("Found multiple full_mesh subdirs inside {}, should only have one".format(workDir))
  fullMeshPath = fullMeshFullPaths[0]

  if dimensionality==1:
    nx = find_num_cells_from_info_file(fullMeshPath, "nx")
    numCells = nx

  elif dimensionality==2:
    nx = find_num_cells_from_info_file(fullMeshPath, "nx")
    ny = find_num_cells_from_info_file(fullMeshPath, "ny")
    numCells = nx*ny

  elif dimensionality==3:
    nx = find_num_cells_from_info_file(fullMeshPath, "nx")
    ny = find_num_cells_from_info_file(fullMeshPath, "ny")
    nz = find_num_cells_from_info_file(fullMeshPath, "nz")
    numCells = nx*ny*nz

  # figure out how many cells I need to keep
  sampleMeshCount = int(numCells*fraction)

  smIndices = random.sample(range(0, numCells), sampleMeshCount)
  smIndices = np.sort(smIndices)
  smIndicesFile = outDir+'/sample_mesh_gids.dat'
  np.savetxt(smIndicesFile, smIndices, fmt='%10d')

  # compute the actual sample mesh files
  _generate_sample_mesh_files(fullMeshPath, outDir, smIndicesFile)



def _make_leverage_scores_sample_mesh_impl(dryRun, workDir, \
                                           setId, dataDirs, \
                                           listOfFractions, \
                                           listOfBetas, \
                                           numDofsPerCell):

  # ensure that all the train runs use the same FULL mesh
  # otherwise doing sample mesh makes no sense
  fullMeshPath = find_mesh_from_run_directory(dataDirs[0])
  for itDir in dataDirs:
    if fullMeshPath != find_mesh_from_run_directory(itDir):
      sys.exit("Some FOM train directories have inconsistent meshes")

  # try to load snapshots of the RHS
  M = load_rhs_snapshot_matrix_if_present(dryRun, dataDirs, numDofsPerCell)
  if M.any() == None:
    print("Cannot find rhs snapshots, so use state")
    M = load_state_snapshot_matrix(dryRun, dataDirs, numDofsPerCell)

  # transpose here so that we have the right shape needed below
  M = M.T
  numCells = int(M.shape[0]/numDofsPerCell)
  print("Snapshot shape  = {}".format(M.shape))
  print("Total num cells = {}".format(numCells))

  assert(M.flags['F_CONTIGUOUS'])

  for fraction in listOfFractions:
    # figure out how many cells I need to keep
    sampleMeshCount = int(numCells*fraction)

    for beta in listOfBetas:
      fractionString = "{:4.3f}".format(fraction)
      betaString     = "{:6.4f}".format(beta)
      smDirFullPath  = workDir+'/sample_mesh_lev_scores_fraction_' + fractionString
      smDirFullPath += '_beta_' + betaString + '_set_'+str(setId)

      if not os.path.exists(smDirFullPath):
        os.mkdir(smDirFullPath)

        # here we have one dof/cell so # of rows = num of cells
        Mv = ptla.MultiVector(M)
        smIndices, pmf = computeNodes(matrix=Mv, beta=beta, \
                                      numSamples=sampleMeshCount, \
                                      dofsPerMeshNode=numDofsPerCell)

        smIndicesFile = smDirFullPath+'/sample_mesh_gids.dat'
        np.savetxt(smIndicesFile, smIndices, fmt='%10d')

        # compute the actual sample mesh files
        _generate_sample_mesh_files(fullMeshPath, smDirFullPath, smIndicesFile)
      else:
        print("{} already exists".format(smDirFullPath))



# def make_sample_mesh_via_p_sampling(fullMeshPath, workDir, U):
#   currentSampleMeshCount = int(U.shape[1])

#   fractionString = "{:4.3f}".format(currentSampleMeshCount)
#   smDirFullPath  = workDir+'/sample_mesh_psampling_fraction_' + fractionString
#   smDirFullPath += '_set'+str(setId)

#   # note that for P sampling we cannot select more cells than # of snapshots
#   if not os.path.exists(smDirFullPath):
#     print("Sample mesh via p-sampling with {} cells".format(currentSampleMeshCount))
#     os.mkdir(smDirFullPath)

#     Q,R,P         = linalg.qr(U[:,0:currentSampleMeshCount].transpose(), pivoting=True )
#     smIndices     = np.sort(P[0:currentSampleMeshCount])
#     smIndicesFile = smDirFullPath+'/sample_mesh_gids.dat'
#     np.savetxt(smIndicesFile, smIndices, fmt='%10d')
#     # compute the actual sample mesh files
#     generate_sample_mesh(fullMeshPath, smDirFullPath, smIndicesFile)
