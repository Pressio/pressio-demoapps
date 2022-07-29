#!/usr/bin/env python

import sys, os, time
import numpy as np
from scipy import linalg

from load_snapshots import \
  load_state_snapshot_matrix, \
  load_rhs_snapshot_matrix_if_present
from extract_from_run_directory import \
  find_mesh_from_run_directory


def _load_basis_from_binary_file(lsvFile):
  nr, nc  = np.fromfile(lsvFile, dtype=np.int64, count=2)
  fullPhi = np.fromfile(lsvFile, offset=np.dtype(np.int64).itemsize*2)
  fullPhi = np.reshape(fullPhi, (nr, nc), order='F')
  return fullPhi


def do_svd_py(dryRun, mymatrix, svaFile, lsvFile):
  timing = np.zeros(1)
  start = time.time()
  U,S,_ = linalg.svd(mymatrix, full_matrices=False, lapack_driver='gesdd')
  end = time.time()
  elapsed = end - start
  timing[0] = elapsed
  print("elapsed ", elapsed)

  # singular values
  print("Writing sing values to file: {}".format(svaFile))
  np.savetxt(svaFile, S)

  assert(U.flags['F_CONTIGUOUS'])

  # left singular vectors
  fileo = open(lsvFile, "wb")
  # write to beginning of file the extents of the matrix
  print("Writing POD modes to file:   {}".format(lsvFile))
  r=np.int64(U.shape[0])
  np.array([r]).tofile(fileo)
  c=np.int64(U.shape[1])
  np.array([c]).tofile(fileo)
  '''
  NOTE: tofile write an array in rowwise, REGARDLESS of the layout of the matrix.
  So here we need to pass U.T to tofile so that tofile writes U in the proper
  way required format for reading in our cpp code
  '''
  UT = np.transpose(U)
  UT.tofile(fileo)
  fileo.close()

  # finally print timings to file
  # find out the directory where we are storing things
  outDir = os.path.dirname(lsvFile)
  np.savetxt(outDir+'/timings.txt', timing)


def _compute_state_pod_and_ref_state_dont_split_vars_impl(dryRun, \
                                                         numDofsPerCell, \
                                                         dataDirs,\
                                                         singValsFile,\
                                                         leftSingVecsFile, \
                                                         referenceStateFile):

  # use first run to find the mesh, and ensure all use the same mesh
  fullMeshPath = find_mesh_from_run_directory(dataDirs[0])
  for itDir in dataDirs:
    assert(fullMeshPath == find_mesh_from_run_directory(itDir))

  M = load_state_snapshot_matrix(dryRun, dataDirs, numDofsPerCell)
  # M is such that:
  # num of rows = num of state instances that the FOM saved
  # num of cols = numCells*numDofsPerCell
  numSnaps, totDofs = M.shape[0], M.shape[1]
  numCells = np.int64(totDofs/numDofsPerCell)

  # write reference state (zero for now)
  np.savetxt(referenceStateFile, np.zeros((totDofs, 1)))

  # note that we need to pass M.T here to have correct shape
  do_svd_py(dryRun, M.T, singValsFile, leftSingVecsFile)



def _compute_rhs_pod_if_needed_dont_split_vars_impl(dryRun, \
                                                    numDofsPerCell,\
                                                    dataDirs,
                                                    singValsFile,\
                                                    leftSingVecsFile):
  # load
  M = load_rhs_snapshot_matrix_if_present(dryRun, dataDirs, numDofsPerCell)
  # some dirs don't have data in it, so skip
  if M.any() == None:
    print("Some dirs don't have RHS snapshots, so skipping computing POD for RHS")

  else:
    # M is such that:
    # num of rows = num of state instances saved
    # num of cols = numSampleMeshCells*numDofsPerCell
    numSnaps, totDofs = M.shape[0], M.shape[1]
    numCells = np.int64(totDofs/numDofsPerCell)

    # note that we need to pass M.T here
    do_svd_py(dryRun, M.T, singValsFile, leftSingVecsFile)



def _compute_state_pod_and_ref_state_split_vars_impl(dryRun, \
                                                     numDofsPerCell, \
                                                     dataDirs, \
                                                     singValsFile,\
                                                     leftSingVecsFile):
  sys.exit(__file__ + " missing impl")

#   # use first run to find the mesh, and ensure all use the same mesh
#   fullMeshPath = find_mesh_from_run_directory(dataDirs[0])
#   for itDir in dataDirs:
#     assert(fullMeshPath == find_mesh_from_run_directory(itDir))

#   M = load_state_snapshot_matrix(dryRun, dataDirs, module.numDofsPerCell)
#   # M is such that:
#   # num of rows = num of state instances that the FOM saved
#   # num of cols = numCells*numDofsPerCell
#   numSnaps, totDofs = M.shape[0], M.shape[1]
#   numCells = np.int64(totDofs/module.numDofsPerCell)

#   # write reference state (zero for now)
#   refStateFile = outDir + "/referenceState"
#   np.savetxt(refStateFile, np.zeros((totDofs, 1)))

#   numVars = module.numDofsPerCell
#   M = np.reshape(M, (numSnaps, numCells, module.numDofsPerCell))
#   phiA = [None]*numVars
#   for i in range(0,numVars):
#     currS = np.rollaxis(M[:,:,i],1)
#     singValsFileName = "sva_dof"+str(i)
#     lsvFileName      = "lsv_dof"+str(i)
#     do_svd_py(dryRun, currS, outDir, singValsFileName, lsvFileName)



def _compute_rhs_pod_if_needed_split_vars_impl(dryRun, \
                                               numDofsPerCell, \
                                               dataDirs, \
                                               singValsFile,\
                                               leftSingVecsFile):

  sys.exit(__file__ + " missing impl")

#   # load
#   M = load_rhs_snapshot_matrix_if_present(dryRun, dataDirs, module.numDofsPerCell)
#   # some dirs don't have data in it, so skip
#   if M.any() == None:
#     print("Some dirs don't have RHS snapshots, so skipping computing POD for RHS")

#   else:
#     # M is such that:
#     # num of rows = num of state instances saved
#     # num of cols = numSampleMeshCells*numDofsPerCell
#     numSnaps, totDofs = M.shape[0], M.shape[1]
#     numCells = np.int64(totDofs/module.numDofsPerCell)

#     numVars = module.numDofsPerCell
#     M = np.reshape(M, (numSnaps, numCells, module.numDofsPerCell))
#     phiA = [None]*numVars
#     for i in range(0,numVars):
#       currS = np.rollaxis(M[:,:,i],1)
#       singValsFile = "sva_dof"+str(i)
#       lsvFile      = "lsv_dof"+str(i)
#       do_svd_py(dryRun, currS, outDir, singValsFile, lsvFile)
