#!/usr/bin/env python

import sys, os, time, shutil, subprocess
import numpy as np
from extract_from_run_directory import \
  find_total_dofs_from_run_directory, \
  find_total_dofs_sample_mesh_from_run_directory

def load_state_snapshot_matrix(dryRun, trainDirs, numDofsPerCell):
  numTotDofs = find_total_dofs_from_run_directory(trainDirs[0], numDofsPerCell)
  M = np.zeros((0,numTotDofs))

  # loop over dirs to read data from
  for targetDirec in trainDirs:
    print("reading data from {}".format(targetDirec))
    # make sure current dir has consisten num of total dofs
    assert( numTotDofs == find_total_dofs_from_run_directory(targetDirec, numDofsPerCell))

    data = np.fromfile(targetDirec+"/fom_snaps_state.bin")
    numTimeSteps = int(np.size(data)/numTotDofs)
    D = np.reshape(data, (numTimeSteps, numTotDofs))
    M = np.append(M, D, axis=0)

  # mean = np.mean(M, axis=0)
  # for i in range(D.shape[0]):
  #   M[i,:] -= mean

  print("state snapshots: shape  : ", M.shape)
  print("state snapshots: min/max: ", np.min(M), np.max(M))
  return M

def load_rhs_snapshot_matrix_if_present(dryRun, trainDirs, numDofsPerCell):
  # first, check that all directories have the data inside
  allBool = [os.path.exists(dirIt+'/fom_snaps_rhs.bin') for dirIt in trainDirs]
  if not all(allBool):
    return np.array([None])
  else:
    # ensure all directories have the data with consistent size
    # all we need to check is that the tot dofs of RHS is consistent
    totRhsSize = find_total_dofs_sample_mesh_from_run_directory(trainDirs[0], numDofsPerCell)
    for dirIt in trainDirs:
      currentDirRhsSize = find_total_dofs_sample_mesh_from_run_directory(dirIt, numDofsPerCell)
      assert( totRhsSize == currentDirRhsSize )

    # if here, then data exists and is consistent, so load it
    M = np.zeros((0, totRhsSize))
    for dirIt in trainDirs:
      print("reading RHS data from {}".format(dirIt))
      data = np.fromfile(dirIt+"/fom_snaps_rhs.bin")
      numInstances = int(np.size(data)/totRhsSize)
      D = np.reshape(data, (numInstances, totRhsSize))
      M = np.append(M, D, axis=0)

    print("rhs snapshots: shape  : ", M.shape)
    print("rhs snapshots: min/max: ", np.min(M), np.max(M))
    return M
