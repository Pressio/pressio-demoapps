
#==============================================================
# imports
#==============================================================

import os, subprocess
import numpy as np

from extract_from_run_directory import \
  find_pod_modes_file_from_run_directory, \
  find_num_modes_from_run_directory, \
  find_reference_state_file_from_run_directory

#==============================================================
# functions
#==============================================================

def _reconstruct_final_state_on_full_mesh_impl(dryRun, workDir, romRunDir):

  # find pod modes file for FULL mesh
  lsvFile = find_pod_modes_file_from_run_directory(romRunDir)
  nModes  = find_num_modes_from_run_directory(romRunDir)

  # load reference state
  refStateFile = find_reference_state_file_from_run_directory(romRunDir)
  refState     = np.loadtxt(refStateFile)

  # extract the FINAL ROM state
  snaps = np.fromfile(romRunDir+"/rom_snaps.bin")
  nt    = int(np.size(snaps)/nModes)
  data  = np.reshape(snaps, (nt, nModes))
  finalRomState = data[-1,:]
  np.savetxt(romRunDir+"/final_state.txt", finalRomState)

  # reconstruct final state
  nr, nc  = np.fromfile(lsvFile, dtype=np.int64, count=2)
  fullPhi = np.fromfile(lsvFile, offset=np.dtype(np.int64).itemsize*2)
  fullPhi = np.reshape(fullPhi, (nr, nc), order='F')
  reconstructed = np.dot(fullPhi[:,:nModes], finalRomState) + refState
  finalReconstructedStateFile = romRunDir+"/final_reconstructed_state.txt"
  np.savetxt(finalReconstructedStateFile, reconstructed)


def _compute_final_state_error_impl(dryRun, fomTestDir, romDir):
  fomFile = fomTestDir+"/final_state.txt"
  romFile = romDir+"/final_reconstructed_state.txt"

  args   = ("python3",
            os.path.dirname(__file__) + "/state_error.py",
            "--skiprows", str(0),
            "--dryrun", str(dryRun),
            "--fomstate", fomFile,
            "--approxstate", romFile)

  if not dryRun:
    print("Running: {}".format(args))
    logFileFullPath = romDir+'/finalError.txt'
    logfile = open(logFileFullPath, 'w')
    p  = subprocess.Popen(args, stdout=logfile, stderr=logfile)
    p.wait()
    logfile.close()
  else:
    print("dryrun:rom_err: {}".format(args))
