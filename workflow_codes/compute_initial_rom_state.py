
#==============================================================
# imports
#==============================================================

import os
import numpy as np
from pod_functions import \
  _load_basis_from_binary_file

#==============================================================
# functions
#==============================================================

def _compute_initial_rom_state_impl(fomTestDir, romDir, \
                                    dicIn, statePodDir, K):

  # set POD, refstate location
  lsvFile      = statePodDir + '/lsv'
  refStateFile = statePodDir + '/referenceState'
  fullPhi      = _load_basis_from_binary_file(lsvFile)
  refState     = np.loadtxt(refStateFile)
  assert(fullPhi.shape[0] == refState.shape[0])

  fomIC = np.loadtxt(fomTestDir+"/initial_state.txt")
  diff  = fomIC - refState
  phi   = fullPhi[:, 0:K]
  romIC = np.dot(phi.transpose(), diff)
  romICFile = romDir + "/romInitialCond"
  np.savetxt(romICFile, romIC)

  dicIn['rom']['romInitialStateFile'] = romICFile
