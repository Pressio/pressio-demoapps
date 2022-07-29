
#==============================================================
# imports
#==============================================================

import numpy as np
from extract_from_run_directory import \
  find_total_dofs_from_run_directory

#==============================================================
# functions
#==============================================================

def _extract_initial_and_final_state_from_run_directory_impl(fomRunDir, \
                                                             numDofsPerCell):

  fomTotDofs = find_total_dofs_from_run_directory(fomRunDir, numDofsPerCell)
  fomSnaps   = np.fromfile(fomRunDir+"/fom_snaps_state.bin")
  nt         = int(np.size(fomSnaps)/fomTotDofs)
  fomData    = np.reshape(fomSnaps, (nt, fomTotDofs))
  np.savetxt(fomRunDir+"/initial_state.txt", fomData[0,:])
  np.savetxt(fomRunDir+"/final_state.txt",   fomData[-1,:])
