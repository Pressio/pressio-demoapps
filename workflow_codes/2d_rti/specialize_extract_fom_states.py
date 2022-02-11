
#==============================================================
# imports
#==============================================================

# from parent module common to all demoapps
from extract_fom_states import \
  _extract_initial_and_final_state_from_run_directory_impl

# specific to this problem
from .constants import numDofsPerCell

#==============================================================
# functions
#==============================================================

def extract_initial_and_final_state_from_run_directory(fomRunDir):
  _extract_initial_and_final_state_from_run_directory_impl(fomRunDir, \
                                                           numDofsPerCell)
