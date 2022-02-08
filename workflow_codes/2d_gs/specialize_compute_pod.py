
#==============================================================
# imports
#==============================================================

from .constants import numDofsPerCell

#==============================================================
# functions (for when we don't split variables)
#==============================================================

def compute_state_pod_and_ref_state_dont_split_vars(dryRun, \
                                                    dataDirs,\
                                                    singValsFile,\
                                                    leftSingVecsFile, \
                                                    referenceStateFile):

  from pod_functions import _compute_state_pod_and_ref_state_dont_split_vars_impl
  _compute_state_pod_and_ref_state_dont_split_vars_impl(dryRun,\
                                                        numDofsPerCell,\
                                                        dataDirs,\
                                                        singValsFile,\
                                                        leftSingVecsFile, \
                                                        referenceStateFile)

def compute_rhs_pod_if_needed_dont_split_vars(dryRun, \
                                              dataDirs,\
                                              singValsFile,\
                                              leftSingVecsFile):
  from pod_functions import _compute_rhs_pod_if_needed_dont_split_vars_impl
  _compute_rhs_pod_if_needed_dont_split_vars_impl(dryRun, \
                                                  numDofsPerCell,\
                                                  dataDirs,\
                                                  singValsFile,\
                                                  leftSingVecsFile)

#==============================================================
# functions (for when we split variables)
#==============================================================

def compute_state_pod_and_ref_state_split_vars(dryRun, \
                                               dataDirs, \
                                               singValsFile,\
                                               leftSingVecsFile):
  # import from parent module common to all demoapps
  from pod_functions import _compute_state_pod_and_ref_state_split_vars_impl
  _compute_state_pod_and_ref_state_split_vars_impl(dryRun, \
                                                   numDofsPerCell,\
                                                   dataDirs, \
                                                   singValsFile,\
                                                   leftSingVecsFile)

def compute_rhs_pod_if_needed_split_vars(dryRun, dataDirs, outDir):
  # import from parent module common to all demoapps
  from pod_functions import _compute_rhs_pod_if_needed_split_vars_impl
  _compute_rhs_pod_if_needed_split_vars_impl(dryRun, \
                                             numDofsPerCell,\
                                             dataDirs,\
                                             singValsFile, \
                                             leftSingVecsFile)
