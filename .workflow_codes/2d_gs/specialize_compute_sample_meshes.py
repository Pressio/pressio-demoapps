
#==============================================================
# imports
#==============================================================

# specific to this problem
from .constants import numDofsPerCell, dimensionality

#==============================================================
# functions
#==============================================================

def make_random_sample_mesh(dryRun, workDir, fraction, outDir):
  # import from parent module common to all demoapps
  from compute_sample_meshes import _make_random_sample_mesh_impl
  _make_random_sample_mesh_impl(dryRun, workDir, \
                                fraction, \
                                outDir, \
                                dimensionality)


def make_leverage_scores_sample_mesh(dryRun, workDir, \
                                     setId, dataDirs, \
                                     listOfFractions, \
                                     listOfBetas):
  # import from parent module common to all demoapps
  from compute_sample_meshes import _make_leverage_scores_sample_mesh_impl
  _make_leverage_scores_sample_mesh_impl(dryRun, workDir, \
                                         setId, dataDirs, \
                                         listOfFractions, \
                                         listOfBetas, \
                                         numDofsPerCell)
