
#==============================================================
# imports
#==============================================================

# LOCAL imports
from .constants import numDofsPerCell

#==============================================================
# functions
#==============================================================

def execute_roms_custom_args(parser):
  parser.add_argument("--exedir",
                      dest="exeDir",
                      required=True,
                      help="Full path with executables. Required!")

  parser.add_argument("--nthreads", "--nth",
                      dest="numThreads", \
                      type=int, \
                      default=1,
                      help="Num threads to use, default = 1")

def create_rom_input_dics(scenario, meshPath, algoName, \
                          modesList, lsvFile, refStateFile, \
                          projectorFile = None,
                          # we need an optional sampleMeshPath because
                          # for masked galerkin the mesh to use is the FULL mesh
                          # but we still need to use gids from the sample mesh
                          sampleMeshPath = None):

  from .impl import _create_rom_input_dics_impl
  return _create_rom_input_dics_impl(scenario, meshPath, \
                                     algoName, modesList, \
                                     lsvFile, refStateFile,\
                                     projectorFile, sampleMeshPath)


def run_galerkin_rom_exe(runDirectory, args):
  # import from parent module common to all demoapps
  from run_exe import _run_exe
  _run_exe(runDirectory, args)


def compute_initial_rom_state(fomTestDir, romDir, dicIn, statePodDir, K):
  # import from parent module common to all demoapps
  from compute_initial_rom_state import _compute_initial_rom_state_impl
  _compute_initial_rom_state_impl(fomTestDir, romDir, dicIn, statePodDir, K)


def make_gappy_projector_dont_split_vars(outputFile, \
                                         sampleMeshPath, \
                                         statePodDir, \
                                         rhsPodDir, \
                                         K):
  # import from parent module common to all demoapps
  from make_projectors import _make_gappy_projector_dont_split_vars_impl
  _make_gappy_projector_dont_split_vars_impl(outputFile, sampleMeshPath, \
                                             statePodDir, rhsPodDir, \
                                             K, numDofsPerCell)
