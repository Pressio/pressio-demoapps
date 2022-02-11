
#==============================================================
# imports
#==============================================================

import os, subprocess

#==============================================================
# functions
#==============================================================

def execute_foms_custom_args(parser):
  parser.add_argument("--exedir",
                      dest="exeDir",
                      required=True,
                      help="Full path with executables. Required!")

  parser.add_argument("--mesh",
                      nargs=2, \
                      dest="meshSize", \
                      type=int, \
                      required=True, \
                      help="Meshing")

  parser.add_argument("--stencilsize",
                      dest="stencilSize", \
                      type=int, \
                      required=True)

  parser.add_argument("--nthreads", "--nth",
                      dest="numThreads", \
                      type=int, \
                      default=1,
                      help="Num threads to use, default = 1")

def run_fom_exe(runDirectory, args):
  # from parent module common to all demoapps
  from run_exe import _run_exe
  _run_exe(runDirectory, args)

def generate_or_set_full_mesh(workingDir, args):
  pdaMeshDir = os.path.dirname(__file__) + "/../../meshing_scripts"

  nx, ny = args.meshSize[0], args.meshSize[1]
  meshDir = workingDir + "/full_mesh" + str(nx) + "x" + str(ny)

  if os.path.exists(meshDir):
    # if mesh exists, do nothing
    print('Mesh {} already exists'.format(meshDir))
  else:
    # generate
    print('Generating mesh {}'.format(meshDir))
    # call script
    owd = os.getcwd()
    args = ("python3", pdaMeshDir + '/create_full_mesh.py',\
            "-n", str(nx), str(ny),\
            "--outDir", meshDir,\
            "--bounds", "0.0", "0.25", "0.0", "1.0",
            "-s", str(args.stencilSize))

    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()
    output = popen.stdout.read();

  return meshDir

def create_fom_train_input_dics(scenario, meshPath):
  from .impl import _create_fom_input_dics_impl
  return _create_fom_input_dics_impl(scenario, meshPath, "train")

def create_fom_test_input_dics(scenario, meshPath):
  from .impl import _create_fom_input_dics_impl
  return _create_fom_input_dics_impl(scenario, meshPath, "test")
