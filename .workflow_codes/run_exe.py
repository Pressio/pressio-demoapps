
import sys, os, subprocess
import numpy as np

def _link_exe(whereToRun, exeDir, exeName):
  # check if exeDir exists and contains the executables
  if not os.path.exists(exeDir):
    print("exeDir {} does not exist".format(exeDir))
    sys.exit(1)

  # check if exe exist
  sourceExe  = exeDir+"/"+exeName
  if not os.path.exists(sourceExe):
    print("exe {} does not exist".format(sourceExe))
    sys.exit(1)

  # delete if already in target dir, and relink
  destExe  = whereToRun+"/"+exeName
  if os.path.exists(destExe):
    os.system("rm -f {}".format(destExe))

  os.system("ln -s {} {}".format(sourceExe, destExe))


def _run_exe(whereToRun, args):
  exeName = "pdaExeWf"
  _link_exe(whereToRun, args.exeDir, exeName)

  my_env = os.environ.copy()
  if (args.numThreads ==1):
    my_env["OMP_NUM_THREADS"] = str(1)
    my_env["OMP_PLACES"]="{0}"
    my_env["OMP_PROC_BIND"]="true"

  else:
    my_env["OMP_NUM_THREADS"] = str(args.numThreads)
    my_env["OMP_PLACES"]      ="threads"
    my_env["OMP_PROC_BIND"]   ="spread"

  args = ("./"+exeName, whereToRun+"/input.yaml")

  # launch subprocess
  print("Running: {}".format(args))
  print("     in: {}".format(whereToRun))
  logfile = open(whereToRun+'/out.log', 'w')
  p = subprocess.Popen(args, stdout=logfile, stderr=logfile, env=my_env, cwd=whereToRun)
  p.wait()
  logfile.close()
  assert(p.returncode==0)
