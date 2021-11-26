#!/usr/bin/env python

from random import randrange, uniform
import re, sys, os, time, yaml, copy
import numpy as np
from argparse import ArgumentParser
import shutil, subprocess

def calculateStencilSize(scheme):
  if scheme == "FirstOrder":
    return 3
  elif scheme == "Weno3":
    return 5
  elif scheme == "Weno5":
    return 7

def generateMesh(meshExeDir, outDir, N, s, prob):
  meshExe = 'create_full_mesh_for.py'
  owd = os.getcwd()
  os.chdir(meshExeDir)
  args = ("python", meshExe,
          "--problem", prob,
          "-n", str(N), str(N),
          "--outDir", outDir)
  print(args)
  popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()
  output = popen.stdout.read(); #print("\n", output)
  os.chdir(owd)

#=======================================================
def runExe(meshDir, loopCount, scheme, nThreads):
  my_env = os.environ.copy()
  my_env["OMP_NUM_THREADS"] = str(nThreads)
  #my_env["OMP_PLACES"]="{0}"
  my_env["OMP_PROC_BIND"]="true"

  args = ("./perf_exe",
          "-m", meshDir,
          "-n", str(loopCount),
          "-s", scheme
          )

  # launch subprocess
  print("Running: {} ".format(args))
  logfile = open('out.log', 'w')
  p = subprocess.Popen(args, stdout=logfile, stderr=logfile, env=my_env)
  p.wait()
  logfile.close()
  assert(p.returncode==0)
  return extractTime("out.log")

def extractTime(logFilePath):
  reg = re.compile(r'elapsed(\D+)\d+.\d+')
  file1 = open(logFilePath, 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return np.float64(strings.group().split()[1])

def schemeToInteger(scheme):
  if scheme == "FirstOrder":
    return 1
  elif scheme == "Weno3":
    return 3
  elif scheme == "Weno5":
    return 5

###############################
if __name__== "__main__":
###############################
  parser = ArgumentParser()
  parser.add_argument("-repoDir", "--repoDir",
                      dest="repoDir", default="empty")
  args = parser.parse_args()
  assert(args.repoDir != "empty")

  meshExeDir = args.repoDir + '/meshing_scripts'

  numTrials = 10 # how many times to run loop in exe
  meshList = [128, 256, 512, 1024]
  threads  = [1, 2, 4, 8]
  schemes = ["FirstOrder", "Weno3", "Weno5"]
  problem = "euler2dsmooth"

  owd = os.getcwd()
  data = []
  for iSch in schemes:
    s = calculateStencilSize(iSch)
    print(iSch, s)

    for N in meshList:
      meshDir = owd + "/mesh_" + problem + "_s" + str(s) + "_" + str(N)
      # check if dir exists, if not generate
      if not os.path.exists(meshDir):
        generateMesh(meshExeDir, meshDir, N, s, problem+"_s"+str(s))

      for numThreads in threads:
        currTime = runExe(meshDir, numTrials, iSch, numThreads)
        data.append([schemeToInteger(iSch), N, numThreads, currTime])

  print(data)
  np.savetxt("data.dat", data,
             fmt=['%d', '%d', '%d', '%.8f'])
