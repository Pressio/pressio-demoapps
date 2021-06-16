
import sys, os, time, subprocess
import numpy as np
from argparse import ArgumentParser

valid = ["linadv1d_s3", \
         "linadv1d_s7", \
         "euler1dsmooth_s3", \
         "euler1dsmooth_s7", \
         "sod1d_s3", \
         "sod1d_s7", \
         "lax1d_s3", \
         "lax1d_s7", \
         #
         "euler2dsmooth_s3", \
         "euler2dsmooth_s7", \
         "sedov2d_s3", \
         "sedov2d_s7", \
         "sedov2dsym_s3", \
         "sedov2dsym_s7", \
         "riemann2d_s3", \
         "riemann2d_s7", \
         #
         "euler3dsmooth_s3", \
         "sedov3dsym_s3"]

#-------------------------------
if __name__== "__main__":
#-------------------------------
  parser = ArgumentParser()
  parser.add_argument(
    "--name", "--problem",
    dest="problemName",
    help="choices:" + str(valid))

  parser.add_argument(
    "--outDir", "--outdir", dest="outDir",
    help="Full path to where to store all the mesh output files.")

  parser.add_argument(
    "-n", "--numCells",
    nargs="*", type=int, dest="numCells",
    help="Num of cells along x and y. \
    If you only pass one value, I assume 1d.\
    If you pass two values, I assume 2d.")

  args = parser.parse_args()
  pName = args.problemName
  print(pName)

  fullMeshScript = "create_full_mesh.py"

  if ("linadv1d" in pName):
    s = int(pName[-1])
    args = ("python", fullMeshScript,
            "-n", str(args.numCells[0]), str(1),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-1.0), str(1.0),
            "--periodic", "true")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("euler1dsmooth" in pName):
    s = int(pName[-1])
    args = ("python", fullMeshScript,
            "-n", str(args.numCells[0]), str(1),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-1.0), str(1.0),
            "--periodic", "true")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("sod1d" in pName):
    s = int(pName[-1])
    args = ("python", fullMeshScript,
            "-n", str(args.numCells[0]), str(1),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-0.5), str(0.5),
            "--periodic", "false")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("lax1d" in pName):
    s = int(pName[-1])
    args = ("python", fullMeshScript,
            "-n", str(args.numCells[0]), str(1),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-5.0), str(5.0),
            "--periodic", "false")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("euler2dsmooth" in pName):
    s = int(pName[-1])
    args = ("python", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-1.), str(1.0), str(-1.0), str(1.0),
            "--periodic", "true")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("sedov2d" in pName):
    s = int(pName[-1])
    args = ("python", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-0.6), str(0.6), str(-0.6), str(0.6),
            "--periodic", "false")
            # "-p", "print", str(0))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("sedov2dsym" in pName):
    s = int(pName[-1])
    args = ("python", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(0.0), str(1.2), str(0.0), str(1.2),
            "--periodic", "false")
            #"-p", "print", str(4))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("riemann2d" in pName):
    s = int(pName[-1])
    args = ("python", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(0.0), str(1.0), str(0.0), str(1.0),
            "--periodic", "false")
            # "-p", "print", str(0))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("euler3dsmooth" in pName):
    s = int(pName[-1])
    if (s!=3):
      print("3D currently only supports 1st order, use _s3")
      sys.exit(1)
    args = ("python", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]), str(args.numCells[2]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-1.0), str(1.0), str(-1.0), str(1.0), str(-1.0), str(1.0),
            "--periodic", "true")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("sedov3dsym" in pName):
    s = int(pName[-1])
    if (s!=3):
      print("3D currently only supports 1st order, use _s3")
      sys.exit(1)
    args = ("python", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]), str(args.numCells[2]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(0.0), str(1.2), str(0.0), str(1.2), str(0.0), str(1.2),
            "--periodic", "false")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  else:
    print("Invalid arg passed to --name, choices are:", valid)
    sys.exit()
