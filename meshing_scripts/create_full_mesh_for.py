
import sys, os, time, subprocess
import numpy as np
from argparse import ArgumentParser
from pathlib import Path

valid = ["linadv1d_s3", \
         "linadv1d_s5", \
         "linadv1d_s7", \
         #
         "diffreac1d", \
         #
         "euler1dsmooth_s3", \
         "euler1dsmooth_s5", \
         "euler1dsmooth_s7", \
         #
         "sod1d_s3", \
         "sod1d_s5", \
         "sod1d_s7", \
         #
         "lax1d_s3", \
         "lax1d_s5", \
         "lax1d_s7", \
         #
         "shuOsher1d_s3", \
         "shuOsher1d_s5", \
         "shuOsher1d_s7", \
         #
         ##############
         ##### 2d #####
         ##############
         #
         "diffreac2d", \
         #
         "grayscott2d", \
         #
         "euler2dsmooth_s3", \
         "euler2dsmooth_s5", \
         "euler2dsmooth_s7", \
         #
         "euler2dKelvinHelmholtz_s3", \
         "euler2dKelvinHelmholtz_s5", \
         "euler2dKelvinHelmholtz_s7", \
         #
         "normalshock2d_s3", \
         "normalshock2d_s5", \
         "normalshock2d_s7", \
         #
         "sedov2d_s3", \
         "sedov2d_s5", \
         "sedov2d_s7", \
         #
         "sedov2dsym_s3", \
         "sedov2dsym_s5", \
         "sedov2dsym_s7", \
         #
         "riemann2d_s3", \
         "riemann2d_s5", \
         "riemann2d_s7", \
         #
         "doublemach2d_s3", \
         "doublemach2d_s5", \
         #
         "swe2dSlipWall_s3",\
         "swe2dSlipWall_s5",\
         "swe2dSlipWall_s7",\
         "swe2dslipwall_s3",\
         "swe2dslipwall_s5",\
         "swe2dslipwall_s7",\
         #
         "burgers2d_periodic_s3", \
         "burgers2d_periodic_s5", \
         "burgers2d_periodic_s7", \
         #
         "burgers2d_dirichlet_s3", \
         "burgers2d_dirichlet_s5", \
         "burgers2d_dirichlet_s7", \
         #
         ##############
         ##### 3d #####
         ##############
         "euler3dsmooth_s3", \
         "euler3dsmooth_s5", \
         #
         "sedov3dsym_s3",
         "sedov3dsym_s5"]

#-------------------------------
if __name__== "__main__":
#-------------------------------
  parser = ArgumentParser()
  parser.add_argument(
    "--name", "--problem",
    dest="problemName",
    help="choices:" + str(valid))

  parser.add_argument(
    "--outDir", "--outdir",
    dest="outDir",
    help="Full path to where to store all the mesh output files.")

  parser.add_argument(
    "-n", "--numCells",
    nargs="*",
    type=int,
    dest="numCells",
    help="Num of cells along x and y. \
    If you only pass one value, I assume 1d.\
    If you pass two values, I assume 2d.")

  args = parser.parse_args()
  # -----------------------------

  pName = args.problemName
  print(pName)

  # we need to do this so that the file is always run properly
  fullMeshScript = str(Path(__file__).parent.absolute())+"/create_full_mesh.py"

  if ("linadv1d" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(1),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-1.0), str(1.0),
            "--periodic", "x")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("diffreac1d" in pName):
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(1),
            "--outDir", args.outDir,
            "--stencilSize", str(3),
            "--bounds", str(0.0), str(1.0))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("euler1dsmooth" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(1),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-1.0), str(1.0),
            "--periodic", "x")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("sod1d" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(1),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-0.5), str(0.5))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("shuOsher1d" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(1),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-5.0), str(5.0))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("lax1d" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(1),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-5.0), str(5.0))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("diffreac2d" in pName):
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(3),
            "--bounds", str(0.0), str(1.0), str(0.0), str(1.0))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("grayscott2d" in pName):
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(3),
            "--bounds", str(-1.25), str(1.25), str(-1.25), str(1.25),
            "--periodic", "x y")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("euler2dsmooth" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-1.), str(1.0), str(-1.0), str(1.0),
            "--periodic", "x y")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("euler2dKelvinHelmholtz" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-5.0), str(5.0), str(-5.0), str(5.0),
            "--periodic", "x y")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("normalshock2d" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(0.0), str(2.0), str(0.0), str(1.0))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("sedov2d" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-1.2), str(1.2), str(-1.2), str(1.2))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("sedov2dsym" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(0.0), str(1.2), str(0.0), str(1.2))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("riemann2d" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(0.0), str(1.0), str(0.0), str(1.0))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("doublemach2d" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(0.0), str(4.0), str(0.0), str(1.0))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("swe2dSlipWall" in pName or "swe2dslipwall" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-5.), str(5.), str(-5.), str(5.))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("euler3dsmooth" in pName):
    s = int(pName[-1])
    if (s!=3 or s!=5):
      print("3D currently only supports 1st order or Weno3, use _s3 or _s5")
      sys.exit(1)
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]), str(args.numCells[2]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-1.0), str(1.0), str(-1.0), str(1.0), str(-1.0), str(1.0),
            "--periodic", "x y z")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("sedov3dsym" in pName):
    s = int(pName[-1])
    if (s!=3 or s!=5):
      print("3D currently only supports 1st order or Weno3, use _s3 or _s5")
      sys.exit(1)
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]), str(args.numCells[2]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(0.0), str(1.2), str(0.0), str(1.2), str(0.0), str(1.2))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("burgers2d_periodic" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-1.0), str(1.0), str(-1.0), str(1.0),
            "--periodic", "x y")
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  elif ("burgers2d_dirichlet" in pName):
    s = int(pName[-1])
    args = ("python3", fullMeshScript,
            "-n", str(args.numCells[0]), str(args.numCells[1]),
            "--outDir", args.outDir,
            "--stencilSize", str(s),
            "--bounds", str(-1.0), str(1.0), str(-1.0), str(1.0))
    popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()

  else:
    print("Invalid arg passed to --name, choices are:", valid)
    sys.exit()
