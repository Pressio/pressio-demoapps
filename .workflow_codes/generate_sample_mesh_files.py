
import os, random, sys, subprocess

def _generate_sample_mesh_files(fullMeshPath, outDir, smIndicesFile):
  meshScriptsDir = os.path.dirname(__file__) + "/../meshing_scripts"

  print('Generating sample mesh in:')
  print(' {}'.format(outDir))
  owd = os.getcwd()
  args = ("python3", meshScriptsDir+'/create_sample_mesh.py',
          "--fullMeshDir", fullMeshPath,
          "--sampleMeshIndices", smIndicesFile,
          "--outDir", outDir)
  popen  = subprocess.Popen(args, stdout=subprocess.PIPE); popen.wait()
  output = popen.stdout.read();
  print(output)
