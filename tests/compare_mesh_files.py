
import numpy as np
import sys, os

if __name__== "__main__":
  print('Argument List:', str(sys.argv))
  tol   = float(sys.argv[1])

  files = ["info.dat", "connectivity.dat", "coordinates.dat"]
  for f in files:
    if not os.path.isfile(f):
      sys.exit(1)

  # info file
  file1 = open('info_gold.dat', 'r')
  file2 = open('info.dat', 'r')
  while True:
    line1 = file1.readline()
    line2 = file2.readline()
    assert(line1==line2)
    if not line1 or not line2:
      break
  file1.close()  
  file2.close()  

  # connectivity
  testD = np.loadtxt("connectivity.dat")
  goldD = np.loadtxt("connectivity_gold.dat")
  assert(np.allclose(testD.shape, goldD.shape))
  assert(np.isnan(testD).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(testD, goldD))

  # coordinates
  testD = np.loadtxt("coordinates.dat")
  goldD = np.loadtxt("coordinates_gold.dat")
  assert(np.allclose(testD.shape, goldD.shape))
  assert(np.isnan(testD).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(testD[:,0], goldD[:,0]))
  assert(np.allclose(testD[:,1:-1], goldD[:,1:-1], atol=tol))
