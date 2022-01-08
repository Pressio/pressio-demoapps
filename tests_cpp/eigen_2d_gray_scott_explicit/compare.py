
import numpy as np
import sys, os, re

def extractN(ns):
  reg = re.compile(r''+ns+'.+')
  file1 = open('info.dat', 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return int(strings.group().split()[1])

if __name__== "__main__":
  nx = extractN('nx')
  ny = extractN('ny')
  totDofs = nx*ny*2

  D = np.fromfile("gs_2d_solution.bin")
  nt = int(np.size(D)/totDofs)
  print(nt)
  D = np.reshape(D, (nt, totDofs))
  y = D[-1, :]
  np.savetxt("field.txt", y)

  goldD = np.loadtxt("gold.txt")
  assert(np.allclose(y.shape, goldD.shape))
  assert(np.isnan(y).all() == False)
  assert(np.isnan(goldD).all() == False)
  assert(np.allclose(y, goldD)) #,rtol=1e-10, atol=1e-12))
