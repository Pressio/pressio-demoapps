
import re, math
import numpy as np
import matplotlib.pyplot as plt

def extractN(ns,meshPath):
  reg = re.compile(r''+ns+'.+')
  file1 = open(meshPath+'/info.dat', 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return int(strings.group().split()[1])

state = np.loadtxt("state_final.txt")
state_x = np.loadtxt("coordinates.dat")[:,1]
state_y = np.loadtxt("coordinates.dat")[:,2]

D = np.loadtxt("grad_result_final.txt")
mask_gx = (D[:,6] == 1)
mask_gy = (D[:,6] == 2)
Xpt_gx = D[mask_gx,0]
Ypt_gx = D[mask_gx,1]
Xpt_gy = D[mask_gy,0]
Ypt_gy = D[mask_gy,1]

for dof in range(2,6):
  comp_gx_V = D[mask_gx,2]
  comp_gy_V = D[mask_gy,2]

  ax1 = plt.figure(dof).add_subplot(projection='3d')
  ax1.view_init(21, -32)
  ax1.scatter3D(state_x, state_y, state[dof-2::4], marker='*', alpha=0.4, label="field")
  ax1.scatter3D(Xpt_gx, Ypt_gx, comp_gx_V, alpha=0.9, label="grad_x")
  ax1.scatter3D(Xpt_gy, Ypt_gy, comp_gy_V, alpha=0.9, label="grad_y")
  ax1.legend()

plt.show()
