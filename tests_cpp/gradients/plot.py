
import re, math
import numpy as np
import matplotlib.pyplot as plt

def fun(x,y):
    return np.sin(math.pi * x * y)

def gradX(x,y):
    return math.pi*y*np.cos(math.pi * x * y)

def gradY(x,y):
    return math.pi*x*np.cos(math.pi * x * y)

def extractN(ns,meshPath):
  reg = re.compile(r''+ns+'.+')
  file1 = open(meshPath+'/info.dat', 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return int(strings.group().split()[1])

D = np.loadtxt("pda_resultsample_mesh_ss_7.txt")
mask_gx = (D[:,4] == 1)
mask_gy = (D[:,4] == 2)
Xpt_gx = D[mask_gx,0]
Ypt_gx = D[mask_gx,1]
Xpt_gy = D[mask_gy,0]
Ypt_gy = D[mask_gy,1]
comp_gx_V = D[mask_gx,2]
comp_gy_V = D[mask_gy,2]
gold_gx_V = D[mask_gx,3]
gold_gy_V = D[mask_gy,3]



meshPath = "./fullmesh_s3"
nx = extractN('nx',meshPath)
ny = extractN('ny',meshPath)
print(nx, ny)

x = np.linspace(0, 1., nx)
y = np.linspace(0, 1., ny)
xv, yv = np.meshgrid(x, y)

f = fun(xv,yv)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.view_init(21, -32)
surf = ax.plot_surface(xv, yv, f, linewidth=0, antialiased=False)
ax.set_title("function")

gx = gradX(xv,yv)
fig, ax1 = plt.subplots(subplot_kw={"projection": "3d"})
ax1.view_init(21, -32)
ax1.plot_surface(xv, yv, gx, linewidth=0, antialiased=False, alpha=0.5)
ax1.set_title("grad_X")

gx = gradY(xv,yv)
fig, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
ax2.view_init(21, -32)
ax2.plot_surface(xv, yv, gx, linewidth=0, antialiased=False, alpha=0.5)
ax2.set_title("grad_Y")

ax1.scatter3D(Xpt_gx, Ypt_gx, comp_gx_V, alpha=0.9)
ax2.scatter3D(Xpt_gy, Ypt_gy, comp_gy_V, alpha=0.9)
ax1.legend()
ax2.legend()

fig, ax3 = plt.subplots()
fig, ax4 = plt.subplots()
ax3.plot(gold_gx_V,'o', label="gold")
ax3.plot(comp_gx_V,'*', label="computed")
ax4.plot(gold_gy_V,'o', label="gold")
ax4.plot(comp_gy_V,'*', label="computed")
ax3.set_title("grad_X")
ax4.set_title("grad_Y")
ax3.legend()
ax4.legend()

plt.show()
