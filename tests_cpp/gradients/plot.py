
import re, math
import numpy as np
import matplotlib.pyplot as plt

def fun(x,y):
    return np.sin(math.pi * x * y)

def gradX(x,y):
    return math.pi*y*np.cos(math.pi * x * y)

def gradY(x,y):
    return math.pi*x*np.cos(math.pi * x * y)

def extractN(ns):
  reg = re.compile(r''+ns+'.+')
  file1 = open('info.dat', 'r')
  strings = re.search(reg, file1.read())
  file1.close()
  assert(strings)
  return int(strings.group().split()[1])

nx = extractN('nx')
ny = extractN('ny')
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


fig, ax3 = plt.subplots()
ax3.set_title("grad_X")
fig, ax4 = plt.subplots()
ax4.set_title("grad_Y")

mylabels = ['computed_ss3']#, 'computed_ss5']
for i in range(1,2):
  myresults = np.loadtxt("fd_results_"+str(i)+".txt")
  mask_gx = (myresults[:,4] == 1)
  mask_gy = (myresults[:,4] == 2)
  comp_gx_X = myresults[mask_gx,1]
  comp_gx_Y = myresults[mask_gx,2]
  comp_gx_V = myresults[mask_gx,3]
  comp_gy_X = myresults[mask_gy,1]
  comp_gy_Y = myresults[mask_gy,2]
  comp_gy_V = myresults[mask_gy,3]
  ax1.scatter3D(comp_gx_X, comp_gx_Y, comp_gx_V, label=mylabels[0], alpha=0.5, s=15)
  ax2.scatter3D(comp_gy_X, comp_gy_Y, comp_gy_V, label=mylabels[0], alpha=0.5, s=15)
  ax3.plot(comp_gx_V,'o', label=mylabels[0])
  ax4.plot(comp_gy_V,'o', label=mylabels[0])

mylabels = ['pda_computed']
myresults = np.loadtxt("pda_result.txt")
mask_gx = (myresults[:,4] == 1)
mask_gy = (myresults[:,4] == 2)
comp_gx_X = myresults[mask_gx,1]
comp_gx_Y = myresults[mask_gx,2]
comp_gx_V = myresults[mask_gx,3]
comp_gy_X = myresults[mask_gy,1]
comp_gy_Y = myresults[mask_gy,2]
comp_gy_V = myresults[mask_gy,3]
ax1.scatter3D(comp_gx_X, comp_gx_Y, comp_gx_V, label=mylabels[0], alpha=0.9)
ax2.scatter3D(comp_gy_X, comp_gy_Y, comp_gy_V, label=mylabels[0], alpha=0.9)
ax3.plot(comp_gx_V,'*', label=mylabels[0])
ax4.plot(comp_gy_V,'*', label=mylabels[0])

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()
plt.show()
