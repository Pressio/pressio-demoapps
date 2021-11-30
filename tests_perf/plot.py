#!/usr/bin/env python

from random import randrange, uniform
import re, sys, os, time, yaml, copy
import numpy as np
from argparse import ArgumentParser
import shutil, subprocess
import matplotlib.pyplot as plt

def extractTimes(d, order, N):
  rows=np.where(d[:,0]==order)
  d2 = d[rows]
  rows = np.where(d2[:,1]==N)
  return d2[rows][0,3]/d2[rows][:,3]

def orderToColor(order):
  if order==1: return 'k'
  elif order==3: return 'r'
  elif order==5: return 'g'

def cellsToMarker(N):
  if N==128: return 'o'
  elif N==256: return '^'
  elif N==512: return 'd'
  elif N==1024: return 's'

def orderToString(order):
  if order==1:   return 'First'
  elif order==3: return 'Weno3'
  elif order==5: return 'Weno5'

if __name__== "__main__":
  '''
  note that this script relies on data.dat
  to be in specific order.
  '''

  d = np.loadtxt("data.dat")
  print(d.shape)

  orders  = np.unique(d[:,0])
  cells   = np.unique(d[:,1])
  threads = np.unique(d[:,2])

  fig = plt.figure(0)
  for i in orders:
    for N in cells:
      G = extractTimes(d, i, N)
      ms  = cellsToMarker(N)
      col = orderToColor(i)
      plt.plot(threads, G, '-', color=col, \
               marker=ms, markersize=4, markerfacecolor='none', \
               linewidth=0.5, \
               label=orderToString(i)+"_N="+str(int(N)))

  plt.plot(threads, threads, '--k', label="Ideal scaling")
  plt.xlabel("Num of threads", fontsize=20)
  plt.ylabel("Speedup", fontsize=20)

  plt.legend(ncol=1, fontsize=8)
  fig.savefig("scaling.pdf", format="pdf", bbox_inches='tight', dpi=450)
  plt.show()
