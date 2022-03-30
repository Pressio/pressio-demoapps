#!/usr/bin/env python

import numpy as np
import sys, os
from argparse import ArgumentParser

def compute_cumulative_energy(s, target):
  sSq = np.square(s)
  den = np.sum(sSq)
  rsum = 0.
  for i in range(0, len(s)):
    rsum += sSq[i]
    ratio = (rsum/den)
    if ratio >= target:
      return i
  return len(s)

if __name__== "__main__":
  parser = ArgumentParser()
  parser.add_argument("--singvals",
                      dest="svFile", default="empty",
                      help="Path to file with sing values.")
  parser.add_argument("--percent", "-p",
                      dest="pct", default=999, type=np.float,
                      help="Target fraction to keep")
  # parse all args
  args = parser.parse_args()
  assert(args.svFile != "empty")

  # convert percentage to decimal
  target = float(args.pct)/100.

  # load data
  sv = np.loadtxt(args.svFile)

  # compute cumulative energy
  if (target == 1.):
    n = len(sv)
    print ("Nv: {}".format(n) )
  else:
    n = compute_cumulative_energy(sv, target)

  print ("Nv: {}".format(n) )
  print ("numBasis: {}".format(n) )
