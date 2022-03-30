#!/usr/bin/env python

import numpy as np
import sys, os
from argparse import ArgumentParser
from numpy import linalg as la

def str2bool(v):
  if isinstance(v, bool):
    return v
  if v.lower() in ('yes', 'true', 't', 'y', '1'):
    return True
  elif v.lower() in ('no', 'false', 'f', 'n', '0'):
    return False
  else:
    raise argparse.ArgumentTypeError('Boolean value expected.')

#======================================================================
def compute_errors_between_fom_and_approximated_state(fomState, fomStateApprox):
  fomL2Norm   = la.norm(fomState)
  fomLinfNorm = la.norm(fomState, np.inf)
  errorField  = fomState-fomStateApprox

  storage = []
  storage.append([ np.min(fomState),              np.max(fomState) ])
  storage.append([ np.min(fomStateApprox),        np.max(fomStateApprox) ])
  storage.append([ np.min(errorField),            np.max(errorField) ])
  storage.append([ la.norm(errorField),           la.norm(errorField)/fomL2Norm ])
  storage.append([ la.norm(errorField, np.inf),   la.norm(errorField, np.inf)/fomLinfNorm ])
  return np.array(storage)

#======================================================================
def loadStates(fomFile, approxFile, numRowsToSkipWhenReading):
  # load states
  if not os.path.exists(fomFile):
    print("fom final state {} does not exist".format(fomFile))
    sys.exit(1)

  if not os.path.exists(approxFile):
    print("approx final state {} does not exist".format(approxFile))
    sys.exit(1)

  print("reading fom    state {}".format(fomFile))
  print("reading approx state {}".format(approxFile))

  fomState       = np.loadtxt(fomFile,    skiprows=numRowsToSkipWhenReading)
  fomStateApprox = np.loadtxt(approxFile, skiprows=numRowsToSkipWhenReading)
  fomNRows       = fomState.shape[0]
  approxNRows    = fomStateApprox.shape[0]
  # if things are correct, this should be a single vector and approx/fom match sizes
  assert( fomNRows == approxNRows )

  return [fomState, fomStateApprox]

###############################
if __name__== "__main__":
###############################
  parser = ArgumentParser()
  parser.add_argument("--dryrun", "--dr",
                      dest="dryRun", type=str2bool, default=True,
                      help="True: creates directory structures/files, does not run. Default=True.")

  parser.add_argument("--fom-state", "--fomstate",
                      dest="fomState", default="empty",
                      help="Full path to fom state. Must be set.")

  parser.add_argument("--approx-state", "--approxstate",
                      dest="approxState", default="empty",
                      help="Full path to approx state. Must be set.")

  parser.add_argument("--dof-name", "--dofname",
                      dest="dofName", default="empty")

  parser.add_argument("--skiprows",
                      dest="skipRows", type=int, default=0,
                      help="If to read data we need to skip rows")

  # parse args
  args = parser.parse_args()
  assert(args.fomState != "empty")
  assert(args.approxState != "empty")
  dofName = args.dofName

  [fomState, fomStateApprox] = loadStates(args.fomState, args.approxState, args.skipRows)

  M = compute_errors_between_fom_and_approximated_state(fomState, fomStateApprox)

  #error = fomState-fomStateApprox
  #fomL2Norm, fomLinfNorm = la.norm(fomState), la.norm(fomState, np.inf)
  #errL2Norm   = [ la.norm(error),         la.norm(error)/fomL2Norm ]
  #errLinfNorm = [ la.norm(error, np.inf), la.norm(error, np.inf)/fomLinfNorm ]

  if dofName!="empty":
    print(" {}_fom_minmax    = {} {}".format(dofName, M[0,0], M[0,1]))
    print(" {}_approx_minmax = {} {}".format(dofName, M[1,0], M[1,1]))
    print(" {}_err_minmax    = {} {}".format(dofName, M[2,0], M[2,1]))
    print(" {}_err_abs_rel_ltwo_norms = {} {}".format(dofName, M[3,0], M[3,1]))
    print(" {}_err_abs_rel_linf_norms = {} {}".format(dofName, M[4,0], M[4,1]))
  else:
    print(" fom_minmax    = {} {}".format(M[0,0], M[0,1]))
    print(" approx_minmax = {} {}".format(M[1,0], M[1,1]))
    print(" err_minmax    = {} {}".format(M[2,0], M[2,1]))
    print(" err_abs_rel_ltwo_norms = {} {}".format(M[3,0], M[3,1]))
    print(" err_abs_rel_linf_norms = {} {}".format(M[4,0], M[4,1]))
