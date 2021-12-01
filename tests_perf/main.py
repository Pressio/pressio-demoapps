#!/usr/bin/env python

import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

from argparse import ArgumentParser
import numpy as np
import time
import pressiodemoapps as pda

def schemeStringToSchemeEnum(s):
  if s=="FirstOrder":
    return pda.InviscidFluxReconstruction.FirstOrder
  elif s=="Weno3":
    return pda.InviscidFluxReconstruction.Weno3
  elif s=="Weno5":
    return pda.InviscidFluxReconstruction.Weno5

# ---------------------------
if __name__ == '__main__':
# ---------------------------
  parser = ArgumentParser()
  parser.add_argument("-m, --mesh", dest="meshDir", default="empty")
  parser.add_argument("-n", dest="loopCount", default=10, type=int)
  parser.add_argument("-s", dest="scheme");
  args = parser.parse_args()

  start = time.time()
  meshPath = str(args.meshDir)
  meshObj  = pda.load_cellcentered_uniform_mesh(meshPath)
  schemeEnum = schemeStringToSchemeEnum(args.scheme)
  probId     = pda.Euler2d.PeriodicSmooth
  appObj     = pda.create_problem(meshObj, probId, schemeEnum)

  yn = appObj.initialCondition()
  V = appObj.createVelocity()
  B = np.ones((len(yn), 25), order='F')
  AJ = appObj.createApplyJacobianResult(B)
  # warmup
  appObj.applyJacobian(yn, 0., B, AJ)
  #appObj.velocity(yn, 0., V)

  print("starting loop")
  start = time.time()
  for i in range(args.loopCount):
    #appObj.velocity(yn, 0., V)
    appObj.applyJacobian(yn, 0., B, AJ)

  end = time.time()
  print("elapsed ", end - start)
