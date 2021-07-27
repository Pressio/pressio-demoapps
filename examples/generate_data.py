
import pathlib, sys
file_path = pathlib.Path(__file__).parent.absolute()

import numpy as np
import time
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import cm
import pressiodemoapps as pda
from argparse import ArgumentParser

def computePressure(rho, u, v, E):
  gamma_ = (5.+2.)/5.
  vel = u**2 + v**2
  return (gamma_ - 1.) * (E - rho*vel*0.5)

class Observer:
  def __init__(self, saveSnapshotsEvery):
    self.stateSnaps = []
    self.rhsSnaps   = []
    self.saveSnapshotsEvery = saveSnapshotsEvery

  def __call__(self, timeStep, state, rhs):
    if timeStep % self.saveSnapshotsEvery == 0:
      self.stateSnaps.append(state.tolist())
      self.rhsSnaps.append(rhs.tolist())
      print("Collect at step = ", timeStep)

  def saveToFile(self):
    np.savetxt("state_snapshots.txt", self.stateSnaps)
    np.savetxt("rhs_snapshots.txt",   self.rhsSnaps)

# -------------------------------
if __name__ == '__main__':
# -------------------------------
  parser = ArgumentParser()
  parser.add_argument(
    "-s", "--sampleFrequency",
    dest="sampleFreq", type=int,
    help="Frequency for sampling state and rhs.")

  parser.add_argument(
    "-T", "--finalTime",
    dest="finalTime", type=float,
    help="Simulation time.")

  parser.add_argument(
    "--dt", "--timeStepSize",
    dest="dt", type=float,
    help="Time step size.")

  args = parser.parse_args()

  meshPath = str(file_path)
  meshObj  = pda.loadCellCenterUniformMesh(meshPath)
  fluxEnum = pda.InviscidFluxReconstruction.Weno3
  probName = pda.Euler2d.DoubleMachReflection
  appObj   = pda.createProblem(meshObj, probName, fluxEnum)
  state    = appObj.initialCondition()

  dt = args.dt
  Nsteps = int(float(args.finalTime)/dt)
  collector = Observer(int(args.sampleFreq))
  pda.advanceSSP3(appObj, state, dt, Nsteps, observer=collector, showProgress=True)
  collector.saveToFile()
