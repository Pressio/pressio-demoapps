"""Root module of your package"""

from pkg_resources import get_distribution

__version__ = get_distribution('pressiodemoapps').version


'''
see this for why this file exists and is done this way
https://stackoverflow.com/questions/47599162/pybind11-how-to-package-c-and-python-code-into-a-single-package?rq=1
'''
try:
  from ._pressiodemoappsimpl import CellCenteredUniformMesh, loadCellCenterUniformMesh
except ImportError:
  raise ImportError("Unable to import mesh classes from _pressiodemoappsimpl")

try:
  from ._pressiodemoappsimpl \
    import InviscidFluxReconstruction, InviscidFluxScheme, \
    ViscousFluxReconstruction, ViscousFluxScheme, \
    Advection1d, DiffusionReaction1d, Euler1d, Euler2d, Euler3d, Swe2d
except ImportError:
  raise ImportError("Unable to import enums from _pressiodemoappsimpl")

try:
  from ._pressiodemoappsimpl import Advection1dProblem
except ImportError:
  raise ImportError("Unable to import Advection classes from _pressiodemoappsimpl")

try:
  from ._pressiodemoappsimpl import DiffusionReaction1dProblem
except ImportError:
  raise ImportError("Unable to import DiffusionReaction classes from _pressiodemoappsimpl")

try:
  from ._pressiodemoappsimpl import Euler1dProblem, Euler2dProblem, Euler3dProblem
except ImportError:
  raise ImportError("Unable to import Euler classes from _pressiodemoappsimpl")

  from ._pressiodemoappsimpl import Swe2dProblem
except ImportError:
  raise ImportError("Unable to import SWE classes from _pressiodemoappsimpl")

try:
  from ._pressiodemoappsimpl import createProblem
except ImportError:
  raise ImportError("Unable to import createProblem from _pressiodemoappsimpl")


# Runge-Kutta4 integrator
def advanceRK4(appObj, state, dt, Nsteps, \
               startTime = 0.0, \
               observer = None, \
               showProgress=False):

  v = appObj.createVelocity()
  tmpState = state.copy()
  half = 0.5
  two  = 2.
  oneOverSix = 1./6.

  time = startTime
  for step in range(1, Nsteps+1):
    if showProgress:
      if step % 50 == 0: print("step = ", step, "/", Nsteps)

    appObj.velocity(state, time, v)
    if observer!= None:
      observer(step-1, state, v)
    k1 = dt * v

    tmpState = state+half*k1
    appObj.velocity(tmpState, time+half*dt, v)
    k2 = dt * v

    tmpState = state+half*k2
    appObj.velocity(tmpState, time+half*dt, v)
    k3 = dt * v

    tmpState = state+k3
    appObj.velocity(tmpState, time, v)
    k4 = dt * v

    state[:] = state + (k1+two*k2+two*k3+k4)*oneOverSix
    time += dt


# SSP3 integrator
def advanceSSP3(appObj, state, dt, Nsteps,\
                startTime = 0.0, \
                observer = None, \
                showProgress=False):

  v = appObj.createVelocity()
  tmpState1 = state.copy()
  tmpState2 = state.copy()

  oneOverThree  = 1./3.
  oneOverFour   = 1./4.
  threeOverFour = 3./4.
  two = 2.

  time = startTime
  for step in range(1, Nsteps+1):
    if showProgress:
      if step % 50 == 0: print("step = ", step, "/", Nsteps)

    appObj.velocity(state, time, v)
    if observer!= None:
      observer(step-1, state, v)
    tmpState1[:] = state + dt * v

    appObj.velocity(tmpState1, time, v)
    tmpState2[:] = threeOverFour*state + oneOverFour*tmpState1 + oneOverFour*dt*v

    appObj.velocity(tmpState2, time, v)
    state[:] = oneOverThree*(state + two*tmpState2 + two*dt*v)

    time += dt
