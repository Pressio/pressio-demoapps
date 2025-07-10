"""Root module of your package"""

from pkg_resources import get_distribution

__version__ = get_distribution('pressiodemoapps').version


'''
see this for why this file exists and is done this way
https://stackoverflow.com/questions/47599162/pybind11-how-to-package-c-and-python-code-into-a-single-package?rq=1
'''

# meshing classes
try:
  from ._pressiodemoappsimpl import \
    CellCenteredUniformMesh, load_cellcentered_uniform_mesh
except ImportError:
  raise ImportError("Unable to import mesh classes from _pressiodemoappsimpl")

# all enums
try:
  from ._pressiodemoappsimpl import \
    InviscidFluxReconstruction, InviscidFluxScheme, \
    ViscousFluxReconstruction, ViscousFluxScheme, \
    Advection1d, DiffusionReaction1d, Euler1d, \
    AdvectionDiffusion2d, DiffusionReaction2d,\
    AdvectionDiffusionReaction2d, Swe2d, Euler2d, Euler3d
except ImportError:
  raise ImportError("Unable to import enums from _pressiodemoappsimpl")

# advection problems
try:
  from ._pressiodemoappsimpl \
    import Advection1dProblem
except ImportError:
  raise ImportError("Unable to import Advection classes from _pressiodemoappsimpl")

# advection-diffusion problems
try:
  from ._pressiodemoappsimpl \
    import AdvectionDiffusion2dProblem
except ImportError:
  raise ImportError("Unable to import AdvectionDiffusion classes from _pressiodemoappsimpl")

# diffusion-reaction problems
try:
  from ._pressiodemoappsimpl \
  import DiffusionReaction1dProblem, DiffusionReaction2dProblem
except ImportError:
  raise ImportError("Unable to import DiffusionReaction classes from _pressiodemoappsimpl")

# advection-diffusion-reaction problems
try:
  from ._pressiodemoappsimpl \
  import AdvectionDiffusionReaction2dProblem
except ImportError:
  raise ImportError("Unable to import AdvectionDiffusionReaction classes from _pressiodemoappsimpl")

# euler problems
try:
  from ._pressiodemoappsimpl \
    import Euler1dProblem #, Euler2dProblem, Euler3dProblem
except ImportError:
  raise ImportError("Unable to import Euler classes from _pressiodemoappsimpl")

# shallow-water problem
  from ._pressiodemoappsimpl \
    import Swe2dProblem
except ImportError:
  raise ImportError("Unable to import SWE classes from _pressiodemoappsimpl")

# functions to create problems
try:
  from ._pressiodemoappsimpl \
    import create_problem, \
    create_linear_advection_1d_problem, \
    create_diffusion_reaction_1d_problem_A,\
    create_burgers_2d_problem,\
    create_diffusion_reaction_2d_problem_A,\
    create_adv_diff_reac_2d_problem_A, \
    create_gray_scott_2d_problem,\
    create_slip_wall_swe_2d_problem,\
    create_cross_shock_problem
except ImportError:
  raise ImportError("Unable to import create_* from _pressiodemoappsimpl")


# Runge-Kutta2 integrator
# https://web.mit.edu/10.001/Web/Course_Notes/Differential_Equations_Notes/node5.html
def advanceRK2(appObj, state, dt, Nsteps, \
               startTime = 0.0, \
               observer = None, \
               showProgress=False):

  v = appObj.createRightHandSide()
  tmpState = state.copy()
  half = 0.5

  time = startTime
  for step in range(1, Nsteps+1):
    if showProgress:
      if step % 50 == 0: print("step = ", step, "/", Nsteps)

    appObj.rightHandSide(state, time, v)
    if observer!= None:
      observer(step-1, state, v)
    k1 = dt * v

    tmpState = state+k1
    appObj.rightHandSide(tmpState, time+dt, v)
    k2 = dt * v

    state[:] = state + k2*half + k1*half
    time += dt

# Runge-Kutta4 integrator
def advanceRK4(appObj, state, dt, Nsteps, \
               startTime = 0.0, \
               observer = None, \
               showProgress=False):

  v = appObj.createRightHandSide()
  tmpState = state.copy()
  half = 0.5
  two  = 2.
  oneOverSix = 1./6.

  time = startTime
  for step in range(1, Nsteps+1):
    if showProgress:
      if step % 50 == 0: print("step = ", step, "/", Nsteps)

    appObj.rightHandSide(state, time, v)
    if observer!= None:
      observer(step-1, state, v)
    k1 = dt * v

    tmpState = state+half*k1
    appObj.rightHandSide(tmpState, time+half*dt, v)
    k2 = dt * v

    tmpState = state+half*k2
    appObj.rightHandSide(tmpState, time+half*dt, v)
    k3 = dt * v

    tmpState = state+k3
    appObj.rightHandSide(tmpState, time+dt, v)
    k4 = dt * v

    state[:] = state + (k1+two*k2+two*k3+k4)*oneOverSix
    time += dt


# SSP3 integrator
def advanceSSP3(appObj, state, dt, Nsteps,\
                startTime = 0.0, \
                observer = None, \
                showProgress=False):

  v = appObj.createRightHandSide()
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

    appObj.rightHandSide(state, time, v)
    if observer!= None:
      observer(step-1, state, v)
    tmpState1[:] = state + dt * v

    appObj.rightHandSide(tmpState1, time+dt, v)
    tmpState2[:] = threeOverFour*state + oneOverFour*tmpState1 + oneOverFour*dt*v

    appObj.rightHandSide(tmpState2, time+dt*0.5, v)
    state[:] = oneOverThree*(state + two*tmpState2 + two*dt*v)

    time += dt
