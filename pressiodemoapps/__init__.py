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
  from ._pressiodemoappsimpl import InviscidFluxReconstruction, InviscidFluxScheme, \
    Advection1d, Euler1d, Euler2d, Euler3d
except ImportError:
  raise ImportError("Unable to import enums from _pressiodemoappsimpl")

try:
  from ._pressiodemoappsimpl import Advection1dProblem
except ImportError:
  raise ImportError("Unable to import Advection classes from _pressiodemoappsimpl")

try:
  from ._pressiodemoappsimpl import Euler1dProblem, Euler2dProblem, Euler3dProblem
except ImportError:
  raise ImportError("Unable to import Euler classes from _pressiodemoappsimpl")

try:
  from ._pressiodemoappsimpl import createProblem
except ImportError:
  raise ImportError("Unable to import createProblem from _pressiodemoappsimpl")
