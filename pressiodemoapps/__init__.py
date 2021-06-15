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
  from ._pressiodemoappsimpl import ReconstructionType, FluxType, Euler1d, Euler2d
except ImportError:
  raise ImportError("Unable to import enums from _pressiodemoappsimpl")

try:
  from ._pressiodemoappsimpl import LinearAdvection1dProblem, createPeriodicLinearAdvection1d
except ImportError:
  raise ImportError("Unable to import LinearAdvection1d from _pressiodemoappsimpl")

try:
  from ._pressiodemoappsimpl import Euler1dProblem, Euler2dProblem
except ImportError:
  raise ImportError("Unable to import Euler classes from _pressiodemoappsimpl")

try:
  from ._pressiodemoappsimpl import createProblem
except ImportError:
  raise ImportError("Unable to import createProblem from _pressiodemoappsimpl")
