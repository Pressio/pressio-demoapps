
'''
see this for why this file exists and is done this way
https://stackoverflow.com/questions/47599162/pybind11-how-to-package-c-and-python-code-into-a-single-package?rq=1
'''
try:
  from ._pressiodemoappsimpl \
  import Euler1dProblem, createEuler1dProblem, \
  		 Euler2dProblem, createEuler2dProblem
except ImportError:
  raise ImportError("Unable to import from _pressiodemoappsimpl")
