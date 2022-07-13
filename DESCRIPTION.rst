pressio-demoapps: Python bindings
=================================

This package provides Python bindings for the C++ library pressio-demoapps (website_).

.. _website: https://pressio.github.io/pressio-demoapps/index.html


Install
-------

You need a C++14-compliant compiler and CMake > 3.18.0:

.. code-block:: bash

  export CXX=<path-to-your-C++-compiler>
  cd pressio-demoapps
  python3 cmake_build.py
  pip install .


You can double check that everything worked fine by doing:

.. code-block:: python

  import pressiodemoapps
  print(pressiodemoapps.__version__)

After installing the library, you can run the regression tests:

.. code-block:: python

  cd pressio-demoapps
  pytest -s


Citations
---------

If you use this package, please acknowledge our work-in-progress:

TBD
