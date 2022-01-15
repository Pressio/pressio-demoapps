.. highlight:: sh

Installation
============

This page walks you through the process of installing pressio-demoapps for both C++ and Python.


C++ library
-----------

The C++ library is header-only so it does not need to be compiled and installed.

To use it, you need a C++14 compiler and you have to:

1. include the ``pressiodemoapps/include`` subdirectory in your compilation line

2. include the Eigen library (whose headers you can find inside ``pressiodemoapps/tpls``).

3. specify the CMake option ``-DPRESSIODEMOAPPS_ENABLE_TPL_EIGEN=ON`` while building your code.

Building the test suite
~~~~~~~~~~~~~~~~~~~~~~~

If you want to build the C++ tests, you need CMake > 3.18.0 and then do::

  git clone --recursive git@github.com:Pressio/pressio-demoapps.git
  export CXX=<path-to-your-CXX-compiler> #must support C++14
  cd pressio-demoapps && mkdir build && cd build
  cmake -DPRESSIODEMOAPPS_ENABLE_TESTS=On -DPRESSIODEMOAPPS_ENABLE_TPL_EIGEN=ON ..
  make -j4
  ctest -j4


Python library
--------------

To build/install the bindings, you need:

- CMake > 3.18.0 and a C++ compiler with C++14 support

- You also need these packages::

    pip3 install build pytest numpy scipy matplotlib

Then, you can do:

.. code-block:: shell

   git clone --recursive git@github.com:Pressio/pressio-demoapps.git
   export CXX=<path-to-your-CXX-compiler> #must support C++14
   cd pressio-demoapps
   python3 build.py
   pip3 install .


This builds/installs pressiodemoapps with default options (build_mode=Release).



..
   To build/install pressiodemoapps with OpenMP and Release mode:
   git clone --recursive git@github.com:Pressio/pressio-demoapps.git
   export CXX=<path-to-your-CXX-compiler> #must support C++14
   python3 build.py --openmp
   pip3 install .
   # to just build do
   python -m build


You can run the tests to verify things::

  cd pressio-demoapps
  pytest -s
