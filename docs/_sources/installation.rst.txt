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

If you want to build the C++ tests, you need CMake > 3.18.0::

  git clone --recursive git@github.com:Pressio/pressio-demoapps.git
  export CXX=<path-to-your-CXX-compiler> #must support C++14
  cd pressio-demoapps && mkdir build && cd build
  cmake -DPRESSIODEMOAPPS_ENABLE_TESTS=On -DPRESSIODEMOAPPS_ENABLE_TPL_EIGEN=ON ..
  make -j4
  ctest -j4


Python library
--------------

Requires:

* CMake > 3.18.0 and C++ compiler with C++14 support

Should be as easy as::

  git clone --recursive git@github.com:Pressio/pressio-demoapps.git
  export CXX=<path-to-your-CXX-compiler> #must support C++14
  cd pressio-demoapps
  python setup.py install

You can run the tests to verify things::

  cd pressio-demoapps
  pytest -s
