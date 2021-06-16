# pressio-demoapps

This repository contains a suite of 1D, 2D and 3D demo problems of varying complexity (from linear advection, to compressible Euler). <br>
The main feature of this work is the built-in support for hyper-reduction. This code, in fact, was originally started as part of the Pressio project to create a suite of benchmark problems to test ROMs and hyper-reduction techniques, but it is being developed to be self-contained, so it can be used for different purposes. For example, you can just use it for doing "standard" simulations, or you can just leverage the Python meshing scripts, or leverage the hyper-reduction capability to study function approximations. Its main scope is to provide a testbed of problems that are well-known to be hard for ROMs, explore the impact of the numerical scheme, and test news ideas for hyper-reduction.

Some features of this code are: uses cell-centered finite volumes, targets uniform structured meshes, the (simple) meshing tool is written in Python and it is separate from the main code, supports well-established shock-capturing schemes and flux forms (which can be easily extended), supports hyper-reduction/sample mesh in 1D, 2D and 3D for varying stencil sizes.

The development is grounded on:
- **simplicity**: we use high-level abstractions and a well-defined API to make this code as simple as possible to use, minimizing the number of steps needed to set up a problem. We follow a *three-step* process: (1) generate the mesh, (2) load the mesh and (3) create a problem instance that has a specific API.

- **generic programming**: the core component is a header-only C++ library based on generic programming, thus allowing us to write the code once and then instantiate it with various data types;

- **prototyping**: fast prototyping is essential for research, so that motivated us to develop Python bindings. This makes the code accessible to Python users, while maintaining the performance of the compiled C++ backend;

- **quality assurance**: we maintain a test suite for both the C++ code and Python bindings to ensure stability and reliability;

- **extensibility**: one can easily add new problems, other numerical schemes (e.g. fluxes, reconstruction methods) and different data types.

Disclaimer/warning: this will **not** become a full solver. The scope is and will stay limited.

# A representative code snippet
```c++
#import <pressiodemoapps/euler1d.hpp>
int main(){
  namespace pda = pressiodemoapps;
  const auto meshObj   = pda::loadCellCenterUniformMeshEigen("/home/myTest");
  constexpr auto order = pda::ReconstructionType::fifthOrderWeno;
  auto problem         = pda::createProblemEigen(meshObj, pda::Euler1d::Sod, order);
}
```

# Usability in more detail

To show our idea of usability more concretely, let's say that we are interested in exploring the Euler1d Sod problem.

**(Step 1)** To create the mesh, we use the meshing scripts as follows:
```py
python ./meshing_scripts/create_full_mesh_for.py --problem sod1d_s7 --outdir ${HOME}/myTest -n 100
```
where the string `sod1d_s7` is composed of two parts: `sod1d` indicates the target problem, and `s7` indicates the stencil size to use, in this case we want a 7-point stencil (the supported choices are discussed later); `${HOME}/myTest` is where all the mesh files are generated, and `-n=100` specifies how many cells we want. This mesh generation step is independent of the language you want to use, and now we show what steps 2 and 3 look like in C++ and Python.

## C++ Synopsis

**(Step 2,3)** Load the mesh and create the problem (here we use Eigen data types which are currently suppoted):
```c++
#import <pressiodemoapps/euler1d.hpp>
// ...
namespace pda = pressiodemoapps;
const auto meshObj   = pda::loadCellCenterUniformMeshEigen("/home/myTest");
constexpr auto order = pda::ReconstructionType::fifthOrderWeno;
auto problem         = pda::createProblemEigen(meshObj, pda::Euler1d::Sod, order);
// ...
```
This will create an instance of the Sod1d problem and selects the 5-th order WENO for doing the edge reconstruction (more details on the schemes and stencils are given below). All C++ problem instances in *pressio-demoapps* meet the following C++ API:
```cpp
class Problem
{
public:
  using scalar_type   = /* type alias of the scalar (e.g. double or float) */;
  using state_type    = /* type alias of the state vector */;
  using velocity_type = /* type alias of the RHS vector */;

  state_type    initialCondition() const;
  velocity_type createVelocity() const;
  void velocity(const state_type & state, scalar_type time, velocity_type & velocity);
}
```
The `initialCondition` returns a state vector initialized for the target problem.
This API allows you to have full control on how to solve the problem. For example, you can use custom time integrators. More details on how to use this will be discussed later.

## Python Synopsis

The Python bindings of pressio-demoapps rely on numpy arrays, and the API is very similar:
```py
import pressiodemoapps as pda
# ...
meshObj = pda.loadCellCenterUniformMesh("where-you-have-mesh-files")
order   = pda.ReconstructionType.fifthOrderWeno
appObj  = pda.createProblem(meshObj, pda.Euler1d.Sod, order)
```
All Python problem instances in *pressio-demoapps* meet the following API:
```py
class Problem:
  def initialCondition(self):          # returns initial conditon
  def createVelocity():                # returns an empty velocity vector
  def velocity(state, time, velocity): # computes velocity for given state and time
```

# List of Supported Problems

|                     | C++ Enum                | Python Enum            | Supported <br> reconstruction scheme |     Flux scheme  |   |
|---------------------|-------------------------|------------------------|:----------------------:|:----------------:|---|
| 1D Linear Advection | Advection1d::PeriodicLinear  | Advection1d.PeriodicLinear                      |        all             |   Rusanov        | [Descrioption](https://github.com/Pressio/pressio-demoapps/wiki/1D-Problems-descriptions#linear-advection)  |
| 1D Euler Smooth     | Euler1d::PeriodicSmooth | Euler1d.PeriodicSmooth |        all             |   Rusanov        | [Descrioption](https://github.com/Pressio/pressio-demoapps/wiki/1D-Problems-descriptions#euler-smooth)  |
| 1D Sod              | Euler1d::Sod            | Euler1d.Sod            |        all             |   Rusanov        | [Descrioption](https://github.com/Pressio/pressio-demoapps/wiki/1D-Problems-descriptions#sod)  |
| 1D Lax              | Euler1d::Lax            | Euler1d.Lax            |        all             |   Rusanov        | [Descrioption](https://github.com/Pressio/pressio-demoapps/wiki/1D-Problems-descriptions#lax)  |
| 2d Euler Smooth     | Euler2d::PeriodicSmooth | Euler2d.PeriodicSmooth |        all             |   Rusanov        | [Descrioption](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#euler-smooth)  |
| 2D Sedov (full)     | Euler2d::SedovFull      | Euler2d.SedovFull      |        all             |   Rusanov        | [Descrioption](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#sedov-version-1-full-version-no-symmetry-use)  |
| 2D Sedov (symmetry) | Euler2d::SedovSymmetry  | Euler2d.SedovSymmetry  |        all             |   Rusanov        | [Descrioption](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#sedov-version-2-exploit-symmetry-to-simulate-only-one-quadrant)  |
| 2D Riemann          | Euler2d::Riemann        | Euler2d.Riemann        |        all             |   Rusanov        | [Descrioption](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#riemann)  |
| 3d Euler Smooth     | Euler3d::PeriodicSmooth | Euler3d.PeriodicSmooth |     1st-order          |   Rusanov        | [Descrioption](https://github.com/Pressio/pressio-demoapps/wiki/3D-Problems-Description#euler-smooth) |
| 3d Sedov (symmetry) | Euler3d::SedovSymmetry  | Euler3d.SedovSymmetry  |     1st-order          |   Rusanov        | [Descrioption](https://github.com/Pressio/pressio-demoapps/wiki/3D-Problems-Description#sedov-exploit-symmetry-to-simulate-only-18)  |

Currently supported reconstruction schemes: `ReconstructionType::firstOrder`, and `ReconstructionType::fifthOrderWeno`.

# Installation

## C++ library
The C++ library is header-only so it does not need to be compiled and installed.
To use the headers, all you have to do is to include in your compilation line the `pressiodemoapps/include` subdirectory, specify this `-DPRESSIODEMOAPPS_ENABLE_TPL_EIGEN=ON` to CMake while building your code, and also include the Eigen library (whose headers you can find inside `pressiodemoapps/tpls`).
Note that you need a C++ compiler with support for C++14.

If you want to build the C++ tests, you need CMake > 3.18.0 and a C++ compiler with C++14 support:

```
git clone --recursive git@github.com:Pressio/pressio-demoapps.git

export CXX=<path-to-your-CXX-compiler> #must support C++14
cd pressio-demoapps && mkdir build && cd build
cmake -DPRESSIODEMOAPPS_ENABLE_TESTS=On -DPRESSIODEMOAPPS_ENABLE_TPL_EIGEN=ON ..
make -j4
ctest -j4
```

## Python bindings
Requires CMake > 3.18.0 and a C++ compiler with C++14 support, but should be as easy as:

```
git clone --recursive git@github.com:Pressio/pressio-demoapps.git
export CXX=<path-to-your-CXX-compiler> #must support C++14
cd pressio-demoapps
python setup.py install

# run the tests to verify things work
pytest -s
```


# Sample Mesh

Todo.

## 1D Example
![Sample mesh example for 1D](https://github.com/Pressio/pressio-demoapps/blob/develop/figures/full.png)




# License and Citation

We are working on a publication of this work. <br>
**In the meantime, if you use this code, we kindly ask that you reference this github page!**

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The full license is available [here](https://pressio-demosapps.github.io/various/license/).


# Questions?
Find us on Slack: https://pressioteam.slack.com or open an issue on [github](https://github.com/Pressio/pressio-tutorials).
