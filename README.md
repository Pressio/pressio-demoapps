# pressio-demoapps

This repository contains a suite of 1D, 2D and 3D demo problems of varying complexity (from linear advection, to compressible Euler). <br>
The main feature of this work is the built-in support for hyper-reduction. This code, in fact, was originally started as part of the Pressio project to create a suite of benchmark problems to test ROMs and hyper-reduction techniques, but it is being developed to be self-contained, so it can be used for different purposes. For example, you can just use it for doing "standard" simulations, or you can just leverage the Python meshing scripts, or leverage the hyper-reduction capability to study function approximations. Its main scope is to provide a testbed of problems that are well-known to be hard for ROMs, explore the impact of the numerical scheme, and test news ideas for hyper-reduction.

Some features of this code are: uses cell-centered finite volumes, targets uniform structured meshes, the (simple) meshing tool is written in Python and it is separate from the main code, we rely on well-established shock-capturing schemes and flux forms (which can be easily extended), supports hyper-reduction/sample mesh in 1D, 2D and 3D for variny stencil size.

The development is grounded on:
- **simplicity**: we use high-level abstractions and a well-defined API to make this code as simple as possible to use, minimizing the number of steps needed to set up a problem: (1) create the mesh, (2) load the mesh and (3) create a problem instance that has a specific API.

- **generic programming**: the core component is a header-only C++ library based on generic programming, allowing us to write the code once and then use it with various data types;

- **prototyping**: fast prototyping is essential, so that is why we have developed Python bindings to make the code accessible to Python users, while maintaining the performance of the compiled C++ backend;

- **quality assurance**: we maintain a test suite for both the C++ and Python bindings to ensure the code remains stable and reliable;

- **extensibility**: one can easily add new problems, other numerical schemes (e.g. fluxes, reconstruction methods) and different data types.

Disclaimer/warning: this will **not** become a full solver. The scope is and will stay limited.


# High-level idea

To show our idea of usability more concretely, let's say that we are interested in exploring the Euler1d Sod problem.

**(Step 1)** Creating the mesh for this problem:
```py
python ./meshing_scripts/create_full_mesh_for.py --name sod1d_s7 --outdir ${HOME}/myTest -n 100
```
where the string `sod1d_s7` is composed of two parts: `sod1d` indicates the problem, and `s7` indicates the stencil size to use, in this case we want a 7-point stencil (the supported choices are discussed later); `${HOME}/myTest` is where all the mesh files are generated, and `-n=100` specifies how many cells we want. This mesh generation step is independent of the language you want to use. We now show what steps 2 and 3 look like in C++ and Python.

## C++ Synopsis

**(Step 2,3)** Load the mesh and create the problem (here we use Eigen data types which are currently suppoted):
```c++
#import <pressiodemoapps/euler1d.hpp>
// ...
namespace pda = pressiodemoapps;
const auto meshObj = pda::loadCellCenterUniformMeshEigen("/home/myTest");
const auto order   = pda::ReconstructionType::fifthOrderWeno;
auto problem       = pda::createProblemEigen(meshObj, pda::Euler1d::Sod, order);
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
More details on how to use this will be discussed later.


## Python Synopsis

The Python bindings of pressio-demoapps rely on numpy arrays, and the API is very similar:
```py
import pressiodemoapps as pda
# ...
meshObj  = pda.loadCellCenterUniformMesh("where-you-have-mesh-files")
order    = pda.ReconstructionType.fifthOrderWeno
appObj   = pda.createProblem(meshObj, pda.Euler1d.Sod, order)
```
All Python problem instances in *pressio-demoapps* meet the following API:
```py
class Problem:
  def initialCondition(self):
    #...
  def createVelocity():
    #...
  def velocity(state, time, velocity):
    # ...
```


# Installation

## Building Tests C++
Requires CMake > 3.18.0 and a C++ compiler with C++14 support:

```
git clone --recursive git@github.com:Pressio/pressio-demoapps.git

export CXX=<path-to-your-CXX-compiler> #must support C++14
cd pressio-demoapps && mkdir build && cd build
cmake -DPRESSIODEMOAPPS_ENABLE_TESTS=On -DPRESSIODEMOAPPS_ENABLE_TPL_EIGEN=ON ..
make -j4
ctest -j4
```
Note that to use the code as a library, you don't need to build anything.
Just point to the include directory.

# Python bindings
Requires CMake > 3.18.0 and a C++ compiler with C++14 support:

```
git clone --recursive git@github.com:Pressio/pressio-demoapps.git

export CXX=<path-to-your-CXX-compiler> #must support C++14
cd pressio-demoapps
python setup.py install

# run the tests to verify things work
pytest -s
```





<!-- ## 1d Problems: Sod, Lax -->

<!-- Mesh generation: -->
<!-- ```py -->
<!-- python ./meshing_scripts/create_full_mesh_for.py --name sod1d_s<3,7> -n <N> -outDir <somewhere> -->
<!-- python ./meshing_scripts/create_full_mesh_for.py --name lax1d_s<3,7> -n <N> -outDir <somewhere> -->
<!-- ``` -->

<!-- C++ problem object syntax: -->
<!-- ```c++ -->
<!-- namespace pda      = pressiodemoapps; -->
<!-- const auto meshObj = pda::loadCellCenterUniformMeshEigen(<mesh-path-string>); -->
<!-- const auto probId  = pda::euler1dproblemsEnum::{sod, lax} -->
<!-- const auto order   = pda::reconstructionEnum::{firstOrder, fifthOrderWeno}; -->
<!-- auto appObj        = pda::createEuler1dEigen(meshObj, order, probId); -->
<!-- ``` -->

<!-- Python object syntax: -->
<!-- ```py -->
<!-- meshO    = loadCellCenterUniformMesh(meshPath) -->
<!-- probId   = euler1d.{sod, lax} -->
<!-- appObj   = createEuler1dProblem(meshO, reconstructWith.fifthOrderWeno, probId) -->
<!-- ``` -->

<!-- ## 2d problems: Sedov, Riemann -->

<!-- Mesh generation: -->
<!-- ```py -->
<!-- python ./meshing_scripts/create_full_mesh_for.py --name sedov2d_s<3,7> -n <Nx> <Ny> -outDir <somewhere> -->
<!-- python ./meshing_scripts/create_full_mesh_for.py --name riemann2d_s<3,7> -n <Nx> <Ny> -outDir <somewhere> -->
<!-- ``` -->

<!-- C++ problem object syntax: -->
<!-- ```c++ -->
<!-- namespace pda      = pressiodemoapps; -->
<!-- const auto meshObj = pda::loadCellCenterUniformMeshEigen(<mesh-path-string>); -->
<!-- const auto probId  = pda::euler2dproblemsEnum::{sedov, riemann}; -->
<!-- const auto order   = pda::reconstructionEnum::{firstOrder, fifthOrderWeno}; -->
<!-- auto appObj        = pda::createEuler1dEigen(meshObj, order, probId); -->
<!-- ``` -->

<!-- Python object syntax: -->
<!-- ```py -->
<!-- meshO    = loadCellCenterUniformMesh(meshPath) -->
<!-- probId   = euler2d.{sedov, riemann} -->
<!-- appObj   = createEuler1dProblem(meshO, reconstructWith.fifthOrderWeno, probId) -->

<!-- ``` -->

# License and Citation

We are working on a publication of this work. <br>
**In the meantime, if you use this code, we kindly ask that you reference this github page!**

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The full license is available [here](https://pressio-demosapps.github.io/various/license/).


# Questions?
Find us on Slack: https://pressioteam.slack.com or open an issue on [github](https://github.com/Pressio/pressio-tutorials).
