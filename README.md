# pressio-demoapps

This repository contains a suite of 1D, 2D and 3D demo problems of varying complexity (from linear advection, to compressible Euler). <br>
The main feature of this work is the built-in support for hyper-reduction. This code, in fact, was originally started as part of the Pressio project to create a suite of benchmark problems to test ROMs and hyper-reduction techniques, but it is being developed to be self-contained, so it can be used for different purposes. For example, you can just use it for doing "standard" simulations, or you can just leverage the Python meshing scripts, or leverage the hyper-reduction capability to study function approximations. Its main scope is to provide a testbed of problems that are well-known to be hard for ROMs, explore the impact of the numerical scheme, and test news ideas for hyper-reduction.

Some features of this code are:
- a cell-centered finite volume discretization on *uniform structured meshes*
- a (simple) meshing tool written in Python that is separate from the main code
- well-established shock-capturing schemes and flux forms (which can be easily extended)
- hyper-reduction/sample mesh in 1D, 2D and 3D for varying stencil sizes.

The development is grounded on:
- **simplicity**: we use high-level abstractions and a well-defined API to make this code as simple as possible to use, minimizing the number of steps needed to set up a problem. We follow a *three-step* process: (1) generate the mesh, (2) load the mesh and (3) create a problem instance that has a specific API.

- **generic programming**: the core component is a header-only C++ library based on generic programming, thus allowing us to write the code once and then instantiate it with various data types;

- **prototyping**: fast prototyping is essential for research, so that motivated us to develop Python bindings. This makes the code accessible to Python users, while maintaining the performance of the compiled C++ backend;

- **quality assurance**: we maintain a test suite for both the C++ code and Python bindings to ensure stability and reliability;

- **extensibility**: one can easily add new problems, other numerical schemes (e.g. fluxes, reconstruction methods) and different data types.

Disclaimer/warning: this will **not** become a full solver. The scope is and will stay limited.

### The main idea in 4 lines
```c++
#import <pressiodemoapps/euler1d.hpp>
int main(){
  namespace pda = pressiodemoapps;
  const auto meshObj   = pda::loadCellCenterUniformMeshEigen("<path-to-mesh-files>");
  constexpr auto order = pda::InviscidFluxReconstruction::Weno5;
  auto problem         = pda::createProblemEigen(meshObj, pda::Euler1d::Sod, order);
}
```

___


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
constexpr auto order = pda::InviscidFluxReconstruction::Weno5;
auto problem         = pda::createProblemEigen(meshObj, pda::Euler1d::Sod, order);
// ...
```
This will create an instance of the Sod1d problem and select the 5-th order WENO for doing the edge reconstruction (more details on the schemes and stencils are given below). All C++ problem instances in *pressio-demoapps* meet the following C++ API:
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
order   = pda.InviscidFluxReconstruction.Weno5;
appObj  = pda.createProblem(meshObj, pda.Euler1d.Sod, order)
```
All Python problem instances in *pressio-demoapps* meet the following API:
```py
class Problem:
  def initialCondition(self):          # returns initial conditon
  def createVelocity():                # returns an empty velocity vector
  def velocity(state, time, velocity): # computes velocity for given state and time
```

___


# Installation

## C++ library
The C++ library is header-only so it does not need to be compiled and installed.
To use it, you need a C++14 compiler and you have to:
1. include the `pressiodemoapps/include` subdirectory in your compilation line
2. include the Eigen library (whose headers you can find inside `pressiodemoapps/tpls`).
3. specify the CMake option `-DPRESSIODEMOAPPS_ENABLE_TPL_EIGEN=ON` while building your code.

If you want to build the C++ tests, you need CMake > 3.18.0:

```
git clone --recursive git@github.com:Pressio/pressio-demoapps.git

export CXX=<path-to-your-CXX-compiler> #must support C++14
cd pressio-demoapps && mkdir build && cd build
cmake -DPRESSIODEMOAPPS_ENABLE_TESTS=On -DPRESSIODEMOAPPS_ENABLE_TPL_EIGEN=ON ..
make -j4
ctest -j4
```

## Python bindings
Requires CMake > 3.18.0, a C++ compiler with C++14 support, and should be as easy as:

```
git clone --recursive git@github.com:Pressio/pressio-demoapps.git
export CXX=<path-to-your-CXX-compiler> #must support C++14
cd pressio-demoapps
python setup.py install

# run the tests to verify things work
pytest -s
```



# List of Supported Problems

|                     | Enum Identifier         | Reconstruct scheme |     Flux scheme  |   |
|---------------------|-------------------------|:----------------------:|:----------------:|---|
| 1D Linear Advection | C++: Advection1d::PeriodicLinear <br> Py &nbsp; : Advection1d.PeriodicLinear  |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/1D-Problems-descriptions#linear-advection)  |
| 1D Euler Smooth     | C++: Euler1d::PeriodicSmooth <br> Py &nbsp; : Euler1d.PeriodicSmooth|        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/1D-Problems-descriptions#euler-smooth)  |
| 1D Sod              | C++: Euler1d::Sod <br> Py &nbsp; : Euler1d.Sod           |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/1D-Problems-descriptions#sod)  |
| 1D Lax              | C++: Euler1d::Lax <br> Py &nbsp; : Euler1d.Lax           |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/1D-Problems-descriptions#lax)  |
| 2d Euler Smooth     | C++: Euler2d::PeriodicSmooth <br> Py &nbsp; : Euler2d.PeriodicSmooth |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#euler-smooth)  |
| 2D Sedov (full)     | C++: Euler2d::SedovFull <br> Py &nbsp; : Euler2d.SedovFull      |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#sedov-version-1-full-version-no-symmetry-use)  |
| 2D Sedov (symmetry) | C++: Euler2d::SedovSymmetry <br> Py &nbsp; : Euler2d.SedovSymmetry  |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#sedov-version-2-exploit-symmetry-to-simulate-only-one-quadrant)  |
| 2D Riemann          | C++: Euler2d::Riemann <br> Py &nbsp; : Euler2d.Riemann        |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#riemann)  |
| 3d Euler Smooth     | C++: Euler3d::PeriodicSmooth <br> Py &nbsp; : Euler3d.PeriodicSmooth |     1st-order          |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/3D-Problems-Description#euler-smooth) |
| 3d Sedov (symmetry) | C++: Euler3d::SedovSymmetry <br> Py &nbsp; :  Euler3d.SedovSymmetry |     1st-order          |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/3D-Problems-Description#sedov-exploit-symmetry-to-simulate-only-18)  |

Reconstruction schemes currently available: `InviscidFluxReconstruction::FirstOrder`, and `InviscidFluxReconstruction::Weno5`.

<!--
# The meshing scripts
The script to generate the mesh are located in `pressiodemoapps/meshing_scripts`.
There are three main scripts:
* `create_full_mesh.py`: script to create the full mesh for an arbitrary domain.
```py
# 1D mesh
python create_full_mesh.py \
  --outDir <where-to-generate-files> \
  -n <numCells> \
  -stencilSize
``` -->

<!-- * `create_full_mesh_for.py`: script to create the full mesh for a specific problem (see above)
* `create_sample_mesh.py`: script to create the sample mesh, more details on this below -->



# Sample Mesh

The sample mesh is a disjoint collection of cells at which the velocity or residual vector are computed. This collection of cells is determined by a hyper-reduction approaches such as the discrete empirical interpolation method ([DEIM](https://doi.org/10.1137/090766498)), or Gauss-Newton with approximate tensors ([GNAT](https://doi.org/10.1016/j.jcp.2013.02.028)). The sample mesh is a crucial part of hyper-reduction implementations, allowing one to create ROMs with a computational cost which *does not* scale with the size of the full model's state vector. 

The sample mesh is used in conjunction with what we refer to as the **stencil mesh**, which contains all cells needed to compute the velocity or residual vector on the sample mesh. Several examples of sample and stencil meshes are shown below, along with the indexing scheme used in this repository. 

### 1D Example
<img src="https://github.com/Pressio/pressio-demoapps/blob/develop/figures/1dfull.png" width="75%">
<!-- ![Sample mesh example for 1D](https://github.com/Pressio/pressio-demoapps/blob/develop/figures/full.png) -->

The first example shows a one-dimensional sample mesh and the corresponding stencil mesh for a first order cell-centered finite volume scheme. The sample mesh cells are shown in green, while the stencil mesh cells are not colored. All cells are labeled with their ID number; note that the sample/stencil mesh cell IDs **do not** match the full mesh IDs, so one must keep track of the mapping between the IDs.  

### 2D Example
<img src="https://github.com/Pressio/pressio-demoapps/blob/develop/figures/2dsm.png" width="75%">

The second example shows a two dimensional sample mesh and a stencil mesh for a first order cell-centered finite volume scheme. The coloring scheme is the same as in the first example. 

### 3D Example
<img src="https://github.com/Pressio/pressio-demoapps/blob/develop/figures/3dsm.png" width="75%">

The third example shows a three dimensional sample mesh and a stencil mesh for a first order cell-centered finite volume scheme. The coloring scheme is the same as in the previous examples. Cell IDs are omitted for clarity. 

# License and Citation

We are working on a publication of this work. <br>
**In the meantime, if you use this code, we kindly ask that you reference this github page!**

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The full license is available [here](https://pressio-demosapps.github.io/various/license/).


# Questions?
Find us on Slack: https://pressioteam.slack.com or open an issue on [github](https://github.com/Pressio/pressio-tutorials).
