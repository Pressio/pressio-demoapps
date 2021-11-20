
# Overview

**pressio-demoapps** is a library of 1D, 2D and 3D demo problems
of varying complexity, ranging from a simple 1D linear advection,
to 2D reaction-diffusion, and 3D compressible Euler, and more.
Key features include:

- support for both C++ and Python
- cell-centered finite volume discretization with various numerical schemes and *exact Jacobians*
- focus on providing self-contained and well-defined problems
- built-in support for a sample mesh: this mean that one can evaluate the residual and Jacobian
at a disjoint subset of the mesh cells (this is useful for intrusive ROMs)

Check the documentation website for more details:

<a href="https://pressio.github.io/pressio-demoapps/index.html" target="_blank">
    <img src='figures/logo-display.svg' width='90%'>
</a>

## Development status

**pressio-demoapps** is planned to be a long-term maintained project.
Therefore, more problems and features will be implemented.
If you are interested in collaborating or would like to see a specific
problem added, please reach out.

## Questions?
Find us on Slack: https://pressioteam.slack.com or open an issue on [github](https://github.com/Pressio/pressio-demoapps).

## License and Citation

While we work on publishing this, **if you use this code, please reference this Github page!**

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The full license is available [here](https://github.com/Pressio/pressio-demoapps/blob/main/LICENSE).



<!-- which is necessary for hyper-reduction, more details below). This code, in fact, was originally started as part of the Pressio project to create a suite of benchmark problems to test ROMs and hyper-reduction techniques, but it is being developed to be self-contained, so it can be used for different purposes. For example, you can just use it for doing "standard" simulations, or you can just use the Python meshing scripts, or leverage the sample mesh capability to study function approximations. The main objective is to provide a testbed of problems that are known to be hard for ROMs, explore the impact of the numerical scheme, and test new sample mesh and hyper-reduction techniques. -->

<!-- Some features of this code are: -->
<!-- - a cell-centered finite volume discretization on *uniform structured meshes* -->
<!-- - a (simple) meshing tool written in Python that also handles sample meshes -->
<!-- - well-established shock-capturing schemes and flux forms (which can be easily extended) -->
<!-- - sample mesh in 1D, 2D and 3D for varying stencil sizes. -->

<!-- The development is grounded on: -->
<!-- - **simplicity**: we use high-level abstractions and a well-defined API to make this code as simple as possible to use, minimizing the number of steps needed to set up a problem. We follow a *three-step* process: (1) generate the mesh, (2) load the mesh and (3) create a problem instance that has a specific API. -->

<!-- - **generic programming**: the core component is a header-only C++ library based on generic programming, thus allowing us to write the code once and then instantiate it with various data types; -->

<!-- - **prototyping**: fast prototyping is essential for research, so that motivated us to develop Python bindings. This makes the code accessible to Python users, while maintaining the performance of the compiled C++ backend; -->

<!-- - **quality assurance**: we maintain a test suite for both the C++ code and Python bindings to ensure stability and reliability; -->

<!-- - **extensibility**: easily add new problems, other numerical schemes (e.g. fluxes, reconstruction methods) and different data types. -->

<!-- Disclaimer/warning: this will **not** become a full solver. The scope is and will stay limited. -->

<!-- ### The main idea in 4 lines -->
<!-- ```c++ -->
<!-- #import <pressiodemoapps/euler1d.hpp> -->
<!-- int main(){ -->
<!--   namespace pda = pressiodemoapps; -->
<!--   const auto meshObj   = pda::loadCellCenterUniformMeshEigen("<path-to-mesh-files>"); -->
<!--   constexpr auto order = pda::InviscidFluxReconstruction::Weno5; -->
<!--   auto problem         = pda::createProblemEigen(meshObj, pda::Euler1d::Sod, order); -->
<!-- } -->
<!-- ``` -->

<!-- <br> -->

<!-- # Usability in more detail -->

<!-- Let's say that we are interested in exploring the Euler1d Sod problem. <br> -->
<!-- As mentioned above, we envision a three-step process: generate the mesh, -->
<!-- instantiate the problem object, and solve. -->

<!-- ## Step 1: Generating the Mesh -->

<!-- We provide a meshing script that works as follows: -->
<!-- ```py -->
<!-- # assuming we are inside the pressio-demoapps dir -->
<!-- python ./meshing_scripts/create_full_mesh_for.py --problem sod1d_s7 --outdir ${HOME}/myTest -n 100 -->
<!-- ``` -->
<!-- where the string `sod1d_s7` is composed of two parts: `sod1d` indicates the target problem, and `s7` indicates the stencil size to use, in this case we want a 7-point stencil (the supported choices are discussed later); `${HOME}/myTest` is where all the mesh files are generated, and `-n=100` specifies how many cells we want. -->
<!-- The advantage of this script is that any information about the problem domain -->
<!-- and other details are encoded in the script, so it only exposes the minimal set of parameters -->
<!-- (e.g. number of cells) to set. -->
<!-- <\!-- This mesh generation step is independent of the language you want to use, and now we show what steps 2 and 3 look like in C++ and Python. -\-> -->


<!-- ## Step 2: Creating a problem instance -->

<!-- This second step uses the mesh generated above and creates an instance of the Sod1d problem.<br> -->
<!-- We now show what this step 2 looks like in C++ and Python. -->

<!-- **C++ Synopsis** -->
<!-- ```c++ -->
<!-- #import <pressiodemoapps/euler1d.hpp> -->
<!-- // ... -->
<!-- namespace pda = pressiodemoapps; -->
<!-- const auto meshObj   = pda::loadCellCenterUniformMeshEigen("/home/myTest"); -->
<!-- constexpr auto order = pda::InviscidFluxReconstruction::Weno5; -->
<!-- auto problem         = pda::createProblemEigen(meshObj, pda::Euler1d::Sod, order); -->
<!-- // ... -->
<!-- ``` -->
<!-- This creates an instance of the Sod1d problem using Eigen data types, and selects the 5-th order WENO for the inviscid flux edge reconstruction (more details on the schemes and stencils are given below). Note that here we explicitly ask for an instance of the problem that is based -->
<!-- on Eigen data types. Eigen is currently the main backend supported, but other ones, -->
<!-- e.g., Kokkos, will be added later. -->

<!-- **Python Synopsis** -->
<!-- ```py -->
<!-- import pressiodemoapps as pda -->
<!-- # ... -->
<!-- meshObj = pda.loadCellCenterUniformMesh("/home/myTest") -->
<!-- order   = pda.InviscidFluxReconstruction.Weno5; -->
<!-- problem = pda.createProblem(meshObj, pda.Euler1d.Sod, order) -->
<!-- ``` -->
<!-- In Python, this creates a problem based on `numpy`. -->
<!-- We remark that we aimed for the Python code to be readable as similarly as possible to the C++. -->


<!-- ## Step 3: Solving the problem -->
<!-- To use a problem instance, you need to know that all *pressio-demoapps* problem instances meet -->
<!-- a specific API.<br> -->
<!-- If you are using C++, this API is: -->
<!-- ```cpp -->
<!-- class Problem -->
<!-- { -->
<!-- public: -->
<!--   using scalar_type   = /* type alias of the scalar (e.g. double or float) */; -->
<!--   using state_type    = /* type alias of the state vector */; -->
<!--   using velocity_type = /* type alias of the RHS vector */; -->

<!--   state_type    initialCondition() const; -->
<!--   velocity_type createVelocity() const; -->
<!--   void velocity(const state_type & state, scalar_type time, velocity_type & velocity); -->
<!-- } -->
<!-- ``` -->
<!-- In Python, the API is: -->
<!-- ```py -->
<!-- class Problem: -->
<!--   def initialCondition(self):          # returns initial conditon -->
<!--   def createVelocity():                # returns an empty velocity vector -->
<!--   def velocity(state, time, velocity): # computes velocity for given state and time -->
<!-- ``` -->

<!-- The API is simple, and exposes three main methods (the computation of the Jacobian is in progress): -->
<!-- - `initialCondition`: returns a state object that is *correctly* initialized for the target problem you want; -->
<!-- - `createVelocity`: returns an empty velocity object; -->
<!-- - `velocity`: computes that system's velocity for a given state and time. -->

<!-- This API allows you to have full control on how to solve the problem. For example, you can use custom time integrators. More details on how to use this will be discussed later. -->

<!-- ### Do you have ready-to-use functions to solve the problem? We do. -->
<!-- We provide ready-to-use functions that you can use to solve the problem. -->

<!-- **C++ Synopsis** -->

<!-- TODO. -->

<!-- **Python Synopsis** -->
<!-- ```py -->
<!-- import pressiodemoapps as pda -->
<!-- ... -->
<!-- problem = # instance of a problem created as in step 2 -->
<!-- state   = problem.initialCondition() -->
<!-- pda.advanceSSP3(problem, state, timeStepSize, numSteps) -->
<!-- ``` -->
<!-- This will use the SSP3 time integrator to solve the problem, -->
<!-- and you need to pass the time step size and number of steps. -->
<!-- At the end, the `state` object will contain the solution at the final time. -->


<!-- <br> -->

<!-- # Installation -->

<!-- ## C++ library -->
<!-- The C++ library is header-only so it does not need to be compiled and installed. -->
<!-- To use it, you need a C++14 compiler and you have to: -->
<!-- 1. include the `pressiodemoapps/include` subdirectory in your compilation line -->
<!-- 2. include the Eigen library (whose headers you can find inside `pressiodemoapps/tpls`). -->
<!-- 3. specify the CMake option `-DPRESSIODEMOAPPS_ENABLE_TPL_EIGEN=ON` while building your code. -->

<!-- If you want to build the C++ tests, you need CMake > 3.18.0: -->

<!-- ``` -->
<!-- git clone --recursive git@github.com:Pressio/pressio-demoapps.git -->
<!-- export CXX=<path-to-your-CXX-compiler> #must support C++14 -->
<!-- cd pressio-demoapps && mkdir build && cd build -->
<!-- cmake -DPRESSIODEMOAPPS_ENABLE_TESTS=On -DPRESSIODEMOAPPS_ENABLE_TPL_EIGEN=ON .. -->
<!-- make -j4 -->
<!-- ctest -j4 -->
<!-- ``` -->

<!-- ## Python bindings -->
<!-- Requires CMake > 3.18.0, a C++ compiler with C++14 support, and should be as easy as: -->

<!-- ``` -->
<!-- git clone --recursive git@github.com:Pressio/pressio-demoapps.git -->
<!-- export CXX=<path-to-your-CXX-compiler> #must support C++14 -->
<!-- cd pressio-demoapps -->
<!-- python setup.py install -->

<!-- # run the tests to verify things work -->
<!-- pytest -s -->
<!-- ``` -->

<!-- <br> -->


<!-- # List of Supported Problems -->

<!-- |                           | Enum Identifier         | Inviscid <br> Reconstruct scheme |  Inviscid <br> Flux scheme  |   | -->
<!-- |---------------------------|-------------------------|:----------------------:|:----------------:|---| -->
<!-- | 1D Linear Advection       | C++: Advection1d::PeriodicLinear <br> Py &nbsp; : Advection1d.PeriodicLinear  |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/1D-Problems-description#linear-advection)  | -->
<!-- | 1D Euler Smooth           | C++: Euler1d::PeriodicSmooth <br> Py &nbsp; : Euler1d.PeriodicSmooth|        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/1D-Problems-description#euler-smooth)  | -->
<!-- | 1D Sod                    | C++: Euler1d::Sod <br> Py &nbsp; : Euler1d.Sod           |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/1D-Problems-description#sod)  | -->
<!-- | 1D Lax                    | C++: Euler1d::Lax <br> Py &nbsp; : Euler1d.Lax           |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/1D-Problems-description#lax)  | -->
<!-- | 2D Euler Smooth           | C++: Euler2d::PeriodicSmooth <br> Py &nbsp; : Euler2d.PeriodicSmooth |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#euler-smooth)  | -->
<!-- | 2D Euler Kelvin Helmholtz | C++: Euler2d::KelvinHelmholtz <br> Py &nbsp; : Euler2d.KelvinHelmholtz |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#kelvin-helmholtz)  | -->
<!-- | 2D Sedov (full)           | C++: Euler2d::SedovFull <br> Py &nbsp; : Euler2d.SedovFull      |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#sedov-version-1-full-version-no-symmetry-use)  | -->
<!-- | 2D Sedov (symmetry)       | C++: Euler2d::SedovSymmetry <br> Py &nbsp; : Euler2d.SedovSymmetry  |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#sedov-version-2-exploit-symmetry-to-simulate-only-one-quadrant)  | -->
<!-- | 2D Riemann                | C++: Euler2d::Riemann <br> Py &nbsp; : Euler2d.Riemann        |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#riemann)  | -->
<!-- | 2D Double Mach Reflection | C++: Euler2d::DoubleMachReflection <br> Py &nbsp; : Euler2d.DoubleMachReflection   |   firstOrder, Weno3  |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#double-mach-reflection)  | -->
<!-- | 2D Shallow Water          | C++: Swe2d::SlipWall <br> Py &nbsp; : Swe2d.SlipWall        |        all             |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/2D-Problems-Description#shallow-water-equations)  | -->
<!-- | 3d Euler Smooth           | C++: Euler3d::PeriodicSmooth <br> Py &nbsp; : Euler3d.PeriodicSmooth |     FirstOrder, Weno3          |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/3D-Problems-Description#euler-smooth) | -->
<!-- | 3d Sedov (symmetry)       | C++: Euler3d::SedovSymmetry <br> Py &nbsp; :  Euler3d.SedovSymmetry |     firstOrder, Weno3          |   Rusanov        | [Details](https://github.com/Pressio/pressio-demoapps/wiki/3D-Problems-Description#sedov-exploit-symmetry-to-simulate-only-18)  | -->

<!-- Reconstruction schemes currently available: -->
<!-- - `InviscidFluxReconstruction::FirstOrder`: needs at least a 3pt stencil -->
<!-- - `InviscidFluxReconstruction::Weno3` [from Jiang and Shu](https://www.sciencedirect.com/science/article/pii/S0021999196901308): needs at least a 5pt stencil -->
<!-- - `InviscidFluxReconstruction::Weno5` [from Jiang and Shu](https://www.sciencedirect.com/science/article/pii/S0021999196901308): needs at least a 7pt stencil -->

<!-- <\!-- -->
<!-- # The meshing scripts -->
<!-- The script to generate the mesh are located in `pressiodemoapps/meshing_scripts`. -->
<!-- There are three main scripts: -->
<!-- * `create_full_mesh.py`: script to create the full mesh for an arbitrary domain. -->
<!-- ```py -->
<!-- # 1D mesh -->
<!-- python create_full_mesh.py \ -->
<!--   --outDir <where-to-generate-files> \ -->
<!--   -n <numCells> \ -->
<!--   -stencilSize -->
<!-- ``` -\-> -->

<!-- <\!-- * `create_full_mesh_for.py`: script to create the full mesh for a specific problem (see above) -->
<!-- * `create_sample_mesh.py`: script to create the sample mesh, more details on this below -\-> -->


<!-- <br> -->

<!-- # Sample Mesh -->

<!-- In practice, the *sample mesh is a disjoint collection of cells at which the velocity or residual vector are computed*. Several methods exist to determine which cells to include, e.g., random sampling, the discrete empirical interpolation method ([DEIM](https://doi.org/10.1137/090766498)), and Gauss-Newton with approximate tensors ([GNAT](https://doi.org/10.1016/j.jcp.2013.02.028)). -->
<!-- The sample mesh is critical to make projection-based ROMs of nonlinear systems practical. -->
<!-- Its at the core of hyper-reduction methods, since it allows one to create ROMs with a computational cost which *does not* scale with the size of the full model's state vector. -->

<!-- The sample mesh is used in conjunction with what we refer to as the **stencil mesh**, which contains all cells needed to compute the velocity or residual vector on the sample mesh. Note that, in general, the sample mesh is a subset of the stencil mesh, because to compute the velocity or residual at a given cell, one also needs the cell-centered values at that target cell. Several examples of sample and stencil meshes are shown below. -->

<!-- ### 1D Example -->
<!-- The figures below show a full 1D mesh for a first order cell-centered finite volume scheme, and right below -->
<!-- a representative sample mesh (yellow-filled cells) and the stencil mesh. -->

<!-- <img src="https://github.com/Pressio/pressio-demoapps/blob/develop/figures/readme_1dmesh.png" width="75%"> -->

<!-- ### 2D Example -->
<!-- The second example shows a two dimensional sample mesh and a stencil mesh for a first order cell-centered finite volume scheme. -->
<!-- The coloring scheme is the same as in the first example. -->
<!-- <img src="https://github.com/Pressio/pressio-demoapps/blob/develop/figures/readme_2dmesh.png" width="75%"> -->

<!-- ### 3D Example -->
<!-- The third example shows a three dimensional sample mesh and a stencil mesh for a first order cell-centered finite volume scheme. The coloring scheme is the same as in the previous examples. Cell IDs are omitted for clarity. -->
<!-- <img src="https://github.com/Pressio/pressio-demoapps/blob/develop/figures/readme_3dmesh.png" width="75%"> -->
