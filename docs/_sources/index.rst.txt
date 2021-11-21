Overview
========


pressio-demoapps is a collection of 1D, 2D and 3D problems
of varying complexity (from linear advection, to diffusion and compressible Euler)
that can be used for a variety of purposes.

Is this *another* physics PDE solver?
-------------------------------------

No! This will **not** become a full solver, since its scope is and will stay limited.

.. Important::
   The main distinguishing feature from other codes is the built-in support
   to handle a sample mesh (which, in brief, is the capability of evaluating
   the residual at a disjoint subset of the mesh cells. More on this below),
   and the focus on providing self-contained and well-defined *problems* to study.

This code, in fact, was originally started as part of the Pressio project
to create a suite of benchmark problems to test ROMs and hyper-reduction techniques,
but it is being developed to be self-contained.
For example, you can just use it for doing "standard" simulations,
or you can just use the Python meshing scripts, or leverage the sample
mesh capability to study function approximations.
One of the objectives is to provide a testbed to explore the impact
of varying numerical schemes, and test new sample mesh
and hyper-reduction techniques.

Some features of this work include:

* support for both C++ and Python

* cell-centered finite volume discretization with various numerical schemes and *exact Jacobians*

* built-in support for sample mesh in 1D, 2D and 3D for varying stencil sizes

* focus on providing self-contained and well-defined problems



Enough text, show me some code!
-------------------------------

In C++, a representative snippet showing our idea of a *problem*:

.. code-block:: c++

   #import <pressiodemoapps/euler1d.hpp>

   int main(){
     namespace pda        = pressiodemoapps;
     const auto meshObj   = pda::load_cellcentered_uniform_mesh_eigen("<path-to-mesh-files>");
     const auto scheme    = pda::InviscidFluxReconstruction::Weno5;
     const auto problemId = pda::Euler1d::Sod;
     auto problem         = pda::create_problem_eigen(meshObj, problemId, scheme);
   }


You don't want to use C++? We have Python bindings too!

.. code-block:: py

   import pressiodemoapps as pda

   mesh      = pda.load_cellcentered_uniform_mesh("<path-to-mesh-files>")
   scheme    = pda.InviscidFluxReconstruction.Weno5
   problemId = pda.Euler1d.Sod
   problem   = pda.create_problem(mesh, problemId, scheme)


Core development principles
---------------------------

* **simplicity**: we use high-level abstractions and well-defined APIs aiming
  to make this code as simple as possible to use, thus minimizing the number
  of steps needed to set up and run a problem.
  The main idea is a three-step process: (1) generate the mesh, (2) load the mesh and (3) create a problem instance that has a specific API.

* **quality assurance**: we maintain a test suite for both the C++ code and Python bindings to ensure stability and reliability;

* **extensibility**: easily add new problems, other numerical schemes (e.g. fluxes, reconstruction methods) and different data types.


..
   * **generic programming**: the core component is a header-only C++ library based on generic programming, thus allowing us to write the code once and then instantiate it with various data types;
   * **fast prototyping**: being able to do fast prototyping is essential for research. This motivated us to develop Python bindings.
     This makes the code accessible to Python users, while maintaining the performance of the compiled C++ backend;


Contents
========

.. toctree::
    :maxdepth: 1

    installation
    howtouse
    your_first_problem
    api
    problems_1d
    problems_2d
    problems_3d
    meshing
    GitHub Repo <https://github.com/Pressio/pressio-demoapps>
    Open an issue/feature req. <https://github.com/Pressio/pressio-demoapps/issues>
    license
