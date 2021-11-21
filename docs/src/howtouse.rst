How to use it in 3 steps!
=========================

To use *any* of the problems supported, you need three steps:
generate the mesh, instantiate the problem object, and solve.

For demonstration purposes, let's say we are interested
in exploring the Euler1d Sod problem.
We show below to handle this problem but the same steps are applicable to any other.

Step 1: Generating the Mesh
---------------------------

We provide a meshing script that works as follows

.. code-block:: shell

   # assuming we are inside the pressio-demoapps dir
   python ./meshing_scripts/create_full_mesh_for.py \
		--problem sod1d_s7 \
		--outdir ${HOME}/myTest \
		-n 100

where the string ``sod1d_s7`` is composed of two parts: ``sod1d`` indicates the target problem,
and ``s7`` indicates the stencil size to use, in this case we want a 7-point stencil
(the supported choices are discussed later); ``${HOME}/myTest`` is where all the mesh files
are generated, and ``-n=100`` specifies how many cells we want.

The advantage of this script is that any information about the problem domain
and other details are encoded in the script, so it only exposes the minimal set of parameters
(e.g. number of cells) to set.


Step 2: Creating a problem instance
-----------------------------------

This second step uses the mesh generated above and creates an instance of the Sod1d problem.<br>
We now show what this step 2 looks like in C++ and Python.

C++ Synopsis
^^^^^^^^^^^^

.. code-block:: c++

   #import <pressiodemoapps/euler1d.hpp>
   // ...
   namespace pda = pressiodemoapps;
   const auto meshObj   = pda::load_cellcentered_uniform_mesh_eigen("/home/myTest");
   constexpr auto order = pda::InviscidFluxReconstruction::Weno5;
   auto problem         = pda::create_problem_eigen(meshObj, pda::Euler1d::Sod, order);
   auto state		= problem.initialCondition();
   // ...

This creates an instance of the Sod1d problem using Eigen data types, and selects
the 5-th order WENO for the inviscid flux edge reconstruction (more details on the schemes
and stencils are given below). Note that here we explicitly ask for an instance
of the problem that is based on Eigen data types. Eigen is currently the main backend
supported, but other ones, e.g., Kokkos, will be added later.


Python Synopsis
^^^^^^^^^^^^^^^

.. code-block:: py

   import pressiodemoapps as pda
   # ...
   meshObj = pda.load_cellcentered_uniform_mesh("/home/myTest")
   order   = pda.InviscidFluxReconstruction.Weno5;
   problem = pda.create_problem(meshObj, pda.Euler1d.Sod, order)
   state   = problem.initialCondition()

In Python, this creates a problem based on ``numpy``.
We remark that we aimed for the Python code to be readable as similarly as possible to the C++.



Step 3: Solving the problem
---------------------------

To use a problem instance, you need to know that all *pressio-demoapps*
problem instances meet a specific API as described in :doc:`this page <api>`.


**!! to do: finish**
