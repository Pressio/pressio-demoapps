How to use it in 3 steps!
=========================

To use *any* of the problems supported, you should always keep
in mind that there are three steps involved: *generate* the mesh,
*instantiate* the problem, and *use the problem object*.

For demonstration purposes, let's say we are interested in the Euler1d Sod problem.
We show below how to use it but similar steps are applicable to any other problem.
So after reading this page, you should be able to understand how to use any problem.

Step 1: Generating the Mesh
---------------------------

To generate the mesh, you should rely on the meshing script as follows

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
	   --outdir $HOME/myTestMesh \
	   --problem sod1d_s7 \
	   -n 100

where we specify the things:

- ``--outdir``: full path to where you want all the mesh files to be generated.
  Here we use ``$HOME/myTestMesh`` but you can obviously set whatever you want.

- ``--problem``: we pass the string ``sod1d_s7`` which has two parts:
  ``sod1d`` indicates the problem name, and ``s7`` indicates
  the stencil size to use, in this case we want a 7-point stencil.
  So if you want a mesh for Sod1d with 3-point stencil then you do ``sod1d_s3``.
  Refer to each problem page for the supported choices.

- ``-n=100`` specifies how many cells to use for the spatial discretization

The advantage of this script is that any information about the problem domain
and other details are encoded in the script, so it only exposes a minimal set of parameters
(e.g. number of cells) to set.

.. Tip::

   For the correct syntax to create the mesh for any problem, refer to each problem page.

|


Step 2: Creating a problem instance
-----------------------------------

After generating a mesh, we need to create an instance of the Sod1d problem.
We now show what this step 2 looks like in both C++ and Python.

C++ Synopsis to create a problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: c++

   #import <pressiodemoapps/euler1d.hpp>
   // ...
   namespace pda	 = pressiodemoapps;
   const auto meshObj    = pda::load_cellcentered_uniform_mesh_eigen("/home/myTestMesh");
   constexpr auto scheme = pda::InviscidFluxReconstruction::Weno5;
   auto problem          = pda::create_problem_eigen(meshObj, pda::Euler1d::Sod, scheme);
   auto state		 = problem.initialCondition();
   // ...

This creates an instance of the Sod1d problem using Eigen data types, and selects
the 5-th order WENO scheme for the inviscid flux reconstruction (to know more details
about which schemes and stencils a problem supports, see the target problem's webpage).
Note that here we explicitly ask for an instance
of the problem that is based on Eigen data types. Eigen is currently the main backend
supported, but other ones, e.g., Kokkos, will be added later.


Python Synopsis to create a problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: py

   import pressiodemoapps as pda
   # ...
   meshObj = pda.load_cellcentered_uniform_mesh("/home/myTestMesh")
   order   = pda.InviscidFluxReconstruction.Weno5;
   problem = pda.create_problem(meshObj, pda.Euler1d.Sod, order)
   state   = problem.initialCondition()

This creates a problem that uses ``numpy`` data structures.
Note how the Python code reads like the C++ one.


|

Step 3: Using the problem
-------------------------

To use a problem instance, you need to know that all *pressio-demoapps*
problem instances meet a specific API as described in :doc:`this page <api>`.

Assuming you read that API page, then it should be clear that the
problem API is complete enough so that you can query all operators
and do something with them.
Here we show some things you can do using the C++ as an example:

.. code-block:: c++

   // ...
   auto problem          = pda::create_problem_eigen(meshObj, pda::Euler1d::Sod, scheme);
   auto state		 = problem.initialCondition();

   // now that we have the problem and initial condition we can
   // create the RHS
   auto rhs = problem.createVelocity()
   // compute the rhs of the discrete system at time=0.0
   problem.velocity(state, 0.0, rhs);

   // create the Jacobian
   auto J = problem.createJacobian()
   // compute J at time=0.0
   problem.velocityAndJacobian(state, 0.0, rhs, J);
   // or we can compute just the Jacobian
   problem.jacobian(state, 0.0, J);



**!! to do: finish**
