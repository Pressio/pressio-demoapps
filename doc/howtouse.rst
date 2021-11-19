How to use it: 3 steps!
=======================

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
   const auto meshObj   = pda::loadCellCenterUniformMeshEigen("/home/myTest");
   constexpr auto order = pda::InviscidFluxReconstruction::Weno5;
   auto problem         = pda::createProblemEigen(meshObj, pda::Euler1d::Sod, order);
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
   meshObj = pda.loadCellCenterUniformMesh("/home/myTest")
   order   = pda.InviscidFluxReconstruction.Weno5;
   problem = pda.createProblem(meshObj, pda.Euler1d.Sod, order)
   state   = problem.initialCondition()

In Python, this creates a problem based on ``numpy``.
We remark that we aimed for the Python code to be readable as similarly as possible to the C++.



Step 3: Solving the problem
---------------------------

To use a problem instance, you need to know that all *pressio-demoapps*
problem instances meet a specific API as described below.


C++ API
^^^^^^^

.. cpp:class:: Problem

   .. cpp:type:: scalar_type

      The scalar type: this is by default ``double``.

   .. cpp:type:: state_type

      Data structure type to store the state: currently, this is an Eigen vector.

   .. cpp:type:: velocity_type

      Data structure type to store the velocity: currently, this is an Eigen vector.

   .. cpp:type:: jacobian_type

      Data structure type to store the Jacobian: currently, this is an Eigen sparse Crs matrix.

   .. cpp:function:: auto totalDofSampleMesh()

      Returns the total number of degrees of freedom on the **sample** mesh.
      Note that, in general, this is not the same as the number of sample mesh cells. 
      When you have multiple dofs/cell (for example Euler equations), 
      then the total # of dofs on sample mesh = # dofs/cell times the # of sample mesh cells. 

   .. cpp:function:: auto totalDofStencilMesh()

      Returns the total number of degrees of freedom on the **stencil** mesh.
      Note that, in general, this is not the same as the number of stencil mesh cells. 
      When you have multiple dofs/cell (for example Euler equations), 
      then the total # of dofs on stencil mesh = # dofs/cell times the # of sample mesh cells. 

   .. cpp:function:: state_type initialCondition()

      Constructs and returns an instance of the initial condition for the target problem.

   .. cpp:function:: velocity_type createVelocity()

      Constructs and returns an instance of the velocity.

   .. cpp:function:: jacobian_type createJacobian()

      Constructs and returns an instance of the jacobian.

   .. cpp:function:: void velocity(const state_type & y, scalar_type time, velocity_type & v)

      Given a state :class:`y` and time :class:`time`,
      evaluates the RHS of the system and overwrites :class:`v`.

   .. cpp:function:: void jacobian(const state_type & y, scalar_type time, jacobian_type & J)

      Given a state :class:`y` and time :class:`time`,
      evaluates the Jacobian of the RHS and stores it into :class:`J`.

   .. cpp:function:: void velocityAndJacobian(const state_type & y, scalar_type time, velocity_type & v, jacobian_type & J)

      Given a state :class:`y` and time :class:`time`,
      evaluates the RHS and its Jacobian. 


Python API
^^^^^^^^^^

.. py:class:: Problem

   .. py:method:: auto totalDofSampleMesh()

      Returns the total number of degrees of freedom on the **sample** mesh.
      Note that, in general, this is not the same as the number of sample mesh cells. 
      When you have multiple dofs/cell (for example Euler equations), 
      then the total # of dofs on sample mesh = # dofs/cell times the # of sample mesh cells. 

      :rtype: integer

   .. py:method:: auto totalDofStencilMesh()

      Returns the total number of degrees of freedom on the **stencil** mesh.
      Note that, in general, this is not the same as the number of stencil mesh cells. 
      When you have multiple dofs/cell (for example Euler equations), 
      then the total # of dofs on stencil mesh = # dofs/cell times the # of sample mesh cells. 

      :rtype: integer

   .. py:method:: initialCondition()

      Constructs and returns an instance of the initial condition for the target problem.

      :rtype: numpy.array

   .. py:method:: createVelocity()

      Constructs and returns an instance of the velocity.

      :rtype: numpy.array


   .. py:method:: createApplyJacobianResult(operand)

      Constructs and returns an instance of the action of the Jacobian applied to :class:`operand`. 
      The result is constructed, and zeroed out before returning it. 

      :param numpy.array operand: rank-1 or rank-2 operand to apply the Jacobian to.
      :rtype: numpy.array


   .. py:method:: velocity(y, time, v)

      Given a state :class:`y` and time :class:`time`,
      evaluates the RHS of the system and stores it into :class:`v`.

      :param numpy.array y: state vector
      :param float time: evaluation time
      :param numpy.array v: velocity to overwrite

   .. py:method:: applyJacobianResult(y, time, operand, result)

      Given a state :class:`y` and time :class:`time`,
      this computes the action of the Jacobian applied to :class:`operand`.

      :param numpy.array y: state vector
      :param float time: evaluation time
      :param numpy.array operand: rank-1 or rank-2 operand to apply the Jacobian to.
      :rtype: numpy.array

.. note::
   Note how the Python interface only supports the Jacobian **action**.
   The main reason behind this is that Pybind11 does not yet allow view semantics for
   Python sparse matrices. 



**!! to do: finish**
