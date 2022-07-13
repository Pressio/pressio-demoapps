
Python API
==========

Basic Python API of the :ref:`py-problem-api` and :ref:`py-mesh-api`.


.. _py-problem-api:

Problem Class
-------------

.. py:class:: Problem

   .. py:method:: initialCondition()

      Constructs and returns an instance of the initial condition for the target problem.

      :rtype: numpy array

   .. py:method:: createVelocity()

      Constructs and returns an instance of the velocity.

      :rtype: numpy array

   .. py:method:: createApplyJacobianResult(operand)

      Constructs and returns an instance of the action of the Jacobian applied to the `operand`.
      The result is constructed, and zeroed out before returning it.

      :param numpy array operand: rank-1 or rank-2 operand to apply the Jacobian to.
      :rtype: numpy array


   .. py:method:: velocity(y, time, v)

      Given a state :math:`y` and time :math:`time`,
      evaluates the RHS of the system and stores it into :math:`v`.

      :param numpy array y: state vector
      :param float time: evaluation time
      :param numpy array v: velocity to overwrite

   .. py:method:: applyJacobian(y, operand, time, result)

      Given a state :math:`y` and time :math:`time`,
      this computes the action of the Jacobian applied to :math:`operand`.

      :param numpy array y: state vector
      :param numpy array operand: rank-1 or rank-2 operand to apply the Jacobian to.
      :param float time: evaluation time
      :param numpy array result: rank-1 or rank-2 operand to apply the Jacobian to.

   .. py:method:: totalDofSampleMesh()

      Returns the total number of degrees of freedom on the **sample** mesh.
      Note that, in general, this is not the same as the number of sample mesh cells.
      When you have multiple dofs/cell (for example Euler equations),
      then the total # of dofs on sample mesh = # dofs/cell times the # of sample mesh cells.

      :rtype: int

   .. py:method:: totalDofStencilMesh()

      Returns the total number of degrees of freedom on the **stencil** mesh.
      Note that, in general, this is not the same as the number of stencil mesh cells.
      When you have multiple dofs/cell (for example Euler equations),
      then the total # of dofs on stencil mesh = # dofs/cell times the # of sample mesh cells.

      :rtype: int

.. note::
   Note how the Python interface only supports the Jacobian **action**.
   The main reason behind this is that Pybind11 does not yet allow view semantics for
   Python sparse matrices.


.. _py-mesh-api:

Cell-Centered Uniform Mesh Class
--------------------------------

A *pressio-demoapps* Python cell-centered mesh class meets the following API

.. py:class:: CellCenteredUniformMesh

   .. py:method:: dimensionality()

      Returns the dimensionality of this mesh: returns 1 for 1d problem, 2 for 2d, etc.

      :rtype: int

   .. py:method:: stencilSize()

      Returns the size of the stencil (connectivity) of this mesh object.

      :rtype: int

   .. py:method:: stencilMeshSize()

      Returns the number of *stencil* cells in the mesh.
      This corresponds to all cells where the state is defined.

      :rtype: int

   .. py:method:: sampleMeshSize()

      Returns the number of *sample* cells in the mesh.
      This corresponds to all cells where the velocity (or RHS) is defined.

      :rtype: int

   .. py:method:: dx()

      Returns the cell size along the x axis.

      :rtype: float

   .. py:method:: dy()

      Returns the cell size along the y axis. This is applicable only
      if the dimensionality is >= 2.

      :rtype: float

   .. py:method:: dz()

      Returns the cell size along the z axis. This is applicable only
      if the dimensionality is == 3.

      :rtype: float

   .. py:method:: viewX()

      Returns a *reference* to the vector of x-coordinates of all *stencil* mesh cells.
      The type of container returned by reference depends on the backend used.

      :rtype: numpy array

   .. py:method:: viewY()

      Returns a *reference* to the vector of y-coordinates of all *stencil* mesh cells
      The type of container returned by reference depends on the backend used.

      :rtype: numpy array

   .. py:method:: viewZ()

      Returns a *reference* to the vector of z-coordinates of all *stencil* mesh cells
      The type of container returned by reference depends on the backend used.

      :rtype: numpy array
