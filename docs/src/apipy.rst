Problem class: Python API
=========================

A *pressio-demoapps* Python problem class meets a specific API as described below.


Python API
----------

.. py:class:: Problem

   .. py:method:: totalDofSampleMesh()

      Returns the total number of degrees of freedom on the **sample** mesh.
      Note that, in general, this is not the same as the number of sample mesh cells.
      When you have multiple dofs/cell (for example Euler equations),
      then the total # of dofs on sample mesh = # dofs/cell times the # of sample mesh cells.

      :rtype: integer

   .. py:method:: totalDofStencilMesh()

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

   .. py:method:: applyJacobian(y, operand, time, result)

      Given a state :class:`y` and time :class:`time`,
      this computes the action of the Jacobian applied to :class:`operand`.

      :param numpy.array y: state vector
      :param numpy.array operand: rank-1 or rank-2 operand to apply the Jacobian to.
      :param float time: evaluation time
      :rtype: numpy.array

.. note::
   Note how the Python interface only supports the Jacobian **action**.
   The main reason behind this is that Pybind11 does not yet allow view semantics for
   Python sparse matrices.
