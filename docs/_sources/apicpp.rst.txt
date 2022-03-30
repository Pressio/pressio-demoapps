.. role:: raw-html-m2r(raw)
   :format: html


Problem class: C++ API
======================

A *pressio-demoapps* C++ problem class meets a specific API as described below.


C++ API
-------

.. cpp:class:: Problem

   .. cpp:type:: scalar_type

      The scalar type: this is by default ``double``.

   .. cpp:type:: state_type

      Data structure type storing the state: this depends on which specific implementation
      you use when you instantiate the problem using the ``create_problem_<>`` API
      (see `Supported Backends And Types`_ for more details about backends and corresponding types).
      For example, if you use ``create_problem_eigen``, the ``state_type`` will be an Eigen vector.

   .. cpp:type:: velocity_type

      Data structure type to store the velocity (or RHS): this depends on which specific implementation
      you use when you instantiate the problem using the ``create_problem_<>`` API
      (see `Supported Backends And Types`_ for more details about backends and corresponding types).
      For example, if you use ``create_problem_eigen``, the ``velocity_type`` will be an Eigen vector.

   .. cpp:type:: jacobian_type

      Data structure type to store the Jacobian: this depends on which specific implementation
      you use when you instantiate the problem using the ``create_problem_<>`` API
      (see `Supported Backends And Types`_ for more details about backends and corresponding types).
      For example, if you use ``create_problem_eigen``, the ``jacobian_type`` will be an Eigen CRS matrix.

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

   .. cpp:function:: void velocityAndJacobian(const state_type & y, \
		                              scalar_type time, \
					      velocity_type & v, \
					      jacobian_type & J)

      Given a state :class:`y` and time :class:`time`,
      evaluates the RHS and its Jacobian.



.. _Supported Backends And Types:

C++ Backends and Corresponding Types
------------------------------------

.. list-table::
   :widths: 5 95
   :header-rows: 1
   :align: left

   * - Backend
     - Type alias

   * - Eigen
     - ``using state_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>`` :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>` ``using velocity_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>`` :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>` ``using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int32_t>;``
