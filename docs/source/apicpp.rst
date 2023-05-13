.. role:: raw-html-m2r(raw)
   :format: html


C++ API
=======

Basic C++ API of the :ref:`cpp-problem-api` and :ref:`cpp-mesh-api`.

.. _cpp-problem-api:

Problem Class
-------------

A *pressio-demoapps* C++ problem class meets the following API

.. cpp:class:: Problem

   .. cpp:type:: scalar_type

      The scalar type: this is by default ``double``.

   .. cpp:type:: independent_variable_type

      This represents "time" and is the same as the ``scalar_type``.

   .. cpp:type:: state_type

      Data structure type storing the state: this depends on which specific implementation
      you use when you instantiate the problem using the ``create_problem_<>`` API
      (see `Supported Backends And Types`_ for more details about backends and corresponding types).
      For example, if you use ``create_problem_eigen``, the ``state_type`` will be an Eigen vector.

   .. cpp:type:: right_hand_side_type

      Data structure type to store the RHS: this depends on which specific implementation
      you use when you instantiate the problem using the ``create_problem_<>`` API
      (see `Supported Backends And Types`_ for more details about backends and corresponding types).
      For example, if you use ``create_problem_eigen``, the ``right_hand_side_type`` will be an Eigen vector.

   .. cpp:type:: rhs_type

      Same as ``right_hand_side_type``

   .. cpp:type:: jacobian_type

      Data structure type to store the Jacobian: this depends on which specific implementation
      you use when you instantiate the problem using the ``create_problem_<>`` API
      (see `Supported Backends And Types`_ for more details about backends and corresponding types).
      For example, if you use ``create_problem_eigen``, the ``jacobian_type`` will be an Eigen CRS matrix.

   .. cpp:function:: state_type initialCondition()

      Constructs and returns an instance of the initial condition for the target problem.

   .. cpp:function:: right_hand_side_type createRightHandSide()

      Constructs and returns an instance of the right hand side.

   .. cpp:function:: right_hand_side_type createRhs()

      Shortcut for ``createRightHandSide()``.

   .. cpp:function:: jacobian_type createJacobian()

      Constructs and returns an instance of the jacobian.

   .. cpp:function:: void operator()(const state_type & y, scalar_type time, right_hand_side_type & v)

      Given a state :math:`y` and time :math:`time`,
      evaluates the RHS of the system overwriting :math:`v`.

   .. cpp:function:: void operator()(const state_type & y, scalar_type time,\
		     right_hand_side_type & v, jacobian_type & J, bool computeJac)

      Given a state :math:`y` and time :math:`time`,
      evaluates the RHS of the system overwriting :math:`v` and  if ``computeJac == true``,
      its Jacobian and stores it into :math:`J`.

   .. cpp:function:: void rhsAndJacobian(const state_type & y, scalar_type time,\
       rhs_type & v, std::optional<jacobian_type *> J)

      Given a state :math:`y` and time :math:`time`,
      evaluates the RHS of the system overwriting :math:`v` and if ``J == true``,
      its Jacobian and stores it into :math:`*J.value()`.

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
     - ``using state_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>`` :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>` ``using right_hand_side_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>`` :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>` ``using jacobian_type = Eigen::SparseMatrix<scalar_type, Eigen::RowMajor, int32_t>;``


.. _cpp-mesh-api:

Cell-Centered Uniform Mesh Class
--------------------------------

A *pressio-demoapps* C++ cell-centered mesh class meets the following API

.. cpp:class:: CellCenteredUniformMesh

   .. cpp:type:: scalar_type

      This is by default ``double``.

   .. cpp:type:: index_type

      This is the type used for all integers for indexing (e.g. cell IDs, etc) so
      basically for all ordinals. Defaults to ``int32_t``, which is big enough for most current cases.

   .. cpp:function:: int dimensionality() const

      Returns the dimensionality of this mesh: returns 1 for 1d problem, 2 for 2d, etc.

   .. cpp:function:: int stencilSize() const

      Returns the size of the stencil (connectivity) of this mesh object.

   .. cpp:function:: index_type stencilMeshSize() const

      Returns the number of *stencil* cells in the mesh.
      This corresponds to all cells where the state is defined.

   .. cpp:function:: index_type sampleMeshSize() const

      Returns the number of *sample* cells in the mesh.
      This corresponds to all cells where the RHS is defined.

   .. cpp:function:: scalar_type dx() const

      Returns the cell size along the x axis.

   .. cpp:function:: scalar_type dy() const

      Returns the cell size along the y axis. This is applicable only
      if the dimensionality is >= 2.

   .. cpp:function:: scalar_type dz() const

      Returns the cell size along the z axis. This is applicable only
      if the dimensionality is == 3.

   .. cpp:function:: auto viewX() const

      Returns a *reference* to the vector of x-coordinates of all *stencil* mesh cells.
      The type of container returned by reference depends on the backend used.

   .. cpp:function:: auto viewY() const

      Returns a *reference* to the vector of y-coordinates of all *stencil* mesh cells
      The type of container returned by reference depends on the backend used.

   .. cpp:function:: auto viewZ() const

      Returns a *reference* to the vector of z-coordinates of all *stencil* mesh cells
      The type of container returned by reference depends on the backend used.
