2D Euler Sedov (with symmetry)
==============================

This problem solves the *2D conservative Euler equations*

.. math::

   \frac{\partial }{\partial t} \begin{bmatrix}\rho \\ \rho u_x \\ \rho u_y\\ \rho E \end{bmatrix} + \frac{\partial }{\partial x} \begin{bmatrix}\rho u_x \\ \rho u_x^2 +p \\ \rho u_x u_y \\ (E+p)u_x \end{bmatrix} \frac{\partial }{\partial y} \begin{bmatrix}\rho u_y  \\ \rho u_x u_y \\ \rho u_y^2 +p \\ (E+p)u_y \end{bmatrix}= 0

where the pressure :math:`p` is related to the conserved quantities through the equation of the state

.. math::

   p=(\gamma -1)(\rho E-\frac{1}{2}\rho (u_x^2 + u_y^2)).


* Initial conditions in primitive variables:

  - a high pressure concentrated small spherical region of radius :math:`R = 3 \min(dx, dy)`

  - :math:`\left\{\begin{matrix}\rho =1, u = 0, v = 0, p = ((\gamma-1)0.851072)/(\pi R^2); & r\leq R \\ \rho =1, u = 0, v = 0, p = 2.5\cdot 10^{-5}; & r>R \end{matrix}\right.`

  - This IC is used to create the corresponding initial conditions in conservative variables.

- By default, :math:`\gamma = 1.4`

- Domain is :math:`[0.0, 1.2]^2` with reflective BC on :math:`y=0` and :math:`x=0` and homogeneous Neumann for :math:`x=1.2` and :math:`y=1.2`

- Typically, integration is performed for :math:`t \in (0, 1)`


* The problem is adapted from `this paper <https://www.researchgate.net/publication/260967068_GENASIS_General_Astrophysical_Simulation_System_I_Refinable_Mesh_and_Nonrelativistic_Hydrodynamics>`_


.. * The problem is adapted from `Reference paper1 <https://reader.elsevier.com/reader/sd/pii/S002199911400477X?token=658F08D28B5C7A6A97E6F4478FD494699F3C8DF23970A256F06E501B7B136F9A6A540EEA749F28AC2AF4A6A7993A8517&originRegion=eu-west-1&originCreation=20210611123033>`_ and `Reference paper2 <http://flash.uchicago.edu/site/flashcode/user_support/flash_ug_devel/node184.html#SECTION010114000000000000000>`_


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem sedov2dsym_s<stencilSize> -n Nx Ny --outDir <destination-path>

where 

- ``Nx, Ny`` is the number of cells you want along :math:`x` and :math:`y` respectively

- ``<stencilSize> = 3 or 5 or 7``: defines the neighboring connectivity of each cell 

- ``<destination-path>`` is where you want the mesh files to be generated.
  The script creates the directory if it does not exist.


.. Important::

  When you set the ``<stencilSize>``, keep in mind the following constraints (more on this below):

  - ``InviscidFluxReconstruction::FirstOrder`` requires ``<stencilSize> >= 3``
 
  - ``InviscidFluxReconstruction::Weno3`` requires ``<stencilSize> >= 5``
  
  - ``InviscidFluxReconstruction::Weno5`` requires ``<stencilSize> >= 7``



C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/euler2d.hpp"

   namespace pda     = pressiodemoapps;

   const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("path-to-mesh");

   const auto probId = pda::Euler2d::SedovSymmetry;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOrder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   meshObj = pda.load_cellcentered_uniform_mesh("path-to-mesh")

   probId  = pda.Euler2d.SedovSymmetry
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()
