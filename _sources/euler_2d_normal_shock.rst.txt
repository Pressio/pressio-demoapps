2D Euler Normal Shock
=====================

This problem solves the *2D conservative Euler equations*

.. math::

   \frac{\partial }{\partial t} \begin{bmatrix}\rho \\ \rho u_x \\ \rho u_y\\ \rho E \end{bmatrix} + \frac{\partial }{\partial x} \begin{bmatrix}\rho u_x \\ \rho u_x^2 +p \\ \rho u_x u_y \\ (E+p)u_x \end{bmatrix} \frac{\partial }{\partial y} \begin{bmatrix}\rho u_y  \\ \rho u_x u_y \\ \rho u_y^2 +p \\ (E+p)u_y \end{bmatrix}= 0

where the pressure :math:`p` is related to the conserved quantities through the equation of the state

.. math::

   p=(\gamma -1)(\rho E-\frac{1}{2}\rho (u_x^2 + u_y^2)).


- Initial condition is a Mach 9 shock *normal* to the x-axis positioned at :math:`x = 1/6`.

- By default, :math:`\gamma = 1.4`

- Domain is :math:`[0, 2]\times[0, 1]` 

- BC are homogeneous Neumann at left (:math:`x=0`) and right (:math:`x=2`) boundaries, 
   and reflective at top (:math:`y=1`) and bottom (:math:`y=0`) boundaries


Mesh
----


.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem normalshock2d_s<stencilSize> -n Nx Ny --outDir <destination-path>

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

   const auto probId = pda::Euler2d::NormalShock;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOrder; //or Weno3 or Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state       = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   meshObj = pda.load_cellcentered_uniform_mesh("path-to-mesh")

   probId  = pda.Euler2d.NormalShock
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3 or Weno5
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()

