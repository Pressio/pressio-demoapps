3D Euler Smooth
===============

This problem solves the *3D Euler equations* in conservative form for a smooth field. The gas dynamics is governed by a system of PDE

.. math::

   \frac{\partial }{\partial t} \begin{bmatrix}\rho \\ \rho u_x \\ \rho u_y \\ \rho u_z\\ \rho E \end{bmatrix} + \frac{\partial }{\partial x} \begin{bmatrix}\rho u_x \\ \rho u_x^2 +p \\ \rho u_x u_y  \\ \rho u_x u_z\\ (\rho E+p)u_x \end{bmatrix} + \frac{\partial }{\partial y} \begin{bmatrix}\rho u_y  \\ \rho u_x u_y \\ \rho u_y^2 +p \\ \rho u_y u_z\\ (\rho E+p)u_y \end{bmatrix} + \frac{\partial }{\partial z} \begin{bmatrix}\rho u_z  \\ \rho u_x u_z  \\ \rho u_y u_z \\ \rho u_z^2 +p\\ (\rho E+p)u_z \end{bmatrix}= 0

where the pressure :math:`p` is related to the conserved quantities through the equation of the state

.. math::

   p=(\gamma -1)(\rho E-\frac{1}{2}\rho (u_x^2 + u_y^2 + u_z^2)).


- Initial conditions (primitive variables):

  - :math:`\rho = 1 + 0.2\sin(\pi (x+y+z))`

  - :math:`u = 1, v = 1, w = 1, p = 1`

  - These are used to create the initial conditions in conservative variables.

- By default, :math:`\gamma = 1.4`

- Domain is :math:`[-1, 1]^3` with periodic BC

- Analytical density as function of time is given as: :math:`\rho(t) = 1 + 0.2\sin(\pi (x+y+z - 3 t))`

- Typically, integration is performed for :math:`t \in (0, 2)`

- The problem is adopted from `this paper <https://www.sciencedirect.com/science/article/pii/S0021999117307830>`_:
  
  - in the original paper, they used :math:`v = -0.5`, while we use :math:`v = 1`


.. Caution::
   Currently, this problem only works for first order and Weno3 inviscid flux reconstruction.

Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem euler3dsmooth_s<stencilSize> -n Nx Ny Nz --outDir <destination-path>

where 

- ``Nx, Ny, Nz`` is the number of cells you want along :math:`x`, :math:`y`, :math:`z` respectively

- ``<stencilSize> = 3 or 5``: defines the neighboring connectivity of each cell 

- ``<destination-path>``: full path to where you want the mesh files to be generated. 
  The script creates the directory if it does not exist.


.. Important::

  When you set the ``<stencilSize>``, keep in mind the following constraints (more on this below):

  - ``InviscidFluxReconstruction::FirstOrder`` requires ``<stencilSize> >= 3``
 
  - ``InviscidFluxReconstruction::Weno3`` requires ``<stencilSize> >= 5``


C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/euler3d.hpp"

   namespace pda = pressiodemoapps;

   const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("path-to-mesh");

   const auto probId = pda::Euler3d::PeriodicSmooth;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOrder; //or Weno3
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   meshObj = pda.load_cellcentered_uniform_mesh("path-to-mesh")

   probId  = pda.Euler3d.PeriodicSmooth
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()
