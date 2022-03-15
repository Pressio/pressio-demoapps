3D Euler Sedov (with symmetry)
==============================

This problem solves the *3D Euler equations* in conservative form. The gas dynamics is governed by a system of PDE

.. math::

   \frac{\partial }{\partial t} \begin{bmatrix}\rho \\ \rho u_x \\ \rho u_y \\ \rho u_z\\ \rho E \end{bmatrix} + \frac{\partial }{\partial x} \begin{bmatrix}\rho u_x \\ \rho u_x^2 +p \\ \rho u_x u_y  \\ \rho u_x u_z\\ (\rho E+p)u_x \end{bmatrix} + \frac{\partial }{\partial y} \begin{bmatrix}\rho u_y  \\ \rho u_x u_y \\ \rho u_y^2 +p \\ \rho u_y u_z\\ (\rho E+p)u_y \end{bmatrix} + \frac{\partial }{\partial z} \begin{bmatrix}\rho u_z  \\ \rho u_x u_z  \\ \rho u_y u_z \\ \rho u_z^2 +p\\ (\rho E+p)u_z \end{bmatrix}= 0

where the pressure :math:`p` is related to the conserved quantities through the equation of the state

.. math::

   p=(\gamma -1)(\rho E-\frac{1}{2}\rho (u_x^2 + u_y^2 + u_z^2)).


- Initial conditions (primitive variables):
    
  - a high pressure concentrated small spherical region of radius :math:`R = 3 \min(dx, dy)`

  - :math:`\left\{\begin{matrix}\rho =1, u = 0, v = 0, p = ((\gamma-1)0.851072)/(4 \pi R^3); & r\leq R \\ \rho =1, u = 0, v = 0, p = 2.5\cdot 10^{-5}; & r>R \end{matrix}\right.`

  - These are used to create the initial conditions in conservative variables.

- By default, :math:`\gamma = 1.4`

- Domain is :math:`[0.0, 1.2]^3` with reflective BC on :math:`x=y=z=0` and homogeneous Neumann on others.

- Typically, integration is performed over :math:`t \in (0, 1.)`.

.. Caution::
   Currently, this problem only works for first order and Weno3 inviscid flux reconstruction.

Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem sedov3dsym_s<stencilSize> -n Nx Ny Nz --outDir <destination-path>

where 

- ``Nx, Ny, Nz`` is the number of cells you want along :math:`x`, :math:`y`, :math:`z` respectively

- ``<stencilSize> = 3 or 5``

- ``<destination-path>`` is where you want the mesh files to be generated


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

   const auto probId = pda::Euler3d::SedovSymmetry;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOrder; //or Weno3
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   meshObj = pda.load_cellcentered_uniform_mesh_eigen("path-to-mesh")

   probId  = pda.Euler3d.SedovSymmetry
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()
