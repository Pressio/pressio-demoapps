2D Euler Smooth
===============

This problem solves the *2D conservative Euler equations* for a smooth field. The gas dynamics is governed by a system of PDE

.. math::

   \frac{\partial }{\partial t} \begin{bmatrix}\rho \\ \rho u_x \\ \rho u_y\\ \rho E \end{bmatrix} + \frac{\partial }{\partial x} \begin{bmatrix}\rho u_x \\ \rho u_x^2 +p \\ \rho u_x u_y \\ (E+p)u_x \end{bmatrix} \frac{\partial }{\partial y} \begin{bmatrix}\rho u_y  \\ \rho u_x u_y \\ \rho u_y^2 +p \\ (E+p)u_y \end{bmatrix}= 0

where the pressure :math:`p` is related to the conserved quantities through the equation of the state

.. math::

   p=(\gamma -1)(E-\frac{1}{2}\rho (u_x^2 + u_y^2)).


* The problem is adapted from `this paper <https://www.proquest.com/openview/ef6ab9a87e7563ad18e56c2f95f624d8/1?pq-origsite=gscholar&cbl=2032364>`_

* Initial conditions: 
  
  - :math:`\rho = 1 + 0.2\sin(\pi (x+y))`
  
  - :math:`u = 1, v = 1, p = 1`
  
* Domain is :math:`[-1, 1]^2` with periodic BC

* Analytical density as function of time is given as :math:`\rho(t) = 1 + 0.2\sin(\pi (x+y - 2 t))`

* Typically, integration is performed for :math:`t \in (0, 2)`


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem euler2dsmooth_s{3,5,7} -n Nx Ny --outDir <destination-path>

where ``Nx, Ny`` are the number of cells you want along :math:`x` and :math:`y` respectively, and ``<stencilSize> = 3 or 5 or 7``,
and ``<destination-path>`` is where you want the mesh files to be generated.

C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/euler2d.hpp"
   // ...
   namespace pda     = pressiodemoapps;
   const auto probId = pda::Euler2d::PeriodicSmooth;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   # ...
   probId  = pda.Euler2d.PeriodicSmooth
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()


Sample Plot
-----------

Representative *density* field at selected time :math:`t=2` using a ``100x100`` mesh with Weno3
and RK4 time integration:

.. image:: ../../figures/wiki_2d_smooth_density.png
  :width: 60 %
  :alt: Alternative text
  :align: center
