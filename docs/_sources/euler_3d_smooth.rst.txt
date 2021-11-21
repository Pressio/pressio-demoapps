3D Euler Smooth
===============

This problem solves the *3D Euler equations* for a smooth field.

- `Reference source <https://www.sciencedirect.com/science/article/pii/S0021999117307830>`_:
  note that in the original paper, they used v = -0.5, while we use v = 1 (see IC below)

- IC: :math:`\rho = 1 + 0.2*sin(\pi (x+y+z)), u = 1, v = 1, w = 1, p = 1`

- Domain is ``[-1, 1]^3`` with periodic BC

- Analytical density at time `t`: :math:`\rho = 1 + 0.2*sin(\pi (x+y+z - 3 t))`

- Typically, integration is performed for :math:`t \in (0, 2.)`.


.. Caution::
   Currently, this problem only works for first order and Weno3 inviscid flux reconstruction.

Mesh
----

.. code-block:: shell

   python create_full_mesh_for.py --problem euler3dsmooth_s{3,5} -n Nx Ny Nz --outDir ...


where ``Nx, Ny`` are the number of cells, and ``<stencilSize> = 3 or 5``.


C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::Euler3d::PeriodicSmooth;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.Euler3d.PeriodicSmooth
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()
