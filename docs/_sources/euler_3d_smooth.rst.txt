3D Euler Smooth
===============

This problem solves the *3D Euler equations* in conserative form for a smooth field.

- `Reference paper <https://www.sciencedirect.com/science/article/pii/S0021999117307830>`_:
  in the original paper, they used :math:`v = -0.5`, while we use :math:`v = 1` (see below)

- IC (primitive variables):

  - :math:`\rho = 1 + 0.2*sin(\pi (x+y+z))`

  - :math:`u = 1, v = 1, w = 1, p = 1`


- Domain is ``[-1, 1]^3`` with periodic BC

- Analytical density at time `t`: :math:`\rho = 1 + 0.2*sin(\pi (x+y+z - 3 t))`

- Typically, integration is performed for :math:`t \in (0, 2.)`.


.. Caution::
   Currently, this problem only works for first order and Weno3 inviscid flux reconstruction.

Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem euler3dsmooth_s{3,5} -n Nx Ny Nz --outDir <destination-path>


where ``Nx, Ny`` are the number of cells, and ``<stencilSize> = 3 or 5``,
and ``<destination-path>`` is where you want the mesh files to be generated.

C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/euler3d.hpp"
   // ...
   namespace pda     = pressiodemoapps;
   const auto probId = pda::Euler3d::PeriodicSmooth;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   # ...
   probId  = pda.Euler3d.PeriodicSmooth
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()
