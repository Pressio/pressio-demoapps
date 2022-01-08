2D Euler Smooth
===============

This problem solves the *2D Euler equations* for a smooth field.

- IC: :math:`\rho = 1 + 0.2*sin(\pi (x+y)), u = 1, v = 1, p = 1`

- Domain is ``[-1, 1]^2`` with periodic BC

- Analytical density at time `t`: :math:`\rho = 1 + 0.2*sin(\pi (x+y - 2 t))`

- Typically, integration is performed for :math:`t \in (0, 2.)`.


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem euler2dsmooth_s{3,5,7} -n Nx Ny --outDir <destination-path>

where ``Nx, Ny`` are the number of cells, and ``<stencilSize> = 3 or 5 or 7``,
and ``<destination-path>`` is where you want the mesh files to be generated.

C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::Euler2d::PeriodicSmooth;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.Euler2d.PeriodicSmooth
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()


Sample Plot
-----------

Representative *density* field at ``t=2`` using a 100x100 mesh with Weno3 and RK4 time integration:

.. image:: ../../figures/wiki_2d_smooth_density.png
  :width: 60 %
  :alt: Alternative text
  :align: center
