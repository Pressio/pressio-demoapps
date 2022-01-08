2D Shallow water equations
==========================

Simulates the 2D shallow water equations with conservative variables.

- IC: :math:`h = 1 + 1/8 \exp( -(x-1)^2 - (y-1)^2), u = v = 0`

- Domain is ``[-5,5]^2`` with slip-wall BCs

**todo: add equations**

Mesh
----

.. code-block:: bash

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem swe2dSlipWall_s{3,5,7} -n Nx Ny --outDir ...

where ``Nx, Ny`` are the number of cells, and ``<stencilSize> = 3 or 5 or 7``.


C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::Swe2d::SlipWall;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.Swe2d::SlipWall;
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()



Sample Plot
-----------

Representative *height* plot as a function of time at `x=y=0`
using a 65x65 mesh with Weno5 and RK4 time integration:

.. image:: ../../figures/wiki_2dswe_height.png
  :width: 60 %
  :alt: Alternative text
  :align: center
