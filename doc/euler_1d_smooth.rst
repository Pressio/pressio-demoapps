1D Euler Smooth
===============

This problem solves the *1D Euler equations* for a smooth field.

* `Reference paper <https://www.proquest.com/openview/ef6ab9a87e7563ad18e56c2f95f624d8/1?pq-origsite=gscholar&cbl=2032364>`_

* IC:

  - :math:`\rho(x, 0) = 1 + 0.2 \sin(\pi x)`

  - :math:`u(x,0) = 1`

  - :math:`p(x,0) = 1`

* Domain is ``[-1,1]`` with periodic BC

* Analytical density at time `t`: :math:`\rho = 1 + 0.2 \sin(\pi (x-t))`

* Typically, integration is performed over :math:`t \in (0, 2)`.


Mesh
----

.. code-block:: shell

   python create_full_mesh_for.py --problem euler1dsmooth_s<stencilSize> -n <N> --outDir ...

where ``N`` is the number of cells you want and ``<stencilSize> = 3 or 5 or 7``.

C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::Euler1d::PeriodicSmooth;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::createProblemEigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.Euler1d.PeriodicSmooth
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.createProblem(meshObj, probId, scheme)
   state   = problem.initialCondition()
