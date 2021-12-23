1D Linear Advection
===================

This problem solves the *1D linear advection*.

* IC: :math:`\phi(x, 0) = \sin(\pi x)`
* Domain is [-1,1] with periodic BC
* Typically, integration is performed over :math:`t \in (0, 2k)` where k is 1 or 2 or 3 ...

Mesh
----

.. code-block:: shell

   cd pressio-demoapps/meshing_scripts
   python create_full_mesh_for.py --problem linadv1d_s<stencilSize> -n <N> --outDir <destination-path>

where ``N`` is the number of cells you want and ``<stencilSize> = 3 or 5 or 7``,
and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::Advection1d::PeriodicLinear;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.Advection1d.PeriodicLinear
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()
