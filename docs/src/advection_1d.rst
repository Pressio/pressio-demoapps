1D Linear Advection
===================

This problem solves the *1D linear advection* with variable advecting velocity.

* Initial condition: :math:`y(x, 0) = \sin(\pi x)`
* Domain is :math:`-1 \leq x \leq 1` with periodic BC
* Integration is typically performed over :math:`t \in (0, 2k)` where k is 1 or 2 or 3 ...
* Default: use a unitary velocity

Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
	  --problem linadv1d_s<stencilSize> -n <N> --outDir <destination-path>

where ``N`` is the number of cells you want and ``<stencilSize> = 3 or 5 or 7``,
and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/advection1d.hpp"

   namespace pda     = pressiodemoapps;

   // 1. using default velocity
   const auto probId = pda::Advection1d::PeriodicLinear;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);

   // 2. specify a custom velocity
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_linear_advection_1d_problem_eigen(meshObj, scheme, /*vel*);


Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   # 1. using default velocity
   probId  = pda.Advection1d.PeriodicLinear
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)

   # 2. specify a custom velocity
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_linear_advection_1d_problem(meshObj, scheme, #your_vel)


Sample Solution
---------------

Representative plot showing initial condition and solution at `t=1.` and `t=2`,
obtained with default velocity, using `dt = 0.001`, Weno5, Runge-Kutta4 integration
with a mesh of 250 cells.

.. image:: ../../figures/wiki_advection_0.001_2_250_weno5_rk4.png
  :width: 60 %
  :align: center
  :alt: Alternative text
