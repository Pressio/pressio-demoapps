1D Linear Advection
===================

This problem solves the *1D linear advection* 

.. math::
   \frac{\partial \phi}{\partial t} + a \frac{\partial \phi}{\partial x} = 0

for a scalar field :math:`\phi` and advection velocity :math:`a`.

* Initial condition: :math:`\phi(x, 0) = \sin(\pi x)`
* Domain is :math:`[-1,1]` with periodic BC
* Integration is typically performed over :math:`t \in (0, 2k)` where :math:`k \in \mathbb{Z}`
* Default setting: :math:`a=1`

Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
	  --problem linadv1d_s<stencilSize> -n <N> --outDir <destination-path>

where:  

- ``N`` is the number of cells you want 

- ``<stencilSize> = 3 or 5 or 7``: defines the neighboring connectivity of each cell 

- ``<destination-path>``: full path to where you want the mesh files to be generated. 
  The script creates the directory if it does not exist.


.. Important::

  When you set the ``<stencilSize>``, keep in mind the following constraints (more on this below):

  - ``InviscidFluxReconstruction::FirstOder`` requires ``<stencilSize> >= 3``
 
  - ``InviscidFluxReconstruction::Weno3`` requires ``<stencilSize> >= 5``
  
  - ``InviscidFluxReconstruction::Weno5`` requires ``<stencilSize> >= 7``



C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/advection1d.hpp"

   namespace pda = pressiodemoapps;

   const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("path-to-mesh");

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

   meshObj = pda.load_cellcentered_uniform_mesh_eigen("path-to-mesh")

   # 1. using default velocity
   probId  = pda.Advection1d.PeriodicLinear
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)

   # 2. specify a custom velocity
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_linear_advection_1d_problem(meshObj, scheme, #your_vel)


Sample Solution
---------------

Representative plot showing initial condition and solution at :math:`t=1` and :math:`t=2`,
obtained using :math:`dt = 10^{-3}`, Weno5, Runge-Kutta4 integration with a mesh of :math:`N=250` cells.


.. image:: ../../figures/wiki_advection_0.001_2_250_weno5_rk4.png
  :width: 60 %
  :align: center
  :alt: Alternative text
