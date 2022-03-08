1D Euler Lax
============

This problem solves the *1D Euler conservative equations* for the 1d Lax problem.

* `Reference paper <https://www.researchgate.net/publication/274407416_Finite_Difference_Hermite_WENO_Schemes_for_Hyperbolic_Conservation_Laws>`_

- IC:

  - :math:`x<=0`: :math:`\rho = 0.445,  u = 0.698, p = 3.528`

  - :math:`x>0`: :math:`\rho = 0.5, u = 0., p = 0.571`

- Domain is ``[-5.0, 5.0]`` with homogeneous Neumann BC

- Time integration is performed over :math:`t \in (0, 1.3)`.


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem lax1d_s<stencilSize> -n <N> --outDir <destination-path>

where ``N`` is the number of cells you want and ``<stencilSize> = 3 or 5 or 7``,
and ``<destination-path>`` is where you want the mesh files to be generated.

C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/euler1d.hpp"
   // ...
   namespace pda     = pressiodemoapps;
   const auto probId = pda::Euler1d::Lax;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   # ...
   probId  = pda.Euler1d.Lax
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)


Sample Solution
---------------

Representative plot for total simulation time ``T=1.3`` showing density at ``t=0, 0.65, 1.3``
obtained using ``dt = 0.0001``, ``Weno5``, SSPRK3 integration with a mesh of 1000 cells.

.. image:: ../../figures/wiki_lax1d_0.0001_1.3_1000_weno5_ssp3.png
  :width: 60 %
  :align: center
  :alt: Alternative text
