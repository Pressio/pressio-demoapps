1D Euler Shu-Osher
==================

This problem solves the *1D Euler conservative equations* for the Shu-Osher 1d problem.

* `Reference paper <https://www.researchgate.net/publication/226065267_Numerical_simulations_of_compressible_mixing_layers_with_a_discontinuous_Galerkin_method>`_

- IC:

  - :math:`x<=-4`: :math:`\rho =27/7, u = 2.629369, p = 31/3`

  - :math:`x>-4`: :math:`\rho =1 + 0.2*sin(5x), u = 0, p = 1`

- Domain is ``[-5, 5]`` with homogeneous Neumann BC

- Typically, integration is performed over :math:`t \in (0, 1.8)`.

Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
	   --problem shuOsher1d_s<stencilSize> -n <N> --outDir <destination-path>

where ``N`` is the number of cells you want and ``<stencilSize> = 3 or 5 or 7``,
and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/euler1d.hpp"

   namespace pda     = pressiodemoapps;
   const auto probId = pda::Euler1d::ShuOsher;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   probId  = pda.Euler1d.ShuOsher
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)


Sample Solution
---------------

Representative plot for total simulation time ``T=1.8`` showing density at ``t=0, 0.9, 1.8``
obtained using ``dt = 0.001``, ``Weno5``, SSPRK3 integration with a mesh of 500 cells.

.. image:: ../../figures/wiki_shuosher1d_0.001_1.8_500_weno5_ssprk3.png
  :width: 60 %
  :align: center
  :alt: Alternative text
