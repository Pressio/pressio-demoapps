1D Euler Sod
============

This problem solves the *1D Euler conservative equations* for the Sod1d problem.

* The problem is adapted from `this paper <https://iopscience.iop.org/article/10.1086/317361>`_

- Initial conditions:

  - :math:`x<=0 :\quad \rho =1, u = 0, p = 1`

  - :math:`x>0 :\quad \rho =1/8, u = 0, p = 0.1`

- Domain is :math:`[-0.5, 0.5]` with homogeneous Neumann BC

- Typically, integration is performed over :math:`t \in (0, 0.2)`


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
	   --problem sod1d_s<stencilSize> -n <N> --outDir <destination-path>

where ``N`` is the number of cells you want and ``<stencilSize> = 3 or 5 or 7``,
and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/euler1d.hpp"

   namespace pda     = pressiodemoapps;
   const auto probId = pda::Euler1d::Sod;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   probId  = pda.Euler1d.Sod
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)


Sample Solution
---------------

Representative plot for total simulation time :math:`T=0.2` showing density at selected time steps :math:`t \in \left \{0, 0.1, 0.2\right \}`
obtained using time step :math:`dt = 10^{-4}`, Weno5, Runge-Kutta4 integration with a mesh of :math:`N=1000` cells.

.. image:: ../../figures/wiki_sod1d_0.0001_0.2_1000_weno5_rk4.png
  :width: 60 %
  :align: center
  :alt: Alternative text
