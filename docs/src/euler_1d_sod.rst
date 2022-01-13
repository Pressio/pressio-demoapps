1D Euler Sod
============

This problem solves the *1D Euler conservative equations* for the Sod1d problem.

* `Reference paper <https://iopscience.iop.org/article/10.1086/317361>`_

- IC:

  - :math:`x<=0`: :math:`\rho =1, u = 0, p = 1`

  - :math:`x>0`: :math:`\rho =0.125, u = 0, p = 0.1`

- Domain is ``[-0.5, 0.5]`` with homogeneous Neumann BC

- Typically, integration is performed over :math:`t \in (0, 0.2)`.


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
   // ...
   namespace pda     = pressiodemoapps;
   const auto probId = pda::Euler1d::Sod;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state        = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   # ...
   probId  = pda.Euler1d.Sod
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()


Sample Solution
---------------

Representative plot for total simulation time ``T=0.2`` showing density at ``t=0, 0.1, 0.2``
obtained using ``dt = 0.0001``, ``Weno5``, Runge-Kutta4 integration with a mesh of 1000 cells.

.. image:: ../../figures/wiki_sod1d_0.0001_0.2_1000_weno5_rk4.png
  :width: 60 %
  :align: center
  :alt: Alternative text
