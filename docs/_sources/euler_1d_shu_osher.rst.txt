1D Euler Shu-Osher
==================

This problem solves the *1D Euler conservative equations* for the Shu-Osher 1D problem.


.. math::

   \frac{\partial }{\partial t} \begin{bmatrix}\rho \\ \rho u\\ \rho E \end{bmatrix} + \frac{\partial }{\partial x} \begin{bmatrix}\rho u \\ \rho u^2 +p\\ u(E+p) \end{bmatrix} = 0

where the pressure :math:`p` is related to the conserved quantities through the equation of the state

.. math::

   p=(\gamma -1)(\rho E-\frac{1}{2}\rho u^2).


- Initial conditions in primivite variables:

  - :math:`x\leq-4: \quad \rho =27/7, u = 2.629369, p = 31/3`

  - :math:`x>-4: \quad \rho =1 + 0.2\sin(5x), u = 0, p = 1`

  - These are used to create the initial conditions in conservative variables.

- By default, :math:`\gamma = 1.4`

- Domain is :math:`[-5, 5]` with homogeneous Neumann BC

- Typically, integration is performed over :math:`t \in (0, 1.8)`.

* The problem is adapted from `this paper <https://www.researchgate.net/publication/226065267_Numerical_simulations_of_compressible_mixing_layers_with_a_discontinuous_Galerkin_method>`_


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
	   --problem shuOsher1d_s<stencilSize> -n <N> --outDir <destination-path>

where 

- ``N`` is the number of cells you want 

- ``<stencilSize> = 3 or 5 or 7``: defines the neighboring connectivity of each cell 

- ``<destination-path>``: full path to where you want the mesh files to be generated. 
  The script creates the directory if it does not exist.


.. Important::

  When you set the ``<stencilSize>``, keep in mind the following constraints (more on this below):

  - ``InviscidFluxReconstruction::FirstOrder`` requires ``<stencilSize> >= 3``
 
  - ``InviscidFluxReconstruction::Weno3`` requires ``<stencilSize> >= 5``
  
  - ``InviscidFluxReconstruction::Weno5`` requires ``<stencilSize> >= 7``


C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/euler1d.hpp"

   namespace pda = pressiodemoapps;

   const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("path-to-mesh");

   const auto probId = pda::Euler1d::ShuOsher;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOrder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   meshObj = pda.load_cellcentered_uniform_mesh("path-to-mesh")

   probId  = pda.Euler1d.ShuOsher
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)


Sample Solution
---------------

Representative plot for total simulation time :math:`T=1.8` showing density at selected time steps :math:`t\in \left \{ 0, 0.9, 1.8 \right \}``
obtained using :math:`dt = 10^{-3}`, Weno5, SSPRK3 integration with a mesh of :math:`N=500` cells.

.. image:: ../../figures/wiki_shuosher1d_0.001_1.8_500_weno5_ssprk3.png
  :width: 60 %
  :align: center
  :alt: Alternative text
