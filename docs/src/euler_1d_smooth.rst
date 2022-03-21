1D Euler Smooth
===============

This problem solves the *1D conservative Euler equations*

.. math::

   \frac{\partial }{\partial t} \begin{bmatrix}\rho \\ \rho u\\ \rho E \end{bmatrix} + \frac{\partial }{\partial x} \begin{bmatrix}\rho u \\ \rho u^2 +p\\ u(E+p) \end{bmatrix} = 0

where the pressure :math:`p` is related to the conserved quantities through the equation of the state

.. math::

   p=(\gamma -1)(\rho E-\frac{1}{2}\rho u^2)


* Initial conditions in primitive variables:

  - :math:`\rho(x, 0) = 1 + 0.2 \sin(\pi x)`

  - :math:`u(x,0) = 1`

  - :math:`p(x,0) = 1`

  - These are used to create the initial conditions in conservative variables.

- By default, :math:`\gamma = 1.4`

* Domain is :math:`[-1,1]` with periodic BC

* Analytical density as function of time :math:`t` is given as :math:`\rho(t) = 1 + 0.2\sin(\pi (x-t))`

* Typically, integration is performed over :math:`t \in (0, 2)`

* The problem is adapted from `this paper <https://www.proquest.com/openview/ef6ab9a87e7563ad18e56c2f95f624d8/1?pq-origsite=gscholar&cbl=2032364>`_



Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
          --problem euler1dsmooth_s<stencilSize> -n <N> --outDir <destination-path>

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

   const auto probId = pda::Euler1d::PeriodicSmooth;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOrder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   meshObj = pda.load_cellcentered_uniform_mesh("path-to-mesh")

   probId  = pda.Euler1d.PeriodicSmooth
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)
