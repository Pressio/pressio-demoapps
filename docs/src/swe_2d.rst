2D Shallow water equations
==========================

Simulates the 2D shallow water equations with conservative variables.

.. math::

   \frac{\partial h}{\partial t} + \frac{\partial hu}{\partial x} + \frac{\partial hv}{\partial y} = 0

   \frac{\partial hu}{\partial t} + \frac{\partial }{\partial x} (hu^2 + 0.5 g h^2) + \frac{\partial hv}{\partial y} = 0

   \frac{\partial hv}{\partial t} + \frac{\partial hu^2}{\partial x} + \frac{\partial }{\partial y} (hv + 0.5 g h^2) = 0


- Domain is ``[-5,5]^2`` with slip-wall BCs

- IC: :math:`h = 1 + \alpha \exp( -(x-1)^2 - (y-1)^2), u = v = 0`

* Default settings:

  - :math:`g = 9.8` (gravity parameter)

  - :math:`\mu = -3.0` (Coriolis parameter)

  - :math:`\alpha = 1/8` (initial pulse magnitude)

- :math:`g`, :math:`\mu`, :math:`\alpha` can be provided to the problem constructor (more below)


Mesh
----

.. code-block:: bash

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem swe2dSlipWall_s{3,5,7} -n Nx Ny --outDir ...

where ``Nx, Ny`` are the number of cells, and ``<stencilSize> = 3 or 5 or 7``.


C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::Swe2d::SlipWall;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5

   // A. constructor for problem using default values
   {
     auto problem = pda::create_problem_eigen(meshObj, probId, scheme);
   }

   // B. constructor for problem specifying two out of three coefficients
   {
     using scalar_type = typename decltype(meshObj)::scalar_t;
     const scalar_type gravity = /* some value */;
     const scalar_type coriolis = /* some value */;

     auto problem = pda::create_problem_eigen(meshObj, probId, scheme, gravity, coriolis);
   }

   // C. constructor for problem specifying all coefficients
   {
     using scalar_type   = typename decltype(meshObj)::scalar_t;
     const scalar_type gravity = /* some value */;
     const scalar_type coriolis = /* some value */;
     const scalar_type initPulseMagnitude = /* some value */;

     auto problem = pda::create_problem_eigen(meshObj, probId, scheme,
					      gravity, coriolis, initPulseMagnitude);
   }

   auto state = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.Swe2d::SlipWall;
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5

   # A. constructor for problem using default values
   problem = pda.create_problem(meshObj, probId, scheme)

   # B. constructor for problem specifying two out of three coefficients
   gravity  = ...
   coriolis = ...
   problem = pda.create_problem(meshObj, probId, scheme, gravity, coriolis)

   # C. constructor for problem specifying all coefficients
   gravity  = ...
   coriolis = ...
   initPulseMagnitude = ...
   problem = pda.create_problem(meshObj, probId, scheme, gravity, coriolis, initPulseMagnitude)

   state   = problem.initialCondition()



Sample Plot
-----------

Representative *height* plot as a function of time at `x=y=0`
using a 65x65 mesh with Weno5 and RK4 time integration:

.. image:: ../../figures/wiki_2dswe_height.png
  :width: 60 %
  :alt: Alternative text
  :align: center
