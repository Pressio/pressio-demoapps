2D Burgers
==========
TODO: Add Initial conditions and parameters. Fix code blocks.

This problem solves the 2D Burgers equations. We consider a two-dimensional nonlinear viscous Burgers equations

.. math::

   \begin{matrix} \frac{\partial u}{\partial t} + u\frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} = \frac{1}{Re}\Big(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}\Big)\\ \frac{\partial v}{\partial t} + u\frac{\partial v}{\partial x} + v \frac{\partial v}{\partial y} = \frac{1}{Re}\Big (\frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2}\Big) \end{matrix}

subject to initial and boundary conditions.

* Domain is :math:`[0,1]^2` with periodic BC

Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem burgers2d -n Nx Ny --outDir <destination-path>

where ``Nx, Ny`` is the number of cells you want along :math:`x` and :math:`y` respectively,
and ``<destination-path>`` is where you want the mesh files to be generated.

C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/advection_diffusion2d.hpp"
   // ...
   namespace pda      = pressiodemoapps;
   const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("path-to-mesh");
   const auto scheme  = pda::InviscidFluxScheme::Rusanov;

   // A. constructor for problem using default values
   {
     const auto probId  = pda::AdvectionDiffusion2d::Burgers;
     auto problem = pda::create_problem_eigen(meshObj, probId, scheme);
   }

   // B. setting custom coefficients
   {
     using scalar_type = typename decltype(meshObj)::scalar_t;
     const auto icPulseMagnitude = static_cast<scalar_t>(0.5);
     const auto icSpread = static_cast<scalar_t>(0.15);
     const auto diffusion = static_cast<scalar_t>(0.00001);
     const auto icCenterX = static_cast<scalar_t>(0.0);
     const auto icCenterY = static_cast<scalar_t>(-0.2);
     auto problem = pda::create_burgers_2d_problem_eigen(meshObj, scheme,
					icPulseMagnitude, icSpread, diffusion,
		                        icCenterX, icCenterY));
   }

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   # ...
   scheme  = pda.InviscidFluxReconstruction.FirstOrder

   # A. constructor for problem using default values
   probId  = pda.AdvectionDiffusion2d::Burgers
   problem = pda.create_problem(meshObj, probId, scheme)

   # B. setting custom coefficients
   icPulseMagnitude = 0.
   icSpread = 0.15
   diffusion = 0.00001
   icCenterX = 0.0
   icCenterY = -0.2
   problem = pda.create_burgers_2d_problem(meshObj, scheme, icPulseMagnitude, icSpread, diffusion, icCenterX, icCenterY)
