1D single-species reaction diffusion
====================================

This problem focuses on the following 1d diffusion reaction PDE:

.. math::

   \frac{\partial y}{\partial t} = D \frac{\partial^2 y}{\partial x^2} + k y^2 + u(x, t)


* Adapted from `this paper <https://arxiv.org/abs/1910.03193>`_

* ``D, k, u(x, t)`` can be provided to the problem constructor (more below)

* IC: :math:`y(x, 0) = 0`

* Default settings:

  - ``D = 0.01``

  - ``k = 0.01``

  - :math:`u(x, t) = \sin(\pi x) x^2 4 \cos(4 \pi x)`

* Domain is ``[0,1]`` with homogenous Dirichlet BC


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
	  --problem diffreac1d -n <N> --outDir <destination-path>

where ``N`` is the number of cells you want, and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/diffusion_reaction1d.hpp"

   namespace pda      = pressiodemoapps;
   const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("path-to-mesh");

   // A. constructor for problem using default values
   {
     const auto probId  = pda::DiffusionReaction1d::ProblemA;
     auto problem = pda::create_problem_eigen(meshObj, probId);
   }

   // B. setting custom coefficients
   {
     using scalar_type = typename decltype(meshObj)::scalar_t;
     const auto diffCoeff = static_cast<scalar_type(0.5);
     const auto reacCoeff = static_cast<scalar_type(0.2);
     auto problem = pda::create_diffusion_reaction_1d_problem_A_eigen(meshObj, diffCoeff, reacCoeff);
   }

   // C. setting custom coefficients and custom source function
   {
     using scalar_type = typename decltype(meshObj)::scalar_t;
     const auto diffCoeff = static_cast<scalar_type(0.5);
     const auto reacCoeff = static_cast<scalar_type(0.2);

     auto mySource = [](const scalar_type & x,
			const scalar_type & t,
			scalar_type & result)
     {
       // x, t are the location and time where the source must be evaluated
       result = /* do whatever you want */;
     };

     auto problem = pda::create_diffusion_reaction_1d_problem_A_eigen(meshObj, mySource,
								      diffCoeff, reacCoeff);
   }


Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   # A. constructor for problem using default values
   probId  = pda.DiffusionReaction1d.ProblemA
   problem = pda.create_problem(meshObj, probId, scheme)

   # B. setting custom coefficients
   myD, myK = 0.2, 0.001
   problem = pda.create_diffusion_reaction_1d_problem_A(meshObj, myD, myK)

   # C. setting custom coefficients and custom source function
   mysource = lambda x, time : np.sin(math.pi*x) *x*x * 4.*np.cos(4.*math.pi*x)
   problem = pda.create_diffusion_reaction_1d_problem_A(meshObj, mysource, 0.2, 0.001)


Notes:
------

.. important::

   Note that this problem does not have advection, so only viscous schemes are applicable.
   Currently, we only support a first order viscous flux
   reconstruction, which leads to a second-order scheme.
