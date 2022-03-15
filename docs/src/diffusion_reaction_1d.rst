1D single-species reaction diffusion
====================================

This problem focuses on the following 1D diffusion reaction PDE:

.. math::

   \frac{\partial \phi}{\partial t} = D \frac{\partial^2 \phi}{\partial x^2} + k \phi^2 + u(x, t)


* The problem is adapted from `this paper <https://arxiv.org/abs/1910.03193>`_

* ``D, k, u(x, t)`` can be provided to the problem constructor (more below)

* Initial condition: :math:`\phi(x, 0) = 0`

* Default settings:

  - ``D = 0.01``

  - ``k = 0.01``

  - :math:`u(x, t) = 4 x^2\sin(\pi x) \cos(4 \pi x)`

* Domain is :math:`[0,1]` with homogenous Dirichlet BC


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
	  --problem diffreac1d -n <N> --outDir <destination-path>

where: 

- ``N`` is the number of cells you want

- ``<destination-path>`` is where you want the mesh files to be generated.
  The script creates the directory if it does not exist.


Notes:
------

.. important::

   For this problem, only viscous schemes are applicable.
   Note also that the stencil cannot be (yet) set because
   the implementation currently only supports a thee-point stencil
   that allows a first order viscous flux reconstruction,
   yielding a second-order scheme. We might relax this in 
   the future if additional scheme are added.


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

   meshObj = pda.load_cellcentered_uniform_mesh_eigen("path-to-mesh")

   # A. constructor for problem using default values
   probId  = pda.DiffusionReaction1d.ProblemA
   problem = pda.create_problem(meshObj, probId, scheme)

   # B. setting custom coefficients
   myD, myK = 0.2, 0.001
   problem = pda.create_diffusion_reaction_1d_problem_A(meshObj, myD, myK)

   # C. setting custom coefficients and custom source function
   myD, myK = 0.55, 0.002
   mysource = lambda x, time : np.sin(math.pi*x) *x*x * 4.*np.cos(4.*math.pi*x)
   problem = pda.create_diffusion_reaction_1d_problem_A(meshObj, probId, scheme, \
							mysource, myD, myK)
