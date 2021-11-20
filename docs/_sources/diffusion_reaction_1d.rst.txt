1D Diffusion Reaction
=====================

This problem focuses on the following 1d diffusion reaction PDE:

.. math::

   \frac{\partial s}{\partial t} = D \frac{\partial^2 s}{\partial x^2} + k s^2 + u(x, t)


* Adapted from `this paper <https://arxiv.org/abs/1910.03193>`_

* ``D, k, u(x, t)`` can be provided to the problem constructor (more below)

* IC: :math:`s(x, 0) = 0`

* Default settings:

  - ``D = 0.01``

  - ``k = 0.01``

  - :math:`u(x, t) = \sin(\pi x) x^2 4 \cos(4 \pi x)`

* Domain is ``[0,1]`` with homogenous Dirichlet BC


Mesh
----

.. code-block:: shell

   python create_full_mesh_for.py --problem diffreac1d -n <N> --outDir ...

where ``N`` is the number of cells you want.


C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::DiffusionReaction1d::ProblemA;
   const auto scheme = pressiodemoapps::ViscousFluxReconstruction::FirstOrder;

   // A. constructor for problem using default values
   {
     auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   }

   // B. setting custom coefficients, you can do
   {
     using scalar_type = typename decltype(problem)::scalar_type;
     const scalar_type diffCoeff = 0.5;
     const scalar_type reacCoeff = 0.2;
     auto problem      = pda::create_problem_eigen(meshObj, probId, scheme, diffCoeff, reacCoeff);
   }

   // C. setting custom coefficients and custom source function, you can do
   {
     using scalar_type = typename decltype(problem)::scalar_type;
     const scalar_type diffCoeff = 0.5;
     const scalar_type reacCoeff = 0.2;

     auto mySource = [](const scalar_type & x,
			const scalar_type & t,
			scalar_type & result)
     {
       // x, t are the location and time where the source must be evaluated
       result = /* do whatever you want */;
     };

     auto problem = pda::create_problem_eigen(meshObj, probId, scheme, mySource, diffCoeff, reacCoeff);
   }


Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.DiffusionReaction1d.ProblemA
   scheme  = pda.ViscousFluxReconstruction.FirstOrder

   # A. constructor for problem using default values
   problem = pda.create_problem(meshObj, probId, scheme)

   # B. setting custom coefficients, you can do
   myD, myK = 0.2, 0.001
   problem = pda.create_problem(meshObj, probId, scheme, myD, myK)

   # C. setting custom coefficients and custom source function, you can do
   mysource = lambda x, time : np.sin(math.pi*x) *x*x * 4.*np.cos(4.*math.pi*x)
   problem = pda.create_problem(meshObj, probId, scheme, 0.2, 0.001)


Notes:
------

.. important::

   Note that this problem does not have advection, so inviscid schemes are not applicable
   but only viscous schemes are. Currently, we only support a first order viscous flux
   reconstruction, which leads to a second-order scheme.
