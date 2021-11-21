2D Diffusion Reaction
=====================

This problem focuses on the following 2d diffusion reaction PDE:

.. math::

   \frac{\partial s}{\partial t} = D \left(\frac{\partial^2 s}{\partial x^2}
   + \frac{\partial^2 s}{\partial y^2} \right) + k s^2 + u(x, y, t)


* ``D, k, u(x, y, t)`` can be provided to the problem constructor (more below)

* IC: :math:`s(x, y, 0) = 0`

* Default settings:

  - ``D = 0.01``

  - ``k = 0.01``

  - :math:`u(x, y, t) = 4 \sin(4 \pi x y) \sin(\pi x (y-0.2))`

* Domain is ``[0,1]^2`` with homogenous Dirichlet BC


Mesh
----

.. code-block:: shell

   python create_full_mesh_for.py --problem diffreac2d -n Nx Ny --outDir ...

where ``Nx, Ny`` is the number of cells you want along x and y.


C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::DiffusionReaction2d::ProblemA;
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
			const scalar_type & y,
			const scalar_type & t,
			scalar_type & result)
     {
       // x, y, t are the location and time where the source must be evaluated
       result = /* do whatever you want */;
     };

     auto problem = pda::create_problem_eigen(meshObj, probId, scheme, mySource, diffCoeff, reacCoeff);
   }


Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.DiffusionReaction2d.ProblemA
   scheme  = pda.ViscousFluxReconstruction.FirstOrder

   # A. constructor for problem using default values
   problem = pda.create_problem(meshObj, probId, scheme)

   # B. setting custom coefficients, you can do
   myD, myK = 0.2, 0.001
   problem = pda.create_problem(meshObj, probId, scheme, myD, myK)

   # C. setting custom coefficients and custom source function, you can do
   mysource = lambda x, y, time : np.sin(math.pi*x) + y * x + time # or whatever
   problem = pda.create_problem(meshObj, probId, scheme, mysource, 0.2, 0.001)


Notes:
------

.. important::

   Note that this problem does not have advection, so inviscid schemes are not applicable
   but only viscous schemes are. Currently, we only support a first order viscous flux
   reconstruction, which leads to a second-order scheme.
