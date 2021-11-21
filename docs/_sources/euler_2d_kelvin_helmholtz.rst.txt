2D Euler Kelvin-Helmholtz
=========================

- Domain is :math:`\Omega = \Omega_1 \cup \Omega_2 = [-5,5]^2` with periodic BC where:

  - :math:`\Omega_1 = [-5,5] \times [-2 + \cos( 0.8 \pi x) , 2 + \cos(0.8 \pi x)]`

  - :math:`\Omega_2`: Rest of domain

- IC is as follows:

  - :math:`\Omega_1`: :math:`\rho = 2, u = 0.5, v = 0, p = 2.5`

  - :math:`\Omega_2`: :math:`\rho = 1, u = -0.5, v = 0, p = 2.5`

- Time integration performed for 2.5 flow through units, `t \in (0, 50)`.

- Mach number in :math:`\Omega_1` : :math:`M_{\infty} = 0.377964`

- Mach number in :math:`\Omega_2` : :math:`M_{\infty} = 0.267261`

- This problem is often *unstable* for a standard Galerkin ROM


Mesh
----

.. code-block:: shell

   python create_full_mesh_for.py --problem euler2dKelvinHelmholtz_s{3,5,7} -n Nx Ny --outDir ...

where ``Nx, Ny`` are the number of cells, and ``<stencilSize> = 3 or 5 or 7``.

C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::Euler2d::KelvinHelmholtz;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.Euler1d.KelvinHelmholtz
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()


Sample Plot
-----------

Representative *density* field at ``t=50`` using a 256x256 mesh with Weno5 and RK4 time integration:

.. image:: ../../figures/wiki_2d_kelvin_helmholtz_density.png
  :width: 60 %
  :alt: Alternative text
