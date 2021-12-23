2D Euler Sedov Full
===================

* `Reference paper <https://www.researchgate.net/publication/260967068_GENASIS_General_Astrophysical_Simulation_System_I_Refinable_Mesh_and_Nonrelativistic_Hydrodynamics>`_

- IC is a high pressure concentrated small spherical region of radius ``R = 2. * min(dx, dy)``

  - For :math:`r<=R`:

    - :math:`\rho =1, u = 0, v = 0, p = (\gamma-1)/(pi*R*R)`

  - For :math:`r>R`:

    - :math:`\rho =1, u = 0, v = 0, p = 5e-5`

  :math:`dx, dy` are the cell widths and the factor of 2 is used to spread the source over 2 cells for numerical reasons, and :math:`r = sqrt(x^2+y^2)` is the distance from the origin.

- Domain is ``[-1.2, 1.2]^2`` with homogeneous Neumann BC on all boundaries

- Typically, integration is performed for :math:`t \in (0, 1.)`.


Mesh
----

.. code-block:: shell

   cd pressio-demoapps/meshing_scripts
   python create_full_mesh_for.py --problem sedov2d_s{3,5,7} -n Nx Ny --outDir <destination-path>

where ``Nx, Ny`` are the number of cells, and ``<stencilSize> = 3 or 5 or 7``, 
and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::Euler2d::SedovFull;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.Euler2d.SedovFull
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()


Sample Plot
-----------

Representative *pressure* plot at `t=0.2` using a 200x200 mesh with Weno3 and SSPRK3 time integration:

.. image:: ../../figures/wiki_2d_sedov_pressure.png
  :width: 60 %
  :alt: Alternative text
