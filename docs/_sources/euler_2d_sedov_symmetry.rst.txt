2D Euler Sedov (with symmetry)
==============================

* `Reference paper1 <https://reader.elsevier.com/reader/sd/pii/S002199911400477X?token=658F08D28B5C7A6A97E6F4478FD494699F3C8DF23970A256F06E501B7B136F9A6A540EEA749F28AC2AF4A6A7993A8517&originRegion=eu-west-1&originCreation=20210611123033>`_

* `Reference paper2 <http://flash.uchicago.edu/site/flashcode/user_support/flash_ug_devel/node184.html#SECTION010114000000000000000>`_

- IC involves a small spherical region of radius `R = 3. * min(dx, dy)` with high pressure

  - :math:`r<=R`: :math:`\rho =1, u = 0, v = 0, p = ((\gamma-1)*0.851072)/(\pi R^2)`

  - :math:`r>R`: :math:`\rho =1, u = 0, v = 0, p = 2.5e-5`

- Domain is [0.0, 1.2]^2 with reflective BC on `y=0` and `x=0` and homogeneous Neumann for `x=1.2` and `y=1.2`

- Typically, integration is performed for :math:`t \in (0, 1.)`.


Mesh
----

.. code-block:: shell

   python create_full_mesh_for.py --problem sedov2dsym_s{3,5,7} -n Nx Ny --outDir ...

where ``Nx, Ny`` are the number of cells, and ``<stencilSize> = 3 or 5 or 7``.


C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::Euler2d::SedovSymmetry;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.Euler1d.SedovSymmetry
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()
