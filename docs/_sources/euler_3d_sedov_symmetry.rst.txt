3D Euler Sedov (with symmetry)
==============================


- IC involves a small spherical region of radius :math:`R = 3. * min(dx, dy)` with high pressure

  - :math:`r<=R`: :math:`\rho = 1, u = 0, v = 0, p = (4 (\gamma-1) 0.851072)/(4 pi R^3)`

  - :math:`r>R`: :math:`\rho = 1, u = 0, v = 0, p = 2.5e-5`

- Domain is ``[0.0, 1.2]^3`` with reflective BC on `x=y=z=0`, and homogeneous Neumann on the others.

- Typically, integration is performed over :math:`t \in (0, 1.)`.

.. Caution::
   Currently, this problem only works for first order and Weno3 inviscid flux reconstruction.

Mesh
----

.. code-block:: shell

   cd pressio-demoapps/meshing_scripts
   python create_full_mesh_for.py --problem sedov3dsym_s{3,5} -n Nx Ny Nz --outDir <destination-path>

where ``Nx, Ny`` are the number of cells, and ``<stencilSize> = 3 or 5``, 
and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::Euler3d::SedovSymmetry;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.Euler3d.SedovSymmetry
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()
