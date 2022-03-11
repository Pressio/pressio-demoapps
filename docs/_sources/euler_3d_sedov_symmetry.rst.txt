3D Euler Sedov (with symmetry)
==============================

This problem solves the *3D Euler equations* in conservative form.

- Initial conditions:
    
  - a high pressure concentrated small spherical region of radius :math:`R = 3 \min(dx, dy)`

  - :math:`\left\{\begin{matrix}\rho =1, u = 0, v = 0, p = ((\gamma-1)0.851072)/(4 \pi R^3); & r\leq R \\ \rho =1, u = 0, v = 0, p = 2.5\cdot 10^{-5}; & r>R \end{matrix}\right.`

- Domain is :math:`[0.0, 1.2]^3` with reflective BC on :math:`x=y=z=0` and homogeneous Neumann on others.

- Typically, integration is performed over :math:`t \in (0, 1.)`.

.. Caution::
   Currently, this problem only works for first order and Weno3 inviscid flux reconstruction.

Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem sedov3dsym_s{3,5} -n Nx Ny Nz --outDir <destination-path>

where ``Nx, Ny, Nz`` are the number of cells you want along :math:`x`, :math:`y` and :math:`z` respectively, and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/euler3d.hpp"
   // ...
   namespace pda     = pressiodemoapps;
   const auto probId = pda::Euler3d::SedovSymmetry;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   # ...
   probId  = pda.Euler3d.SedovSymmetry
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()
