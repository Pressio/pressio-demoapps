1D Euler Sod
============

This problem solves the *1D Sod problem* for the Euler equations.

* `Reference paper <https://iopscience.iop.org/article/10.1086/317361>`_

- IC:

  - :math:`x<=0`: :math:`\rho =1, u = 0, p = 1`

  - :math:`x>0`: :math:`\rho =0.125, u = 0, p = 0.1`

- Domain is ``[-0.5, 0.5]`` with homogeneous Neumann BC

- Typically, integration is performed over :math:`t \in (0, 0.2)`.


Mesh
----

.. code-block:: shell

   cd pressio-demoapps/meshing_scripts
   python create_full_mesh_for.py --problem sod1d_s<stencilSize> -n <N> --outDir <destination-path>

where ``N`` is the number of cells you want and ``<stencilSize> = 3 or 5 or 7``, 
and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::Euler1d::Sod;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state        = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.Euler1d.Sod
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()
