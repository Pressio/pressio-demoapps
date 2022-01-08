2D Euler Double Mach Reflection
===============================

- `Reference paper <http://flash.uchicago.edu/site/flashcode/user_support/flash_ug_devel/node184.html#SECTION010116000000000000000>`_

- IC is a Mach 10 shock in air tilted by an angle, see paper above.

- Domain is [0.0, 4.0]x[0.0, 1.0]; for BC see link above.

- Typically, integration is performed for `t \in (0, 0.25)`.


.. Caution::
   Currently, this problem only works for first order and Weno3 inviscid flux reconstruction.


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem doublemach2d_s{3,5} -n Nx Ny --outDir <destination-path>

where ``Nx, Ny`` are the number of cells, and ``<stencilSize> = 3 or 5 or 7``,
and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   namespace pda = pressiodemoapps;
   const auto probId = pda::Euler2d::DoubleMachReflection;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
   auto state	     = problem.initialCondition();

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   probId  = pda.Euler2d.DoubleMachReflection
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()


Sample Plot
-----------

Representative *density* plot at `t=0.25` using a 600x150 mesh with Weno3 and SSPRK3 time integration:

.. image:: ../../figures/wiki_2d_dmr_density.png
  :width: 80 %
  :alt: Alternative text
  :align: center
