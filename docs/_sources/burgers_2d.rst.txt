2D Burgers
==========

This problem solves the 2d Burgers equations:

.. math::

   \frac{\partial u}{\partial t} ...


* Domain is ``[0,1]^2`` with periodic BC


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem burgers2d -n Nx Ny --outDir <destination-path>

where ``Nx, Ny`` is the number of cells you want along x and y,
and ``<destination-path>`` is where you want the mesh files to be generated.

C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/advection_diffusion2d.hpp"
   // ...
   namespace pda      = pressiodemoapps;
   const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("path-to-mesh");

   FINISH

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   # ...
   scheme  = pda.InviscidFluxReconstruction.FirstOrder

   FINISH
