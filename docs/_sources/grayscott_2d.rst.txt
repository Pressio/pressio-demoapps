2D Gray Scott reaction-diffusion
================================

This problem focuses on the following 2d diffusion reaction PDE:

.. math::

   \frac{\partial A}{\partial t} &=
   D_a \left(\frac{\partial^2 A}{\partial x^2} + \frac{\partial^2 A}{\partial y^2}\right)
   - AB^2 + F(1-A)

   \frac{\partial B}{\partial t} &=
   D_b \left(\frac{\partial^2 B}{\partial x^2} + \frac{\partial^2 B}{\partial y^2}\right)
   + AB^2 - (F+K)B


* :math:`D_a, D_b, F, K` can be provided to the constructor (more below)

* IC:

  - :math:`A=0.5, B=0.25` if :math:`|x| < 0.1` and :math:`|y|< 0.1`

  - :math:`A=1, B=0`, otherwise

* Default settings:

  - :math:`D_a = 0.0002, D_b = D_a/4`

  - :math:`F=0.042, K=0.062`

* Domain is ``[-1.25,1.25]^2`` with periodic BCs

* See the following references: `link1 <https://itp.uni-frankfurt.de/~gros/StudentProjects/Projects_2020/projekt_schulz_kaefer/>`_, `link2 <https://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/>`_.


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem grayscott2d -n Nx Ny --outDir <destination-path>

where ``Nx, Ny`` is the number of cells you want along x and y,
and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/diffusion_reaction.hpp"
   // ...

   namespace pda      = pressiodemoapps;
   const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("path-to-mesh");
   const auto probId  = pda::DiffusionReaction2d::GrayScott;
   const auto scheme  = pda::ViscousFluxReconstruction::FirstOrder;

   // A. constructor for problem using default values
   {
     auto problem = pda::create_problem_eigen(meshObj, probId, scheme);
   }

   // B. setting custom coefficients
   {
     using scalar_type = typename decltype(meshObj)::scalar_t;
     const scalar_type Da       = /*some_value*/;
     const scalar_type Db       = /*some_value*/;
     const scalar_type feedRate = /*some_value*/;
     const scalar_type killRate = /*some_value*/;
     auto problem      = pda::create_problem_eigen(meshObj, probId, scheme,
						   Da, Db, feedRate, killRate);
   }


Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   # ...
   probId  = pda.DiffusionReaction2d.GrayScott
   scheme  = pda.ViscousFluxReconstruction.FirstOrder

   # A. constructor for problem using default values
   problem = pda.create_problem(meshObj, probId, scheme)

   # B. setting custom coefficients
   Da, Db, F, K = ...
   problem = pda.create_problem(meshObj, probId, scheme, Da, Db, F, K)


Notes:
------

.. important::

   Note that this problem does not have advection, so inviscid schemes are not applicable
   but only viscous schemes are. Currently, we only support a first order viscous flux
   reconstruction, which leads to a second-order scheme.


Sample Plot
-----------

Representative plots at ``t=1000`` obtained using ``dt=0.5``, Runge-Kutta4 integration,
a mesh of ``160x160``, and defaults values for ``Da, Db, F, K``.

.. image:: ../../figures/wiki_grayscott_2d_0.25_1000_rk4.png
  :width: 75 %
  :alt: Alternative text
  :align: center
