2D Euler Riemann
================

This problem solves the *2D conservative Euler equations*.

- Domain is :math:`[0, 1]^2` with homogeneous Neumann on all boundaries

- Various initial conditions are supported:

  - ``icId==1``:

    - See page 15 of `paper1 <https://www.researchgate.net/publication/269636534_A_Compact_Third-Order_Gas-Kinetic_Scheme_for_Compressible_Euler_and_Navier-Stokes_Equations>`_

      :math:`\left\{\begin{matrix}\rho = 0.5313, u = 0, v = 0, p = 0.4; & x\geq 1/2, y\geq 1/2\\ \rho = 1, u = 0.7276, v = 0, p = 1; & x<1/2, y\geq 1/2 \\ \rho = 4/5, u = 0, v = 0, p = 1; & x<1/2, y<1/2 \\ \rho = 1, u = 0, v = 0.7276, p = 1;& x>1/2, y<1/2 \end{matrix}\right.`

    - This IC is used to create the corresponding initial conditions in conservative variables.

    - Time integration is performed for :math:`t \in (0, 1/5)`


  - ``icId==2``:

    - See configuration 3 of `paper2 <http://www.amsc-ouc.ac.cn/Files/Papers/2016_Don_Hybrid%20Compact-WENO%20finite%20difference%20scheme%20with%20conjugate%20Fourier%20shock%20detection%20algorithm%20for%20hyperbolic%20conservation%20laws.pdf>`_
  
      :math:`\left\{\begin{matrix}\rho = 1.5, u = 0, v = 0, p = 1.5; & x\geq 4/5, y\geq 4/5\\ \rho = 0.5323, u = 1.206, v = 0, p = 0.3; & x<4/5, y\geq 4/5 \\ \rho = 0.138, u = 1.206, v = 1.206, p = 0.029; &x<4/5, y<4/5 \\ \rho = 0.5323, u = 0, v = 1.206, p = 0.3;& x>4/5, y<4/5 \end{matrix}\right.`

    - This IC is used to create the corresponding initial conditions in conservative variables.

    - Time integration is performed for :math:`t \in (0, 4/5)`


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem riemann2d_s{3,5,7} -n Nx Ny --outDir <destination-path>

where ``Nx, Ny`` are the number of cells you want along :math:`x` and :math:`y` respectively, and ``<stencilSize> = 3 or 5 or 7``,
and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/euler2d.hpp"
   // ...
   namespace pda     = pressiodemoapps;
   const auto probId = pda::Euler2d::Riemann;
   const auto scheme = pda::InviscidFluxReconstruction::FirstOder; //or Weno3, Weno5
   auto problem      = pda::create_problem_eigen(meshObj, probId, scheme [, icId]);
   auto state	     = problem.initialCondition();

Where the ``icId`` is an integer identifying the initial condition above.

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda
   # ...
   probId  = pda.Euler2d.Riemann
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3, Weno5
   problem = pda.create_problem(meshObj, probId, scheme [, icId])
   state   = problem.initialCondition()

Where the ``icId`` is an integer identifying the initial condition above.


Sample Plot
-----------

Representative *density* plot at :math:`t=4/5` using ``icId=2`` initial conditions with Weno5,
SSPRK3 time integration:

.. image:: ../../figures/wiki_2d_riemann_density.png
  :width: 60 %
  :alt: Alternative text
  :align: center
