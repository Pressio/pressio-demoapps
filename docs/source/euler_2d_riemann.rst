2D Euler Riemann
================

This problem solves the *2D conservative Euler equations*

.. math::

   \frac{\partial }{\partial t} \begin{bmatrix}\rho \\ \rho u_x \\ \rho u_y\\ \rho E \end{bmatrix} + \frac{\partial }{\partial x} \begin{bmatrix}\rho u_x \\ \rho u_x^2 +p \\ \rho u_x u_y \\ (E+p)u_x \end{bmatrix} \frac{\partial }{\partial y} \begin{bmatrix}\rho u_y  \\ \rho u_x u_y \\ \rho u_y^2 +p \\ (E+p)u_y \end{bmatrix}= 0

where the pressure :math:`p` is related to the conserved quantities through the equation of the state

.. math::

   p=(\gamma -1)(\rho E-\frac{1}{2}\rho (u_x^2 + u_y^2)).


- Various initial conditions are supported, with various default parameterizations:

  - ``icId==1``:

    - Accepts the following user-specified parameters

      - ``riemannTopRightPressure``: pressure in the upper right quadrant (default: 0.4)

    - Default setting follows those on page 15 of `paper1 <https://www.researchgate.net/publication/269636534_A_Compact_Third-Order_Gas-Kinetic_Scheme_for_Compressible_Euler_and_Navier-Stokes_Equations>`_

      :math:`\left\{\begin{matrix}\rho = 0.5313, u = 0, v = 0, p = 0.4; & x\geq 1/2, y\geq 1/2\\ \rho = 1, u = 0.7276, v = 0, p = 1; & x<1/2, y\geq 1/2 \\ \rho = 4/5, u = 0, v = 0, p = 1; & x<1/2, y<1/2 \\ \rho = 1, u = 0, v = 0.7276, p = 1;& x>1/2, y<1/2 \end{matrix}\right.`

    - This IC is used to create the corresponding initial conditions in conservative variables.

    - Time integration is performed for :math:`t \in (0, 1/5)`


  - ``icId==2``:

    - Accepts the following user-specified parameters

      - ``riemannTopRightPressure``: pressure in upper right quadrant (default: 1.5)

      - ``riemannTopRightXVel``: x-velocity in upper right quadrant (default: 0.0)

      - ``riemannTopRightYVel``: y-velocity in upper right quadrant (default: 0.0)

      - ``riemannTopRightDensity``: density in upper right quadrant (default: 1.5)

      - ``riemannBotLeftPressure``: pressure in lower right quadrant (default: 0.029)

    - The remaining, unspecified primitive state values are computed following the compatibility relations provided by Configuration 3 of `paper2 <https://epubs.siam.org/doi/pdf/10.1137/0524006>`_

    - Default settings follow those in configuration 3 of `paper3 <http://www.amsc-ouc.ac.cn/Files/Papers/2016_Don_Hybrid%20Compact-WENO%20finite%20difference%20scheme%20with%20conjugate%20Fourier%20shock%20detection%20algorithm%20for%20hyperbolic%20conservation%20laws.pdf>`_

      :math:`\left\{\begin{matrix}\rho = 1.5, u = 0, v = 0, p = 1.5; & x\geq 4/5, y\geq 4/5\\ \rho = 0.5323, u = 1.206, v = 0, p = 0.3; & x<4/5, y\geq 4/5 \\ \rho = 0.138, u = 1.206, v = 1.206, p = 0.029; &x<4/5, y<4/5 \\ \rho = 0.5323, u = 0, v = 1.206, p = 0.3;& x>4/5, y<4/5 \end{matrix}\right.`

    - This IC is used to create the corresponding initial conditions in conservative variables.

    - Time integration is performed for :math:`t \in (0, 4/5)`


- By default, :math:`\gamma = 1.4`

- Domain is :math:`[0, 1]^2` with homogeneous Neumann on all boundaries


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem riemann2d_s<stencilSize> -n Nx Ny --outDir <destination-path>

where

- ``Nx, Ny`` is the number of cells you want along :math:`x` and :math:`y` respectively

- ``<stencilSize> = 3 or 5 or 7``: defines the neighboring connectivity of each cell

- ``<destination-path>`` is where you want the mesh files to be generated.
  The script creates the directory if it does not exist.


.. Important::

  When you set the ``<stencilSize>``, keep in mind the following constraints (more on this below):

  - ``InviscidFluxReconstruction::FirstOrder`` requires ``<stencilSize> >= 3``

  - ``InviscidFluxReconstruction::Weno3`` requires ``<stencilSize> >= 5``

  - ``InviscidFluxReconstruction::Weno5`` requires ``<stencilSize> >= 7``


C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/euler2d.hpp"

   int main(){
     namespace pda     = pressiodemoapps;

     const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("path-to-mesh");

     const auto probId = pda::Euler2d::Riemann;
     const auto scheme = pda::InviscidFluxReconstruction::FirstOrder; //or Weno3, Weno5
     auto problem      = pda::create_problem_eigen(meshObj, probId, scheme [, icId, userParams]);
     auto state	     = problem.initialCondition();
   }

Where the ``icId`` is an integer identifying the initial condition above, and ``userParams`` is an ``std::unordered_map<std::string, decltype(meshObj)::scalar_type>`` mapping user-specified parameter names to their designated values (see above for valid parameter names and their default values for each ``icId``).

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   meshObj = pda.load_cellcentered_uniform_mesh("path-to-mesh")

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
