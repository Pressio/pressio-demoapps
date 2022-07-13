2D Euler Double Mach Reflection
===============================

This problem solves the *2D conservative Euler equations*

.. math::

   \frac{\partial }{\partial t} \begin{bmatrix}\rho \\ \rho u_x \\ \rho u_y\\ \rho E \end{bmatrix} + \frac{\partial }{\partial x} \begin{bmatrix}\rho u_x \\ \rho u_x^2 +p \\ \rho u_x u_y \\ (E+p)u_x \end{bmatrix} \frac{\partial }{\partial y} \begin{bmatrix}\rho u_y  \\ \rho u_x u_y \\ \rho u_y^2 +p \\ (E+p)u_y \end{bmatrix}= 0

where the pressure :math:`p` is related to the conserved quantities through the equation of the state

.. math::

   p=(\gamma -1)(\rho E-\frac{1}{2}\rho (u_x^2 + u_y^2)).


- By default, :math:`\gamma = 1.4`

- Initial condition is a Mach 10 shock tilted by an angle, see reference paper above.

- Domain is :math:`[0, 4]\times[0, 1]`. For BC see link above.

- Typically, integration is performed for :math:`t \in (0, 0.25)`.

- The problem is adopted from `this paper <http://flash.uchicago.edu/site/flashcode/user_support/flash_ug_devel/node184.html#SECTION010116000000000000000>`_


.. Warning::
   Currently, this problem only works for first order and Weno3 inviscid flux reconstruction.


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem doublemach2d_s<stencilSize> -n Nx Ny --outDir <destination-path>

where

- ``Nx, Ny`` is the number of cells you want along :math:`x` and :math:`y` respectively

- ``<stencilSize> = 3 or 5``: defines the neighboring connectivity of each cell

- ``<destination-path>`` is where you want the mesh files to be generated.
  The script creates the directory if it does not exist.


.. Important::

  When you set the ``<stencilSize>``, keep in mind the following constraints (more on this below):

  - ``InviscidFluxReconstruction::FirstOrder`` requires ``<stencilSize> >= 3``

  - ``InviscidFluxReconstruction::Weno3`` requires ``<stencilSize> >= 5``


C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/euler2d.hpp"

   int main(){
     namespace pda     = pressiodemoapps;

     const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("path-to-mesh");

     const auto probId = pda::Euler2d::DoubleMachReflection;
     const auto scheme = pda::InviscidFluxReconstruction::FirstOrder; //or Weno3
     auto problem      = pda::create_problem_eigen(meshObj, probId, scheme);
     auto state	     = problem.initialCondition();
   }

Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   meshObj = pda.load_cellcentered_uniform_mesh("path-to-mesh")

   probId  = pda.Euler2d.DoubleMachReflection
   scheme  = pda.InviscidFluxReconstruction.FirstOrder # or Weno3
   problem = pda.create_problem(meshObj, probId, scheme)
   state   = problem.initialCondition()


Sample Plot
-----------

Representative *density* plot at :math:`t=0.25` using a ``600x150`` mesh with Weno3
and SSPRK3 time integration:

.. image:: ../../figures/wiki_2d_dmr_density.png
  :width: 80 %
  :alt: Alternative text
  :align: center
