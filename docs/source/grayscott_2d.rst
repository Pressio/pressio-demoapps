2D Gray Scott reaction-diffusion
================================

This problem focuses on the following 2D reaction-diffusion system of PDE:

.. math::

   \frac{\partial A}{\partial t} &=
   D_a \left(\frac{\partial^2 A}{\partial x^2} + \frac{\partial^2 A}{\partial y^2}\right)
   - AB^2 + F(1-A)

   \frac{\partial B}{\partial t} &=
   D_b \left(\frac{\partial^2 B}{\partial x^2} + \frac{\partial^2 B}{\partial y^2}\right)
   + AB^2 - (F+K)B


* :math:`D_a, D_b, F, K` can be provided to the constructor (more below)

* Initial conditions:

  :math:`\left\{\begin{matrix}A=1/2, B=1/4; & \text{if }|x| < 1/10 \text{ and } |y|< 1/10 \\ A=1, B=0; & \text{otherwise} \end{matrix}\right.`

* Default settings:

  - :math:`D_a = 2\cdot 10^{-4}, D_b = D_a/4`

  - :math:`F=0.042, K=0.062`

* Domain is :math:`[-5/4,5/4]^2` with periodic BCs

* For more details on the problem, see the following references: `link1 <https://itp.uni-frankfurt.de/~gros/StudentProjects/Projects_2020/projekt_schulz_kaefer/>`_, `link2 <https://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/>`_.


Mesh
----

.. code-block:: shell

   python3 pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
           --problem grayscott2d -n Nx Ny --outDir <destination-path>

where ``Nx, Ny`` is the number of cells you want along :math:`x` and :math:`y` respectively,
and ``<destination-path>`` is where you want the mesh files to be generated.


C++ synopsis
------------

.. code-block:: c++

   #include "pressiodemoapps/diffusion_reaction2d.hpp"

   int main(){
     namespace pda      = pressiodemoapps;
     const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen("path-to-mesh");
     const auto scheme  = pda::ViscousFluxReconstruction::FirstOrder;

     // A. constructor for problem using default values
     {
       const auto probId  = pda::DiffusionReaction2d::GrayScott;
       auto problem = pda::create_problem_eigen(meshObj, probId, scheme);
     }

     // B. setting custom coefficients
     {
       using scalar_type = typename decltype(meshObj)::scalar_t;
       const scalar_type Da  = /*some_value*/;
       const scalar_type Db  = /*some_value*/;
       const scalar_type F   = /*some_value*/;
       const scalar_type K   = /*some_value*/;
       auto problem      = pda::create_gray_scott_2d_problem(meshObj, scheme, Da, Db, F, K);
     }
   }


Python synopsis
---------------

.. code-block:: py

   import pressiodemoapps as pda

   scheme  = pda.ViscousFluxReconstruction.FirstOrder

   # A. constructor for problem using default values
   probId  = pda.DiffusionReaction2d.GrayScott
   problem = pda.create_problem(meshObj, probId, scheme)

   # B. setting custom coefficients
   Da, Db, F, K = ...
   problem = pda.create_gray_scott_2d_problem(meshObj, scheme, Da, Db, F, K)


Notes:
------

.. important::

   Note that this problem does not have advection, so inviscid schemes are not applicable
   but only viscous schemes are. Currently, we only support a first order viscous flux
   reconstruction, which leads to a second-order scheme.


Sample Plot
-----------

Representative plots at selected time :math:`t=1000` obtained using time step :math:`dt=0.5`, Runge-Kutta4 integration,
a mesh of ``160x160`` and default values for :math:`D_a, D_b, F, K`.

.. image:: ../../figures/wiki_grayscott_2d_0.25_1000_rk4.png
  :width: 75 %
  :alt: Alternative text
  :align: center
