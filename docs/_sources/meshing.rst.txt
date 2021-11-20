Meshing script
==============

How to use it
-------------

Missing!



Sample Mesh
-----------

In practice, the *sample mesh is a disjoint collection of cells at which the velocity or residual vector are computed*.

Several methods exist to determine which cells to include, e.g., random sampling,
the discrete empirical interpolation method ([DEIM](https://doi.org/10.1137/090766498)),
and Gauss-Newton with approximate tensors ([GNAT](https://doi.org/10.1016/j.jcp.2013.02.028)).

The sample mesh is critical to make projection-based ROMs of nonlinear systems practical.
Its at the core of hyper-reduction methods, since it allows one to create ROMs with a computational
cost which *does not* scale with the size of the full model's state vector.

The sample mesh is used in conjunction with what we refer to as the **stencil mesh**,
which contains all cells needed to compute the velocity or residual vector on the sample mesh.
Note that, in general, the sample mesh is a subset of the stencil mesh, because to compute
the velocity or residual at a given cell, one also needs the cell-centered values at that target cell.
Several examples of sample and stencil meshes are shown below.

1D Example
^^^^^^^^^^

The figures below show a full 1D mesh for a first order cell-centered finite volume scheme, and right below
a representative sample mesh (yellow-filled cells) and the stencil mesh.

.. image:: ../../figures/readme_1dmesh.png
  :width: 75 %
  :alt: Alternative text


2D Example
^^^^^^^^^^

The second example shows a two dimensional sample mesh and a stencil mesh for a first order cell-centered finite volume scheme.
The coloring scheme is the same as in the first example.

.. image:: ../../figures/readme_2dmesh.png
  :width: 75 %
  :alt: Alternative text


3D Example
^^^^^^^^^^

The third example shows a three dimensional sample mesh and
a stencil mesh for a first order cell-centered finite volume scheme.
The coloring scheme is the same as in the previous examples. Cell IDs are omitted for clarity.

.. image:: ../../figures/readme_3dmesh.png
  :width: 75 %
  :alt: Alternative text
