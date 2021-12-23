Full mesh
=========

The ``pressio-demoapps/meshing`` subdirectory contains the following scripts:

- ``create_full_mesh.py``: for generating a full mesh

- ``create_full_mesh_for.py``: for generating a full mesh for a **specific problem**

For all the above, you can see the corresponding info by typing ``script-name.py --help``.
Below, we describe the scripts in more details.


Script Usage
------------

.. code-block:: shell

   python create_full_mesh.py \
		--outDir <path-to-where-you-want-mesh-to-be-generated> \
		-n Nx [Ny Nz] \
		--bounds xMin xMax [yMin yMax zMin zMaz] \
		-s stencilSize        # valid choices are 3,5,7 \
		--periodic true/false # default=false \
		--debug true/false    # default=false

- ``--outDir``: full path to directory where you want all mesh files to be generated;

- ``-n``: number of cells along each axis. If you only pass one value, script assumes 1d domain.
  If you pass two values, script assumes 2d. If you pass three values, assumes 3d;

- ``-s``: target stencil size;


**finish**


Problem-specific mesh script
----------------------------

.. code-block:: shell

   python create_full_mesh_for.py --problem lax1d_s<stencilSize> -n <N> --outDir ...

where ``N`` is the number of cells you want and ``<stencilSize> = 3 or 5 or 7``.


**finish**
