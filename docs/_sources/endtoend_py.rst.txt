Step-by-step in Python
======================

This page shows how to setup, run and visualize your first problem using Python.

For demonstration purposes, we show how to run the **1D Sod** problem.
The *same* process applies to *every* other problem.

.. Hint::
   You can copy each snippet below by simply moving your mouse over
   the snippet's box, and clicking the copy icon
   appearing on the top-right corner.


Step 1: Prepare
---------------

You need CMake > 3.18.0.
Let's make a directory to run things and export the C++ compiler:

.. code-block:: shell

  export MYTEST=/home/myDemoTest
  mkdir $MYTEST
  export CXX=<path-to-your-CXX-compiler> #must support C++14

  cd $MYTEST
  git clone --recursive git@github.com:Pressio/pressio-demoapps.git
  cd pressio-demoapps
  python setup.py install


After this step, you should have demoapps installed.


Step 2: Generate the mesh for Sod1D
-------------------------------------

.. code-block:: shell

   cd $MYTEST
   python ./pressio-demoapps/meshing_scripts/create_full_mesh_for.py \
		--problem sod1d_s7 \
		--outdir ${MYTEST}/mesh \
		-n 500

where the string ``sod1d_s7`` indicates the target problem
and that we want a 7-point stencil, ``500`` is the number of cells we want.
The mesh files are generated inside ``${MYTEST}/mesh``.


Step 3: Main file and run
-------------------------

Make a ``main.py`` and copy the following code in it:

.. code-block:: py

   import pathlib, sys, numpy as np
   import matplotlib.pyplot as plt
   import pressiodemoapps as pda
   file_path = pathlib.Path(__file__).parent.absolute()

   if __name__ == '__main__':
     meshPath = str(file_path) + "/mesh"
     meshObj  = pda.load_cellcentered_uniform_mesh(meshPath)
     problem  = pda.Euler1d.Sod
     scheme   = pda.InviscidFluxReconstruction.Weno5
     appObj   = pda.create_problem(meshObj, problem, scheme)

     yn = appObj.initialCondition()
     dt = 0.001
     Nsteps = 200
     pda.advanceRK4(appObj, yn, dt, Nsteps)

     x = meshObj.viewX()
     # plot only density
     plt.plot(x, yn[0:-1:3])
     plt.show()
