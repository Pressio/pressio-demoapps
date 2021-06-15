# pressio-demoapps

Suite of demo/mini apps.
While the primary purpose of this repository is to contain mini-apps
to use for testing ROMs in Pressio, this repository is self-contained in the sense
that the code here does *not* depend on Pressio.

# Building Tests C++
Requires CMake > 3.18.0 and a C++ compiler with C++14 support:

```
git clone --recursive git@github.com:Pressio/pressio-demoapps.git

export CXX=<path-to-your-CXX-compiler> #must support C++14
cd pressio-demoapps && mkdir build && cd build
cmake -DPRESSIODEMOAPPS_ENABLE_TESTS=On -DPRESSIODEMOAPPS_ENABLE_TPL_EIGEN=ON ..
make -j4
ctest -j4
```
Note that to use the code as a library, you don't need to build anything.
Just point to the include directory.

# Python bindings
Requires CMake > 3.18.0 and a C++ compiler with C++14 support:

```
git clone --recursive git@github.com:Pressio/pressio-demoapps.git

export CXX=<path-to-your-CXX-compiler> #must support C++14
cd pressio-demoapps
python setup.py install
pytest -s # this runs the tests, located in side tests_py
```


# Synopsis

Coming soon.

<!-- ## 1d Problems: Sod, Lax -->

<!-- Mesh generation: -->
<!-- ```py -->
<!-- python ./meshing_scripts/create_full_mesh_for.py --name sod1d_s<3,7> -n <N> -outDir <somewhere> -->
<!-- python ./meshing_scripts/create_full_mesh_for.py --name lax1d_s<3,7> -n <N> -outDir <somewhere> -->
<!-- ``` -->

<!-- C++ problem object syntax: -->
<!-- ```c++ -->
<!-- namespace pda      = pressiodemoapps; -->
<!-- const auto meshObj = pda::loadCellCenterUniformMeshEigen(<mesh-path-string>); -->
<!-- const auto probId  = pda::euler1dproblemsEnum::{sod, lax} -->
<!-- const auto order   = pda::reconstructionEnum::{firstOrder, fifthOrderWeno}; -->
<!-- auto appObj        = pda::createEuler1dEigen(meshObj, order, probId); -->
<!-- ``` -->

<!-- Python object syntax: -->
<!-- ```py -->
<!-- meshO    = loadCellCenterUniformMesh(meshPath) -->
<!-- probId   = euler1d.{sod, lax} -->
<!-- appObj   = createEuler1dProblem(meshO, reconstructWith.fifthOrderWeno, probId) -->
<!-- ``` -->

<!-- ## 2d problems: Sedov, Riemann -->

<!-- Mesh generation: -->
<!-- ```py -->
<!-- python ./meshing_scripts/create_full_mesh_for.py --name sedov2d_s<3,7> -n <Nx> <Ny> -outDir <somewhere> -->
<!-- python ./meshing_scripts/create_full_mesh_for.py --name riemann2d_s<3,7> -n <Nx> <Ny> -outDir <somewhere> -->
<!-- ``` -->

<!-- C++ problem object syntax: -->
<!-- ```c++ -->
<!-- namespace pda      = pressiodemoapps; -->
<!-- const auto meshObj = pda::loadCellCenterUniformMeshEigen(<mesh-path-string>); -->
<!-- const auto probId  = pda::euler2dproblemsEnum::{sedov, riemann}; -->
<!-- const auto order   = pda::reconstructionEnum::{firstOrder, fifthOrderWeno}; -->
<!-- auto appObj        = pda::createEuler1dEigen(meshObj, order, probId); -->
<!-- ``` -->

<!-- Python object syntax: -->
<!-- ```py -->
<!-- meshO    = loadCellCenterUniformMesh(meshPath) -->
<!-- probId   = euler2d.{sedov, riemann} -->
<!-- appObj   = createEuler1dProblem(meshO, reconstructWith.fifthOrderWeno, probId) -->

<!-- ``` -->

# License and Citation

**If you use this code, we only ask that you reference this github page!**

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The full license is available [here](https://pressio-demosapps.github.io/various/license/).


# Questions?
Find us on Slack: https://pressioteam.slack.com or open an issue on [github](https://github.com/Pressio/pressio-tutorials).
