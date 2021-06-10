# pressio-demoapps

Suite of demo/mini apps.
While the primary purpose of this repository is to contain mini-apps
to use for testing ROMs in Pressio, this repository is self-contained in the sense
that the code here does *not* depend on Pressio.

# Building Tests
Requires CMake > 3.18.0 and a C++ compiler with C++11 support:

```
git clone --recursive git@github.com:Pressio/pressio-demoapps.git

export CXX=<path-to-your-CXX-compiler>
cd pressio-demoapps && mkdir build && cd build
cmake -DPRESSIODEMOAPPS_ENABLE_TESTS=On ..
make -j4
ctest -j4
```


# Problems Synopsis 

## 1d Problems: Sod, Lax

Mesh generation: 
```py
python create_full_mesh_for.py --name sod1d_s<3,7> -n <N> -outDir <somewhere>
python create_full_mesh_for.py --name lax1d_s<3,7> -n <N> -outDir <somewhere>
```

C++ problem object syntax:
```c++
const auto meshObj = pda::loadCellCenterUniformMeshEigen(<mesh-path-string>);
const auto probId  = pda::euler1dproblemsEnum::{sod, lax}
const auto order   = pda::reconstructionEnum::{firstOrder, fifthOrderWeno};
auto appObj        = pda::createEuler1dEigen(meshObj, order, probId);
```

Python object syntax: 
```py
meshO    = loadCellCenterUniformMesh(meshPath)
probId   = euler1d.{sod, lax}
appObj   = createEuler1dProblem(meshO, reconstructWith.fifthOrderWeno, probId)
```

## 2d problems: Sedov, Riemann

Mesh generation: 
```py
python create_full_mesh_for.py --name sedov2d_s<3,7> -n <Nx> <Ny> -outDir <somewhere>
python create_full_mesh_for.py --name riemann2d_s<3,7> -n <Nx> <Ny> -outDir <somewhere>
```

C++ problem object syntax:
```c++
const auto meshObj = pda::loadCellCenterUniformMeshEigen(<mesh-path-string>);
const auto probId  = pda::euler2dproblemsEnum::{sedov, riemann};
const auto order   = pda::reconstructionEnum::{firstOrder, fifthOrderWeno};
auto appObj        = pda::createEuler1dEigen(meshObj, order, probId);
```

Python object syntax: 
```py
meshO    = loadCellCenterUniformMesh(meshPath)
probId   = euler2d.{sedov, riemann}
appObj   = createEuler1dProblem(meshO, reconstructWith.fifthOrderWeno, probId)

```



# License and Citation
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The full license is available [here](https://pressio-demosapps.github.io/various/license/).

# Questions?
Find us on Slack: https://pressioteam.slack.com or open an issue on [github](https://github.com/Pressio/pressio-tutorials).
