# pressio-demoapps

Suite of demo/mini apps.
While the primary purpose of this repository is to contain mini-apps
to use for testing ROMs in Pressio, this repository is self-contained in the sense
that the code here does *not* depend on Pressio.

# Building Tests
Requires CMake > 3.18.0 and a C++ compiler with C++11 support:

```
git clone git@github.com:Pressio/pressio-demoapps.git

export CXX=<path-to-your-CXX-compiler>
cd pressio-demoapps && mkdir build && cd build
cmake -DPRESSIODEMOAPPS_ENABLE_TESTS=On ..
make -j4
ctest -j4
```

# License and Citation
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The full license is available [here](https://pressio-demosapps.github.io/various/license/).

# Questions?
Find us on Slack: https://pressioteam.slack.com or open an issue on [github](https://github.com/Pressio/pressio-tutorials).
