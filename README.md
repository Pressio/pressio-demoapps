
# Overview

**pressio-demoapps** is a library of 1D, 2D and 3D demo problems
of varying complexity, ranging from a simple 1D linear advection,
to 2D reaction-diffusion, and 3D Euler equations, and more.

Key features include:

- support for both C++ and Python
- cell-centered finite volume discretization with various numerical schemes and *exact Jacobians*
- focus on providing self-contained and well-defined problems
- built-in support for a sample mesh: this mean that one can evaluate the residual and Jacobian
at a disjoint subset of the mesh cells (this is useful for intrusive ROMs, but in general for other purposes too)

Click below to check the documentation for more details:

<a href="https://pressio.github.io/pressio-demoapps/index.html" target="_blank">
    <img src='figures/logo-display.svg' width='90%'>
</a>

## Development status

**pressio-demoapps** is planned to be maintained in the long-term.
Hopefully, more problems and features will be implemented.
If you are interested in collaborating or would like to see a specific problem added, please reach out.

## Questions?
Find us on Slack: https://pressioteam.slack.com or open an issue on [github](https://github.com/Pressio/pressio-demoapps).

## License and Citation

While we work on publishing this, **if you use this code, please reference this Github page!**

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

The full license is available [here](https://github.com/Pressio/pressio-demoapps/blob/main/LICENSE).
