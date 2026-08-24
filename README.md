# TORTOISE

A high-performance C++20 library for solving the far-from-equilibrium Boltzmann transport equation, including the full collision operator.

---

## What it does

Most treatments of the Boltzmann transport equation approximate the collision integral — typically through a relaxation-time approximation — because evaluating it directly means integrating over a high-dimensional space at every point of the domain and at every time step. That approximation is reasonable close to equilibrium and breaks down far from it, which is precisely the regime that matters for ultrafast dynamics.

TORTOISE solves the equation with the collision operator treated directly, making far-from-equilibrium electron, phonon and spin dynamics tractable. It does this through:

- **Discontinuous Galerkin discretisation** over structured and unstructured meshes, in arbitrary dimension
- **Adaptive time stepping** (Runge–Kutta 4, Dormand–Prince 5(4))
- **Expression-template metaprogramming**, so that mathematical expressions on fields are written in natural notation at the call site and evaluated lazily, without temporaries or runtime abstraction cost
- **High-dimensional collision integrals** evaluated with dedicated quadrature schemes

The numerical schemes are described in the papers listed under [References](#references).

## Requirements

| | |
|---|---|
| **Compiler** | C++20 (`-std=c++20`), built at `-O3` — the library's design depends on the optimiser |
| **[Eigen](https://eigen.tuxfamily.org)** | Required. Must be on the compiler's header search path |
| **[GNUPlot](http://www.gnuplot.info)** | Required for plotting. Its full path must be set in `TORTOISE/include/Configuration` |

## Getting started

TORTOISE is a library, so you supply the `main` file. Include the whole library with:

```cpp
#include <TORTOISE>
```

To build your own project against it:

1. Add `TORTOISE/include` and your Eigen directory to the header search paths.
2. Compile every file under `TORTOISE/src` as part of your target.
3. Build with `-std=c++20 -O3`.

A Makefile is provided in `TORTOISE/doc` if you would rather not configure an IDE.

Set GNUPlot's path in `TORTOISE/include/Configuration` before running anything that plots.

## Learning the library

- **`TORTOISE/tutorials/tutorial1` … `tutorial5`** — -guided introduction, in order. Start here.
- **`TORTOISE/tutorials/developers`** — internals, for extending the library rather than using it.
- **`TORTOISE/examples/main.cpp`** — a complete worked problem.

## Repository layout

```
TORTOISE/
├── include/            Headers — include <TORTOISE> pulls in everything
│   ├── Generics/       Containers, algorithms (RK4, Dormand–Prince, Nelder–Mead), utilities
│   ├── Geometry/       Meshes, regions, function spaces, plotting
│   │   ├── Structured/
│   │   └── Unstructured/
│   └── Physics/        Materials, bands, scattering channels, DG integrators (Vlasov)
├── src/                Implementation — all of this must be compiled into your target
├── plugins/            Band structures and optional algorithms
├── examples/           A complete worked problem
├── tutorials/          Numbered tutorials, plus developer-facing internals
└── doc/                Makefile and documentation
```

## References

The numerical methods implemented here are described in:

1. M. Wais, K. Held, M. Battiato, *Numerical solver for the time-dependent far-from-equilibrium Boltzmann equation*, **Computer Physics Communications 264, 107877 (2021)**
2. I. Wadgaonkar, M. Wais, M. Battiato, *Numerical solver for the out-of-equilibrium time-dependent Boltzmann collision operator: application to 2D materials*, **Computer Physics Communications 271, 108207 (2022)**
3. I. Wadgaonkar, R. Jain, M. Battiato, *Numerical scheme for the far-out-of-equilibrium time-dependent Boltzmann collision operator: 1D second-degree momentum discretisation and adaptive time stepping*, **Computer Physics Communications 263, 107863 (2021)**

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).

## Author

Marco Battiato — <battiato.marco@gmail.com>
