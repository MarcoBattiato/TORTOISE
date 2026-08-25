# TORTOISE

A set of C++ libraries for computing **strongly out-of-equilibrium dynamics in realistic materials**.

TORTOISE solves the time evolution of the full Boltzmann equation — transport and scattering — under **no close-to-equilibrium approximation**, with every population left free to change in time. It conserves particle number, energy and momentum **exactly, to machine precision**.

---

## What it does

The Boltzmann scattering equation (also known as the quantum Fokker–Planck equation) is usually solved by approximating the collision term, most often through a relaxation-time approximation, and by treating part of the system — typically phonons — as a bath that absorbs energy. Both simplifications are reasonable near equilibrium and fail far from it.

TORTOISE removes both:

- **No close-to-equilibrium approximation.** The full scattering term is solved as it stands.
- **No bath required.** All quasiparticle populations are propagated explicitly, so proper thermalisation and partial thermalisation can be described rather than assumed. Baths can still be included if you want them.
- **Any dimensionality**, on structured or unstructured meshes.

### Why exact conservation is the point

This is the feature that makes the rest possible, and it is easy to underrate as a technical detail.

When part of a system is treated as a bath, energy leaves the simulation and small errors are washed out. When it is not — when you actually want to propagate every quasiparticle in the system — there is nothing to absorb the error. A tiny violation of energy conservation at each time step accumulates over many steps, and the numerical solution converges to the wrong equilibrium, or fails to converge at all. Predictions spanning several timescales become meaningless.

Every TORTOISE library conserves particle number, momentum and energy to machine precision. That is what allows thermalisation processes occurring over widely separated timescales to be computed in a single propagation.

### What it does not do

TORTOISE does **not** calculate dispersions or scattering amplitudes — it requires them as input. They can be supplied as analytic expressions, read from file as *ab initio* data, or computed on the fly.

---

## What using it looks like

Setting up graphene: a hexagonal Brillouin zone, two meshes around the Dirac points, and four bands defined by analytic dispersions.

```cpp
const int           nDim = 2;                  // dimensionality of the material
const realValueType protonCharge = 1.;

Point<nDim>  originBZ({0.0, 0.0}),
             sideBZ0({1.0, 0.0}),
             sideBZ1({-0.5, std::sqrt(3.0)/2});

Region<nDim>   brillouinZone(originBZ, {sideBZ0, sideBZ1});
Material<nDim> graphene("graph", brillouinZone);

const int    resolution = 10;
Mesh<nDim>   meshK0(brillouinZone, {0.1,0.1}, {0.3,0.3}, {resolution, resolution});
Mesh<nDim>   meshK1(brillouinZone, {0.6,0.6}, {0.3,0.3}, {resolution, resolution});

Point<nDim>  diracPoint0(originBZ + 0.25 * sideBZ0 + 0.25 * sideBZ1);
Point<nDim>  diracPoint1(originBZ + 0.75 * sideBZ0 + 0.75 * sideBZ1);

// conduction and valence bands at each Dirac point, from analytic dispersions
graphene.addBand("elC0", -protonCharge, -1.0, meshK0,
    [diracPoint0](Point<nDim> k){ return  std::sqrt(std::pow(k(0)-diracPoint0(0), 2.0)
                                                  + std::pow(k(1)-diracPoint0(1), 2.0)); });
graphene.addBand("elV0", -protonCharge, -1.0, meshK0,
    [diracPoint0](Point<nDim> k){ return -std::sqrt(std::pow(k(0)-diracPoint0(0), 2.0)
                                                  + std::pow(k(1)-diracPoint0(1), 2.0)); });

graphene.plot();
```

---

## Requirements

| | |
|---|---|
| **Compiler** | C++20 (`-std=c++20`), built at `-O3` — the library's design depends on the optimiser |
| **[Eigen](https://eigen.tuxfamily.org)** | Required. Must be on the compiler's header search path |
| **[GNUPlot](http://www.gnuplot.info)** | Required for plotting. Its full path must be set in `TORTOISE/include/Configuration` |

## Getting started

TORTOISE is a library, so you supply the `main` file. Include everything with:

```cpp
#include <TORTOISE>
```

To build your own project against it:

1. Add `TORTOISE/include` and your Eigen directory to the header search paths.
2. Compile every file under `TORTOISE/src` as part of your target.
3. Build with `-std=c++20 -O3`.

A Makefile is provided in `TORTOISE/doc` if you would rather not configure an IDE. Set GNUPlot's path in `TORTOISE/include/Configuration` before running anything that plots.

## Learning the library

- **[`TORTOISE/tutorials`](TORTOISE/tutorials)** — numbered tutorials, in order. Start here.
- **[`TORTOISE/examples/main.cpp`](TORTOISE/examples)** — a complete worked problem.
- **[`TORTOISE/tutorials/developers`](TORTOISE/tutorials/developers)** — internals, for extending the library rather than using it.

## Implementation

Discontinuous Galerkin discretisation over structured and unstructured meshes, adaptive time stepping (Runge–Kutta 4, Dormand–Prince 5(4)), and dedicated quadrature for the high-dimensional collision integrals. The interface is built on expression templates repurposed to produce lazily-evaluated functors, so expressions on fields read as mathematics at the call site without temporaries or runtime abstraction cost.

## Repository layout

```
TORTOISE/
├── include/            Headers — include <TORTOISE> pulls in everything
│   ├── Generics/       Containers, algorithms (RK4, Dormand–Prince, Nelder–Mead), utilities
│   ├── Geometry/       Meshes, regions, function spaces, plotting
│   └── Physics/        Materials, bands, scattering channels, DG integrators
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

## Development version

A version under active development, with further features, is available on request — write to <battiato.marco@gmail.com>.

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).

## Author

Marco Battiato — <battiato.marco@gmail.com>
