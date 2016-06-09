

![lujia-logo](./docs/logo-small.png)

# Lujia-Light

A Julia library for numerical analysis experiments with
multi-scale algorithms for materials science, based on an ongoing work on
a book manuscript. The package implements some basic interatomic potentials
and a minimalistic electronic structure model (tight-binding).
It is *not* designed for materials science applications,
but purely for experimentation with new multi-scale algorithms
in a simplified setting.

For serious molecular simulation see
[Atoms.jl](https://github.com/libAtoms/Atoms.jl), but note that
this is also still under heavy development.

**Warning:** This package is under heavy development. At the
moment, the following notebooks can be followed to look at some
completed parts of the library:

* `Atm Examples`: a few examples running a basic atomistic supercell model
         with clamped boundary conditions
* `FEM`: brief intro how to use the `LujiaLt.FEM` module
* `Introduction to Approximation`: implementation of some elementary toy problems
* `3DGraphics`: collection of code snippets to generate some visualisations
         used in a book manuscript on which this package is based.

## Installation

In the Julia REPL:
```jl
Pkg.add("PyCall")
Pkg.add("Compose")
Pkg.add("PGFPlots")
Pkg.clone("https://github.com/cortner/LujiaLt.git")
```

## Quick Example

Copy-paste the following into an IJulia notebook (plotting doesn't
work from the REPL at the moment)

```jl
using LujiaLt
using LujiaLt.Potentials, LujiaLt.Plotting, LujiaLt.Solve
at = Atm(V=ToyEAMPotential(), Ra=8.1, defect=:vacancy)
Ysol = solve(at, display_result = true)
plot(at, X=Ysol, axis = [-5.1, 5.1, -3.1, 3.1])
```

<!-- ## Examples

For now, please see the `IJulia` notebooks in the
`./notebooks` folder. -->
