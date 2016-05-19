

![lujia-logo](./docs/logo-small.png)

# Lujia-Light

A Julia library for numerical analysis experiments with
multi-scale algorithms for materials science. The package
implements some basic interatomic potentials
and a minimalistic electronic structure model (tight-binding).
It is *not* designed for materials science applications,
but purely for experimentation with new multi-scale algorithms.

**Warning:** as of 19 May, this is under heavy development.

## Installation

In the Julia REPL:
```{.julia}
Pkg.add("PyCall")
Pkg.add("Compose")
Pkg.clone("https://github.com/cortner/LujiaLt.git")
```

## Examples

For now, please see the `IJulia` notebooks in the
`./notebooks` folder.
