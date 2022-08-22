# KineticsCalc

Notebooks and scripts for [CytCKinetics.jl][pkg]. Actual scripts used and data format used for
analysis of kinetic parameters.

## Download and use

If the user has `git` installed, it is possible to clone the repository in the same directory where
`CytCKinetics.jl` is _cloned_. Julia might raise an error upon package instantiation with
`Pkg.instantiate()` because it expects the package to be registered. As a workaround:
+ Remove `CytCKinetics` as a dependency: `pkg> rm CytCKinetics`
+ Add `CytCKinetics` as a dependency through relative path: `pkg> add ../CytCKinetics#v0.3.1`

Since Pluto is a dependency of the this repository, it is not necessary to add it on your own
`@v1.#` environment. Furthermore, the notebooks use the repo environment (`Manifest.toml`) for
uniformity, so it is important to `Pkg.instantiate()` beforehand.

[pkg]: https://github.com/bldamalla/CytCKinetics.jl

