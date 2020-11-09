# Multitaper.jl Documentation

*Documentation for Multitaper.jl*

A package for multitaper frequency-domain analysis.

## Package Features

When doing exploratory analysis of time series, frequency domain methods, that is,
statistical methods that display information about the temporal correlations of one
or more time series in terms of frequencies, can be used to infer physical mechanisms
for underlying process dynamics in e.g. geophysical time series, medical time series,
etc.). The [multitaper method](https://en.wikipedia.org/wiki/Multitaper), which
leverages Slepian functions to estimate power spectral densities, coherences, and so
forth, is implemented here for application to univariate, multivariate, and
higher-dimensional (e.g. space-time) processes.

### Univariate time series

* Discrete Prolate Spheroidal sequences (Slepian 1974) and with gaps (Chave 2019) 

* Multitaper spectra that use dpss (multispec) tapers, where time series can have
  gaps. 

* Jackknifing for estimation of the confidence intervals.

* F-test for line components (Thomson, 1982).

* Complex demodulation

* Multitaper spectrum estimation for time series with equal temporal spacing except
  with gaps (Chave 2019)

### Multivariate time series

* Magnitude squared coherence, implemented in the spirit of Thomson and Chave, 1991.

* Jackknife estimates of phase, unwrapped, similar to R's implementation. 

* T-squared test (Thomson, "Some comments on possibly cyclostationary series",
  Asilomar Conference Proceedings) for simultaneous line components. 

### 2 Dimensional

* 2 Dimensional Cartesian tapers supported on the square in space and 2D disk in time
  as described in Simons and Wang, 2011 

### Plots Recipes

Plots recipes are provided for numerous function ouputs. 

The [Examples](@ref) provides tutorial IJulia.jl notebooks showing how to get 
started using Multitaper.

See the [Index](@ref main-index) for the complete list of documented functions and
types.

## Installation

This package is unregistered, so please install with

```
Pkg> add https://bitbucket.org/clhaley/Multitaper.jl.git
```

This package runs on julia v 1.4.2 and above.

## Testing

In order to test this package, use

```@julia-repl
Pkg> test Multitaper
```

This will run the tests in `test/runtests.jl`.

## Manual

```@contents
```

### [Index](@id main-index)

```@index
Pages = ["lib/public.md"]
```
