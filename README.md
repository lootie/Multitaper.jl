
# Multitaper.jl

## Installation

This package is unregistered, so please install with

```

Pkg> add git://bitbucket.org/clhaley/Multitaper.jl.git

```

## Quick Start

If you want a quick-and-dirty univariate spectrum estimate, simply issue

``` 

julia> using Multitaper 

julia> S = multispec(x) 

```

at the prompt, where x is a real data vector. The default time-bandwidth product selection is `NW =
4.0` while the number of tapers `K = 6`. If Plots.jl is being used, the resulting struct can be
plotted via the recipe

```

julia> using Plots

julia> plot(S)

```  

For further in-depth reading, consult the docs. For useful examples, see the notebooks in the
Examples folder (requires IJulia.jl).

## Quick Synopsis of Capabilities

As of version 0.0, 01/2020

# Univariate

* Discrete Prolate Spheroidal sequences and sine tapers

* Multitaper spectra that use dpss (multispec) tapers. 

* Jackknifing for estimation of the confidence intervals.

* F-test for line components (Thomson, 1982).

* Expected number of upcrossings of extreme levels of a univariate spectrum estimate on white noise
  data. See Thomson and Haley, 2014. 

# Multivariate

* Magnitude squared coherence, implemented in the spirit of Thomson and Chave, 1991.

* Jackknife estimates of phase, unwrapped, similar to R's implementation. 

* T-squared test (Thomson, "Some comments on possibly cyclostationary series", Asilomar Conference
  Proceedings) for simultaneous line components. 

## Plots Recipes

You will find plots recipes in the separate file ./Examples/plotsrecipes.jl. These provide a
convenient Plots.jl format for plotting multispec objects.

## Main References 

The functions described below are those mainly given in 

@Article{Thomson82,
  author        = {D. J. Thomson},
  journal       = {Proceedings of the IEEE},
  pages         = {1055--1096},
  title         = {Spectrum estimation and harmonic analysis},
  volume        = {70},
  number        = {9},
  year          = {1982}
}  

@InCollection{ThomsonChave91,
  address       = {Upper Saddle River, NJ},
  author        = {D. J. Thomson and A. D. Chave},
  booktitle     = {Advances in Spectrum Analysis and Array Processing},
  chapter       = {2},
  editor        = {Simon Haykin},
  pages         = {58--113},
  publisher     = {Prentice-Hall},
  title         = {Jackknifed error estimates for spectra, coherences, and
                  transfer functions},
  volume        = {1},
  year          = {1991}
}

## License

This software is distributed under the GNU GPL v2 license.


