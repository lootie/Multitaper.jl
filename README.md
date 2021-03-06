
![](docs/Multitaper_white.png)

[![Build Status](https://github.com/lootie/Multitaper.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/lootie/Multitaper.jl/actions)
[![Coverage Status](https://coveralls.io/repos/github/lootie/Multitaper.jl/badge.svg?branch=master)](https://coveralls.io/github/lootie/Multitaper.jl?branch=master)

When doing exploratory analysis of time series, frequency domain methods, that is,
statistical methods that display information about the temporal correlations of one
or more time series in terms of frequencies, can be used to infer physical mechanisms
for underlying process dynamics in e.g. geophysical time series, medical time series,
etc.). The [multitaper method](https://en.wikipedia.org/wiki/Multitaper), which
leverages Slepian functions to estimate power spectral densities, coherences, and so
forth, is implemented here for application to univariate, multivariate, and
higher-dimensional (e.g. space-time) processes.

See documentation below. 

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://lootie.github.io/Multitaper.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://lootie.github.io/Multitaper.jl/dev)

## Installation

Please install using 

```

Pkg> add Multitaper

```

This package runs on julia v 1.4.2 and above.

## Paper

If you make use of Multitaper.jl, please cite the following paper: [![DOI](https://joss.theoj.org/papers/10.21105/joss.02463/status.svg)](https://doi.org/10.21105/joss.02463). A previous, unregistered version of Multitaper.jl was hosted on bitbucket.

## Python

Python users may wish to make use of the python wrapper for this package
![multitaperpy](https://github.com/lootie/multitaperpy/). A demo for its useage
is available
![here](https://github.com/lootie/multitaperpy/blob/main/Notebooks/00_GettingStarted.ipynb).

We also highly recommend the multitaper python package by G. A. Prieto, available ![here](https://github.com/gaprieto/multitaper).

## Contributing

We welcome input of any kind via [issues](https://github.com/lootie/Multitaper.jl/issues)
 or by pull requests.
Support requests can be directed to haley@anl.gov.
