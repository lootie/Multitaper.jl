---
title: 'Multitaper.jl: A Julia package for frequency domain analysis of time series'
tags:
  - Julia
  - time series
  - statistics
  - spatial statistics
  - space-time processes
authors:
  - name: Charlotte L. Haley
    orcid: 0000-0003-3996-773X
    affiliation: "1"
    email: "haley@anl.gov" 
  - name: Christopher J. Geoga
    affiliation: "1, 2"
affiliations:
 - name: Argonne National Laboratory
   index: 1
 - name: Rutgers University
   index: 2
date: 2 July, 2020
bibliography: paper.bib

---

# Summary

Spectral analysis is widely used to infer physical mechanisms for underlying
process dynamics from a realization of a stationary time series. The multitaper
method is a nonparametric technique for estimating the power spectrum of a discrete
time series that simultaneously controls bias and variance of the estimator by
premultiplying the data by a set of orthogonal sequences---discrete prolate
spheroidal sequences that are optimally concentrated in both time and frequency.
These have the effect of stabilizing the bias and variance of the spectrum estimator
[@T82]. While multitaper codes have been introduced in multiple languages, including
Julia, the `Multitaper.jl` package that we present offers functionality beyond
univariate and bivariate time series analysis and provides routines for a number of
(selected) research-level topics not found in other packages.

# Statement of need

`Multitaper.jl` is a Julia package for spectrum analysis of time series, multivariate
time series, and spatial or space-time processes. `Multitaper.jl` was designed to be
useful to researchers in diverse fields, including geophysics (climate, seismology,
limnology, and stratigraphy), cognitive radio, space science (solar physics), speech
processing, astronomy, and biomedicine. For example, a researcher might
want to compute the multitaper spectrum of a time series so he or she can identify
which periodic components contribute the most to signal variance, and do so with
jackknifed error bounds on the oscillation amplitudes. 

The high-level character of Julia allows for widely readable and extendible codes,
while the low-level functionality provides speed and efficiency. The `Multitaper.jl`
package provides a user-friendly implementation of many of the basic concepts such as
spectrum analysis, F-testing for harmonic analysis [@T82], coherence and phase,
jackknifed variance estimates [@TC91], and complex demodulation [@T07]; more advanced
techniques such as dual-frequency spectra, cepstrum, multitaper for time series
containing gaps [@chave2019multitaper], $T^2$ tests for multiple line components
[@T11] implementations of higher-dimensional Slepian tapers on Cartesian domains
[@SimonsWang2011] [@Geoga2018]; and others [@ThomsonHaley2014] [@HaleyAnitescu2017].
In addition, we provide tutorial-style notebooks to allow accessibility to those new
to these concepts or to Julia in general.

`Multitaper.jl` has been used in graduate courses to provide fast spectrum estimates
of unequally sampled time series. Early versions of this code have also been used to
compute statistics for research publications.

# Other software

Multitaper spectrum analysis is implemented in Julia in the `DSP.jl` package, but
it is limited to estimation of the spectrum. In the R programming language the `R
multitaper` package gives an R-wrapped Fortran 77 implementation, which provides
fast spectrum estimates, F-tests, jackknifing, coherences, and complex demdodulation
[@Rahim]. A C-subroutine for multitaper methods was introduced in
[@lees1995multiple], and the derivative Python implementation `pymutt` [@Smith] is
a wrapped version of the former.  The Fortran 90 library of multitaper methods
[@Prieto] provides for spectrum analysis, F-testing, spectral reshaping, coherences,
and quadratic inverse spectrum estimation. While wrapped versions of low-level codes
run rapidly, they can be more difficult to extend by the research community. The
Matlab Signal Processing Toolbox implements discrete prolate spheroidal sequences and
multitaper spectrum estimation, while standalone contributions on the Matlab file exchange
describe the extension to multitaper methods with gaps [@chave2019multitaper],
higher-dimensional spectrum estimation on Cartesian domains [@SimonsWang2011] and the
sphere [@simons2006] and the freely available `jlab` codes [@Lilly].

# To contribute

We welcome input of any kind via bitbucket issues or by pull requests.
Support requests can be directed to haley@anl.gov.

# Acknowledgements

We acknowledge contributions from Mihai Anitescu, David J. Thomson, and
Sally Dodson-Robinson during the writing of these codes. We also gratefully
acknowledge the help of our reviewers, Andy Nowacki and Robert Guggenberger, in 
editing the code repository.

This work was supported by the U.S. Department of Energy, Office of Science, Advanced
Scientific Computing Research, under contract number DE-AC02-06CH11357.

The submitted manuscript has been created by UChicago Argonne, LLC, Operator of
Argonne National Laboratory (“Argonne”). Argonne, a U.S. Department of Energy Office
of Science laboratory, is operated under Contract No. DE-AC02-06CH11357. The U.S.
Government retains for itself, and others acting on its behalf, a paid-up
nonexclusive, irrevocable worldwide license in said article to reproduce, prepare
derivative works, distribute copies to the public, and perform publicly and display
publicly, by or on behalf of the Government. The Department of Energy will provide
public access to these results of federally sponsored research in accordance with the
DOE Public Access Plan. http://energy.gov/downloads/doe-public-access-plan

# References
