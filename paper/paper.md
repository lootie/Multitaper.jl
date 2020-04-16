---
title: 'Multitaper.jl: A Julia package for power spectrum analysis'
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
  - name: Christopher J. Geoga
    affiliation: "1, 2"
affiliations:
 - name: Argonne National Laboratory
   index: 1
 - name: Rutgers University
   index: 2
date: 16 April, 2020
bibliography: paper.bib

---

# Summary

Spectral analysis is widely used to infer physical mechanisms for underlying process dynamics from a
realization of a stationary time series.  The multitaper method is a nonparametric technique for the
estimation of the Fourier spectrum of a discrete time series which simultaneously controls bias and
variance of the estimator by premultiplying the data by a set of orthogonal sequences, discrete
prolate spheroidal sequences, having optimality properties chosen by adjusting a bandwidth
parameter [@T82]. While codes are readily available in multiple languages for time-series analysis, spatial
and higher dimensional processes require numerical integration on fine grids that generally cannot
be accommodated for by higher level languages. 

`Multitaper.jl` is a Julia package for spectrum analysis of time series, multivariate time series,
and space-time processes. The high-level character of Julia allows for widely readable and
extendible codes while the low-level functionality provides speed and efficiency. The
`Multitaper.jl` package was built to provide a user-friendly implementation of many of the basic
concepts such as spectrum analysis, coherence and phase, dual-frequency spectra, cepstrum, and
others, as well as tutorial notebooks on research topics such as crossings and multiple comparisons
[@ThomsonHaley2014], and finally implementations of higher dimensional spectrum analysis on Cartesian
domains and the sphere [@SimonsWang2011].  

`Multitaper.jl` was designed to be useful both by statisticians and researchears in diverse fields,
namely geophysics (climate, seismology, limnology, and stratigraphy), cognitive radio,
space science (solar physics), speech processing, astronomy, and biomedicine. It has been used in
graduate courses in astronomy to provide fast spectrum estimates of unequally sampled time series.
Early versions of this code have been used to compute figures for research publications [@Geoga2018].  

# Acknowledgements

We acknowledge contributions from Mihai Anitescu, David J. Thomson, and Alan Chave during 
the genesis of this project.

# References
