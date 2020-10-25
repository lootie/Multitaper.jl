
## Univariate time series functionality for Multitaper.jl

The univariate functionality of this package computes multitaper estimates of the
power spectrum, autocorrelation, cepstrum, and complex demodulation of a time series.
One can also estimate the power spectrum when there is missing data.  

An estimate of the power spectral density contains information about the sample
variance of the time series, decomposed in terms of frequency. When the multitaper
method is used, multiple tapering functions ensure that the statistical bias-variance
tradeoff is optimized.  Discrete prolate spheroidal sequences, which are used as
tapering functions, are finite-length sequences whose Fourier transforms contain most
of their power on a narrow bandwidth. This property makes them optimal for reducing
spectral leakage. The orthogonality of the tapering functions also means that one
can compute empirical confidence intervals for the spectrum estimate by jackknifing,
and estimate the significance of discrete sinusoids using an F-test.  Finally, one
can compute the quantities above when time series are sampled on a regular grid for
which there are missing data points. 

Cepstrum analysis is widely used in speech processing applications and for
deconvolution. Complex demodulation of a time series allows for examination of the
amplitude and phase of a sinusoid at a frequency of interest to the time series being
studied. 

## Quick Synopsis of Capabilities

Univariate

* Discrete Prolate Spheroidal sequences

* Multitaper spectra that use dpss (multispec) tapers

* Jackknifing

## Univariate Multitaper Estimation

### multispec

This command does a lot of the univariate stuff. Its signature is:

```
function multispec(S1::Union{Vector{T}}; 
                   NW::Real = 4.0, 
                   K::Int = 6, 
                   dt::Float64 = 1.0, 
                   ctr::Bool = true, 
                   pad::Union{Int,Float64} = 1.0, 
                   dpVec::Union{Vector{Float64},Matrix{Float64},Nothing} = nothing,
                   egval::Union{Vector{Float64},Nothing} = nothing,
                   guts::Bool = false, 
                   a_weight::Bool = true, 
                   Ftest::Bool = false, 
                   highres::Bool = false,
                   jk::Bool = false, 
                   Tsq::Union{Vector{Float64},Vector{Vector{Float64}},
                      Vector{Int64},Vector{Vector{Int64}},Nothing} = nothing, 
                   alph::Float64 = 0.05
                  ) where T<:Number
```
S1 is the time series, and you have the following keyword options:

  * NW, the choice of time-bandwidth product
  * K, the number of tapers to use
  * dt, the temporal sampling frequency (in, say, seconds)
  * ctr, whether or not to remove the mean from the data series, default is true
  * pad, default is not to pad, but if this is a float then the padded length will be
    pad times the length of S1 and if it is an int greater than the length of S1 that
    will be the FFT length. 
  * dpVec, you can choose to supply the dpss's or not (speeds things up if you're
    calling the function many times)
  * egval, eigenvalues associated with the vectors above, optional
  * guts, whether you'd like the eigencoefficients as output. They will come in an
    eigencoefficient struct with the field coef, wts where coef contains the
    eigencoefficients and wts will contain the adaptive weights
  * a_weight, whether to use the adaptive weighting scheme described in Thomson, 1982
  * Ftest, whether to compute the harmonic F-test described in Thomson, 1982
  * highres, whether the estimate should be a high-resolution estimate, see Thomson 1982
  * jk, jackknifing to give a confidence interval for the spectrum
  * Tsq, T-squared test for multiple line components (Thomson Asilomar conference
    proceedings)
  * alph, confidence level for jackknife confidence intervals and Tsq tests.

The output of this command is a MtSpec struct which contains the following fields (in
the following order):

  * frequency (f), as a LinRange
  * spectrum (S), a vector giving half the spectrum up to the Nyquist if the input is
    real
  * phase (optional), 
  * chosen values of the multitaper time bandwidth product etc of type MtParams
    (params) This makes its own parameter struct that contains NW, K, N, dt, M (padded
  length), nsegments (number of segments of data to averae), overlap (if the sample was
  divided into overlapping chunks) and it gets carried around for future reference and
  for plotting purposes
  * eigencoefficients (coef, optional), 
  * Ftest values (Fpval, optional), 
  * jackknife output (jkvar, optional), and
  * Tsquared test results (Tsq_pval, optional). 

Some of the fields will contain nothing, or if you requested one of the jacknifed
confidence interval (`jk = true`), F-test for line components (`Ftest = true`),
multivariate T^2 test for line components (`Tsq != nothing`), or the
eigencoefficients and possibly adaptive weights (`guts=true`), then the relevant
fields in the output struct will be filled. 

Note that `pad` here can be given in two different ways: if it is a float, then it is
the factor by which to multiply the length of the original data sequence, and if it
is an integer larger than the length of the original data sequence, it will be the
number of data points in the full FFT.

If dpVec is given, you have supplied pre-computed Slepians, which will speed things
up if the function is going to be called many times.  The option `a_weight` uses
adaptive weighting. 

A note on plotting: if you are using Plots.jl there are pre-loaded recipes that make
plotting of MtSpec structs completely trivial. Simply plot your MtSpec structs as if
they were vectors, and you'll get a bunch of preformatting for free. Consult the
jupyter notebooks for examples of some of the recipes.

# Welch

The welch estimate of the spectrum is one that is an average of multitaper spectra
computed on overlapping data blocks. You'd call the function like this

```
welch(S1::Union{Vector{Float64}, Vector{ComplexF64}, Matrix{Float64}, 
               Matrix{ComplexF64}}, 
               nsegments::Int64, overlap::Float64 = 0.5, ws::Symbol = :welch; 
               outp::Symbol = :spec,
               NW::Real = 4.0, K::Int = 6, dt::Float64 = 1.0,
               ctr::Bool = true, pad::Union{Int,Float64} = 1.0,
               dpVec::Union{Vector{Float64},Matrix{Float64},Nothing} = nothing, 
               egval::Union{Vector{Float64},Nothing} = nothing, 
               guts::Bool = false, a_weight::Bool = true, Ftest::Bool = false, 
               jk::Bool = false, 
               Tsq::Union{Array{Int64,1},Array{Array{Int64,1},1},Nothing}=nothing,
               alph::Float64 = 0.05) 
```

The call is very similar to the multispec call above, except you will enter the
parameters nsegments and overlap, that is, the number of data blocks to make, and the
overlap between them. Since there is some fudging these numbers, you will get the
number of segments and overlap back in the params field of the output struct. The ws
toggle is by default set to :welch mode, but this code is also the backbone for a
spectrogram, which is why this keyword argument is here. 

## Time domain stuff

# mt_acvf

This function computes multitaper estimates of the covariance, correlation, and
cepstrum, by way of inverse-FFT of a multitaper spectrum estimate. Its signature is
either

``` 
function mt_acvf(S::MtSpec; typ::Symbol = :acvf)
```

or

```
function mt_acvf(S1::Union{Vector{T}}; 
                 typ::Symbol = :acvf, 
                 NW::Real = 4.0, 
                 K::Int = 6, 
                 dt::Float64=1.0, 
                 ctr::Bool = true, 
                 pad::Union{Int,Float64} = 1.0, 
                 dpVec::Union{Vector{Float64},Matrix{Float64},Nothing} = nothing,
                 egval::Union{Vector{Float64},Nothing} = nothing,
                 a_weight::Bool = true, 
                 reshape::Union{Bool, Float64, Int64, Vector{Float64}, 
                          Vector{Int64}} = false,
                 ) where T<:Number
```

The first method is used if you've already computed the spectrum and have the MtSpec
struct handy.  If not, you can use the second version, putting in the time series,
and using any other relevant input arguments mentioned in the multispec call above.


The only other toggle is `typ` which can take values in (`:acvf`, `:acf`, and
`:ceps`) with `:acvf` being the default value. Depending on the value of typ, you
will get one of three different structs

  * MtAcf: Contains lags, autocorrelation function, and a params struct (mentioned
    above) that carries around the relevant multitaper options. 

  * MtAcvf: Contains lags, autocovariance function, and a params struct.

  * MtCeps: Contains lags (quefrency), a cepstrum estimate, and a params struct. The
    cepstrum is the inverse-FFT (or cosine transform, when the signal is real) of the
  logarithm of the spectrum. 

when you plot one of the `MtAcf`, `MtAcvf`, or `MtCeps` structs using the recipe,
you'll get a stem plot. 

# demodulate

This function computes the multitaper estimate of the complex demodulate. Note that
there are several published implementations of this, but this one uses a single
zeroth order slepian taper as a filter. Note that this is identical to the
implementation in the R multitaper package. The function signature is 

```
function demodulate(x::Vector{Float64}, f0::Float64, NW::Float64, blockLen::Int64, 
                    wrapphase::Bool = true, dt::Float64, basetime::Float64)
```

The inputs are the following:
  * x is the data vector to be demodulated
  * dt is the sampling rate (in years, for example)
  * f0 is the center frequency (e.g. one cycle per year)
  * NW is the time-bandwidth product for the filter (typically narrow, use a float)
  * blockLen is a subjective length to use for the data filter, note that the output 
    will be shortened
    by this amount, so a short filter is sometimes better. 
  * wrapphase is an optional boolean which tells the algorithm whether or not to 
    unwrap the phase for pretty plotting. 
  * basetime is the time for the first time index

The output will be a `Demodulate` struct containing time (time), magnitude (mag), and
phase (phase).  For easy plotting there is a Plots.jl recipe. See also the notebook
in the Examples directory for usage.  

# Data with gaps

## Spectrum estimation for data with gaps

Chave wrote a paper in GJI that shows how to compute dpss's on data with
gaps and we have rewritten the code in Julia so that it can be used here. The
function signature is 

```
function mdmultispec(tt::Union{Vector{Int64},Vector{Float64}}, x::Vector{Float64}; 
                bw=5/length(tt), k=Int64(2*bw*size(x,1)-1), 
                lambdau::Union{Tuple{Array{Float64,1},
                               Array{Float64,2}},Nothing} = nothing,
                dt=tt[2]-tt[1], nz=0, Ftest=true, jk=true,
                dof=false)
```

The inputs are the following:
  * tt -- real vector of time (required)
  * x -- real vector of data (required)
  * bw -- bandwidth of estimate, 5/length(t) default
  * k -- number of slepian tapers, must be <=2 bw length(x), 2 bw length(x)-1 default
  * lambdau -- missing data Slepian tapers and their concentrations, if precomputed
  * dt -- sampling in time
  * nz -- zero padding factor, 0 default
  * Ftest -- whether or not to compute the F-test p-value at all frequencies
  * jk -- whether or not to compute jackknife variance estimates
  * dof -- whether or not to output the degrees of freedom of the estimate

The outputs are simply a list of the following 
  * sxx -- MtSpec spectrum 
  * nu1 -- Degrees of freedom, if dof is set to true

## Missing-data Slepian tapers

One can also simply generate Slepian tapers using the mdslepian function. Its
function signature is 

```
function mdslepian(w, k, t)
``` 

The inputs are the following:
  * w -- the bandwidth of the taper
  * k -- the number of Slepian tapers
  * t -- the vector containing the time indicees

The outputs are
  * lambda -- the concentrations of the missing-data Slepians
  * u -- length(t) times k matrix containing the missing-data Slepians

## Generalized prolate spheroidal sequences

Following the work of Bronez in 1988, one can compute optimally concentrated data
tapers on a general, uneven, temporal grid. The function with the following signature
computes these

```
function gpss(w::Float64, k::Int64, t::Union{Vector{Int64},Vector{Float64}}, 
        f::Float64; beta::Float64 = 0.5)
```

The inputs are the following:
  * w -- the bandwidth of the taper
  * k -- the number of Slepian tapers
  * t -- the vector containing the time indicees
  * f -- the frequency of interest, between 0 and beta
  * beta -- the unequal sampling equivalent to the Nyquist rate

The outputs are
  * lambda -- the concentrations of the missing-data Slepians
  * u -- length(t) times k matrix containing the missing-data Slepians
  * R -- Cholesky factor for the generalized eigenvalue problem

Note that this function is unexported, so one needs to use `Multitaper.gpss(...)` to
call it. This function is known to occasionally have numerical errors in computing
the Cholesky factors of the genealized eigenvalue problem, so use at own risk. For
the cases we tested, however, this function generally returns the expected result.


If you make use of the functions in this "Data with gaps" section, kindly
cite: Chave, Alan D. "A multitaper spectral estimator for
time-series with missing data." Geophysical Journal International 218.3 (2019):
2165-2178.


