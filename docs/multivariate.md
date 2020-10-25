
## Multivariate time series analysis in Multitaper.jl

The multivariate capabilities of this package include cross-spectrum and coherence
analysis, which contain information about the covariance properties (respectively,
correlations) between two time series. As in the univariate case, multiple tapers are
used to de-bias the estimates, and this results in individual estimates that can be
used to get a jackknifed confidence interval for either quantity of interest.
Trivially, one can obtain multitaper estimates for the cross-covariance and
cross-correlation functions by inverse Fourier transforming the multitaper
cross-spectrum and coherence. 

## Quick Synopsis of Capabilities

Multivariate

* Magnitude squared coherence  

* Jackknife estimates of phase

## Multivariate spectrum estimation 

### multispec

The multispec command does a lot of the multivariate stuff as well. You will
recognize most of the keyword arguments from the univariate version of this call.
However, there are a few differences now. You can either use the first, simpler
version, of the multivariate call: 

```
function multispec(S1::Union{Vector{T},Ecoef}, S2::Union{Vector{T},Ecoef}; 
                   outp=:coh, NW=4.0, K=6, offset=0, dt=1.0, ctr=true, pad=1.0,
                   dpVec=nothing, guts=false, jk=false, Tsq=nothing, alph=0.05) where{T}
```

In which

 * `S1` is now the first input time series, and S2 is the second, OR, S1 and S2 can be
  eigencoefficient structs from a previous computation (for speed if you are calling
these functions often) 

  * `outp` is the desired output of the computation which now can be `:spec` for
  cross-spectrum, `:coh` for coherence and phase, or `:transf` for transfer
function.
  * `NW`, the choice of time-bandwidth product
  * `K`, the number of tapers to use
  * `offset`, the frequency offset, if desired (else 0.0)
  * `dt`, the temporal sampling frequency (in, say, seconds)
  * `ctr`, whether or not to remove the mean from the data series, default is true
  * `pad`, default is not to pad, but if this is a float then the padded length will be
    pad times the length of S1 and if it is an int greater than the length of S1 that
    will be the FFT length. 
  * `dpVec`, you can choose to supply the dpss's or not (speeds things up if you're
    calling the function many times)
  * `guts`, whether you'd like the eigencoefficients as output. They will come in an
    eigencoefficient struct with the field coef, wts where coef contains the
    eigencoefficients and wts will contain the adaptive weights
  * `jk`, jackknifing to give a confidence interval for the spectrum
  * `Tsq`, T-squared test for multiple line components (Thomson Asilomar conference
    proceedings)
  * `alph`, confidence level for jackknife confidence intervals and Tsq tests.

The output struct you get will be determined by what `outp` is, namely one of 

- `MtSpec` struct which will contain the following fields (in
the following order):

  * frequency (f), as a LinRange
  * cross-spectrum (S), a vector giving half the spectrum up to the Nyquist if the input is
    real
  * phase, 
  * chosen values of the multitaper time bandwidth product etc of type MtParams
    (params) This makes its own parameter struct that contains NW, K, N, dt, M (padded
  length), nsegments (number of segments of data to averae), overlap (if the sample was
  divided into overlapping chunks) and it gets carried around for future reference and
  for plotting purposes
  * eigencoefficients (coef, optional), 
  * Ftest values (Fpval, optional), 
  * jackknife output (jkvar, optional), and
  * Tsquared test results (Tsq_pval, optional).This struct was described
  earlier, in this case the `phase` output field will be filled in.

* `MtCoh` coherence struct. Its fields are 

- frequency (f)
- coh, a vector giving the squared coherence up to the Nyquist if the input is real
- phase, 
- chosen values of the multitaper time bandwidth product etc of type MtParams
  (params) This makes its own parameter struct that contains NW, K, N, dt, M (padded
  length), nsegments (number of segments of data to averae), overlap (if the sample 
  was divided into overlapping chunks) and it gets carried around for future reference 
  and for plotting purposes
- eigencoefficients (coef, optional), 
- jackknife output (jkvar, optional), and
- Tsquared test results (Tsq_pval, optional). 

* `MtTransf` transfer function struct. Its fields are
- frequency (f)
- transf, a vector giving the transfer function up to the Nyquist if the input is
  real
- phase, 
- chosen values of the multitaper time bandwidth product etc of type MtParams
  (params) This makes its own parameter struct that contains NW, K, N, dt, M (padded
length), nsegments (number of segments of data to averae), overlap (if the sample was
divided into overlapping chunks) and it gets carried around for future reference and
for plotting purposes
- eigencoefficients (coef, optional), 
- jackknife output (jkvar, optional), and
- Tsquared test results (Tsq_pval, optional). 

Now, you can also use the batch-version of this call: 

```
function multispec(S1::Matrix{T}; outp=:coh, NW=4.0, K=6, dt=1.0, ctr=true,
                   pad=1.0, dpVec=nothing, guts=false, a_weight=true, jk=false,
                   Ftest=false, Tsq=nothing, alph=0.05) where{T}
```

  * `S1` is now just a matrix with, say, p columns and N rows, meaning that there are p
  input time series. 
  * `outp` is the desired output of the computation which now can be `:cross` for
  cross-spectrum, `:coh` for coherence and phase, or `:justspecs` for only the
  spectra.
  * `NW`, the choice of time-bandwidth product
  * `K`, the number of tapers to use
  * `offset`, the frequency offset, if desired (else 0.0)
  * `dt`, the temporal sampling frequency (in, say, seconds)
  * `ctr`, whether or not to remove the mean from the data series, default is true
  * `pad`, default is not to pad, but if this is a float then the padded length will be
    pad times the length of S1 and if it is an int greater than the length of S1 that
    will be the FFT length. 
  * `dpVec`, you can choose to supply the dpss's or not (speeds things up if you're
    calling the function many times)
  * `guts`, whether you'd like the eigencoefficients as output. They will come in an
    eigencoefficient struct with the field coef, wts where coef contains the
    eigencoefficients and wts will contain the adaptive weights
  * `jk`, jackknifing to give a confidence interval for the spectrum
  * `Tsq`, T-squared test for multiple line components (Thomson Asilomar conference
    proceedings)
  * `alph`, confidence level for jackknife confidence intervals and Tsq tests.

The output of this command, depending on the desired output type, is one of three
things:

* `outp = :justspecs` is a vector of MtSpec structs containing the spectra alone. 

* `outp = :coh` is a tuple containing the spectra, a matrix filled with coherences on
  the super-diagonal (so if the result was called out, you'd access it by using
`out[2][1,2]` to get the coherence between the first and second series.), and finally
the result of the T-squared test, if you asked for it. 

* `outp = :cross` is a tuple containing the spectra, a matrix filled with cross
  spectra on the super-diagonal, and finally the result of the T-squared test, if you
asked for it. 

A note on plotting: if you are using Plots.jl there recipes to directly plot the
output of the above multivariate calculations, especially the tuples, in a gridded
plot format. See the second jupyter notebook for more details. 

# Missing-data coherences

Extending the missing-data spectrum estimation of Chave from the univariate case to
the bivariate case, one can compute coherences using the function with signature

```
function mdmultispec(t::Union{Vector{Int64}, Vector{Float64}}, 
                x::Vector{Float64},
                y::Vector{Float64};
                bw::Float64 = 5/length(t),
                k::Int64    = Int64(2*bw*size(x,1) - 1),
                dt::Float64 = 1.0, jk::Bool = true,
                nz::Union{Int64,Float64}   = 0, 
                Ftest::Bool = false,
                lambdau::Union{Tuple{Array{Float64,1},
                               Array{Float64,2}},Nothing} = nothing)
```

The inputs are the following:
  * t -- real vector of time
  * x -- first missing-data time series
  * y -- second missing-data time series
  * bw -- bandwidth of estimate, 5/length(t) default
  * k -- number of slepian tapers, must be <=2 bw length(x), 2 bw length(x)-1 default
  * dt -- sampling in time
  * jk -- whether or not to compute jackknife variance estimates
  * nz -- zero padding factor, 0 default
  * Ftest -- whether or not to compute the F-test p-value at all frequencies
  * lambdau -- missing data Slepian tapers and their concentrations, if precomputed

The output is 
  * sxx -- MtCoh coherence estimate 

## Time domain:

# mt_ccvf

This function computes multitaper estimates of the cross-covariance and
cross-correlation by way of inverse-FFT of a multitaper spectrum estimate. Its
signature is either

``` function mt_ccvf(S::MtSpec; typ::Symbol = :ccvf) ```

or

```
mt_ccvf(S1::Vector{T}, S2::Vector{T}; typ::Symbol = :ccvf, NW::Real = 4.0, K::Int = 6, 
                   dt::Float64=1.0, ctr::Bool = true, 
                   pad::Union{Int,Float64} = 1.0, 
                   dpVec::Union{Vector{Float64},Matrix{Float64},Nothing} = nothing,
                   guts::Bool = false, 
                   jk::Bool = false, 
                   Tsq::Union{Vector{Float64},Vector{Vector{Float64}},Vector{Int64},
                   Vector{Vector{Int64}},Nothing}=nothing, 
                   alph::Float64 = 0.05) where T<:Number
```

In the first, we assume that you have the MtSpec struct handy (must be cross-spectra,
coherences won't work), and in the second you give the two time series, similar to
above.  The `typ` kwarg can take values in (`:ccvf` and `:ccf`) with `:ccvf` being
the default value. Depending on the value of typ, you will get one of two different
structs

* MtCcf: Contains lags, cross correlation function, and a params struct (mentioned
  above) that carries around the relevant multitaper options. 

* MtCcvf: Contains lags, cross covariance function, and a params struct.

when you plot one of the `MtCcf` or `MtCcvf` structs using the recipe, you'll get a
stem plot. 



