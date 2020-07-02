
## Multivariate time series analysis in Multitaper.jl

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
                    outp::Symbol = :coh, ...) where T<:Number
```

In which

- S1 is now the first input time series, and S2 is the second, OR, S1 and S2 can be
  eigencoefficient structs from a previous computation (for speed if you are calling
these functions often) 

- The ellipsis replaces all of the keyword arguments from the univariate version of
  multispec (already discussed in the other doc)

- `outp` is the desired output of the computation which now can be `:spec` for
  cross-spectrum, `:coh` for coherence and phase, or `:transf` for transfer
function.

The output struct you get will be determined by what `outp` is, namely one of 

* `MtSpec` struct which contains the cross spectrum. This struct was described
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
function multispec(S1::Matrix{T}; outp::Symbol=:coh, 
                   ...) where T<:Number
```

- S1 is now just a matrix with, say, p columns and N rows, meaning that there are p
  input time series. 

- The ellipsis replaces all of the keyword arguments from the univariate version of
  multispec (already discussed in the other doc)

- `outp` is the desired output of the computation which now can be `:cross` for
  cross-spectrum, `:coh` for coherence and phase, or `:justspecs` for only the
spectra.

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

## Time domain:

# mt_ccvf

This function computes multitaper estimates of the cross-covariance and
cross-correlation by way of inverse-FFT of a multitaper spectrum estimate. Its
signature is either

``` function mt_acvf(S::MtSpec; typ::Symbol = :ccvf) ```

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



