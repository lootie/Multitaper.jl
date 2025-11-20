# Similar to mdmultispec, but adaptive weighting has been turned off
"""
   mdmultispec_noadapt(tt, x; <keyword arguments>)

Computes multitaper spectrum for missing-data time series without adaptive weighting. Similar to mdmultispec().

...

# Arguments

- `tt::Vector{T} where T<:Real`: the vector containing the time indices. Must have dt = 1 except across gaps corresponding to missing data.

- `x::Vector{P} where P<:Number`: the vector containing the measurements

- `bw=5/length(tt)`: multitaper bandwidth, default corresponds to NW=5

- `k=Int64(2*bw*size(x,1)-1)`: number of tapers, defaults to 2NW-1

- `lambdau::Union{Tuple{Array{Float64,1}, Array{Float64,2}},Nothing} = nothing`: set to tuple of (spectral concentrations, tapers) if both are precomputed, otherwise function will compute them

- `dt=1`: set if the timestep is something other than 1

- `nz=0`: zero-padding factor. Frequency grid will have length(tt)*(nz+1) entries.

- `Ftest=true`: set to false to turn off the F test for harmonic components

- `jk=true`: set to false to turn of computation of jackknife confidence intervals
...

...

# Outputs

- `MTSpectrum` structure containing the spectrum
...

See also: [`mdmultispec`](@ref), [`mdslepian`](@ref), [`MTSpectrum`](@ref), [`multispec`](@ref)

"""

function mdmultispec_noadapt(tt::Vector{T}, x::Vector{P}; 
                bw=5/length(tt), k=Int64(2*bw*size(x,1)-1), 
                lambdau::Union{Tuple{Array{Float64,1}, Array{Float64,2}},Nothing} = nothing,
                dt=tt[2]-tt[1], nz=0, Ftest=true, jk=true) where{T<:Real,P<:Number}

  lambda,u = (lambdau == nothing) ? Multitaper.mdslepian(bw, k, tt) : lambdau
  s2    = var(x)
  n, nfft, nfft2 = Multitaper._pregap(tt, x, nz)
  ak = Multitaper.multispec_coef(tt, x, u, n, nfft, nfft2)
  sxx   = zeros(size(ak,1))
  d = Multitaper.multispec_aweight(ak, lambda, s2)
  Threads.@threads for j in 1:nfft2
        sxx[j] = Multitaper._dot(d[j,:],abs2.(view(ak,j,:)))
  end
  sxx[2:(nfft2-1)] .*= 2
  coefswts = Multitaper.EigenCoefficient(ak, nothing)
  jv       = jk ? Multitaper.jknife(coefswts,nothing,:spec)[2] : nothing
  # F-test without adaptive weighting
  if Ftest
    dpsw0 = sum(u,dims=1)'
    dpsw02= sum(dpsw0.^2)
    mu    = ak*dpsw0/dpsw02
    num   = (2*k .- 2).*abs2.(mu)*dpsw02
    denom = 2*sum(abs2.(ak - mu*dpsw0'), dims = 2)
    ft    = vec((num./denom))
    Fpval = 1.0 .- cdf.(FDist.(2,2*k .- 2), ft)
  end
  Tv = nothing
  # Package up the outputs
  pkg = Multitaper.MTSpectrum((1/dt)*range(0,1,length=nfft+1)[1:length(sxx)], dt*sxx, nothing, 
              MTParameters(bw*n, k, n, tt[2]-tt[1], 2*(nfft2-1), 1, nothing),
              coefswts, Fpval, jv, Tv)

  return pkg
end
