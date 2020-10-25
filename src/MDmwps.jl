# MT routines by Chave, adapted for Julia from Matlab by C Haley in 2020 Original
# Credit is due to 
#
# Chave, Alan D. "A multitaper spectral estimator for time-series with missing data."
# Geophysical Journal International 218.3 (2019): 2165-2178.
#
# Coherence routine is original
#

using Statistics, Arpack, LinearAlgebra, Distributions, NFFT

"""
Find all of the necessary lengths
"""
function _pregap(tt, xx, nz)
  n    = length(tt)
  nfft = length(xx)
  if mod(nfft,2) != 0 
    nfft = nfft + 1
  end
  nfft  = Int64(round((nz + 1)*nfft))
  nfft2 = Int64(round(nfft/2)) + 1
  return n, nfft, nfft2 
end

""" Computes multitaper adaptive weighting """
function multispec_aweight(eigcoefs, egval, cvar, maxit = 15, tol = 0.05)
  halffreq,k   = size(eigcoefs)
  dsq = (egval == nothing) ? nothing : zeros(size(eigcoefs))
  if egval != nothing
    (old, new) = vcat(ones(2).*0.5, zeros(k-2)), zeros(k)     
      Threads.@threads for j in 1:halffreq
      @inbounds dsq[j,:] = aweighted(old, new, k, abs2.(view(eigcoefs,j,:)), 
                                   egval, cvar, maxit, tol)
    end
  end
  return dsq
end

""" Computes multitaper eigencoefs """
function multispec_coef(tt, x, u, n, nfft, nfft2)
  x     .-= mean(x)
  (length(x) != length(tt)) && error("The vector of data and the vector of times must have the same lengths.")
  # plan the fft, compute eigencoefficients 
  cent = (vcat(tt, tt[end] .+ collect(1:(nfft-n))) .- nfft/2)/nfft
  p = NFFTPlan(cent, nfft)
  return mapreduce(slep -> nfft_adjoint(p, vcat(slep.*x, zeros(nfft-n)) .+ 0.0im)[(end-nfft2+1):end], hcat, eachcol(u))
end

"""
    mdmultispec(tt, x; <keyword arguments>)

Multitaper power spectrum estimation for time series with missing data (gaps)

...
# Arguments
 - `tt::Vector{T} where T<:Float64`: the vector containing the time indices
 - `x::Vector{T}`: data vector
 - `bw::Float64 = 5/length(tt)`: bandwidth of estimate
 - `k::Int64 = 2*bw*length(x)-1`: number of slepian tapers, must be <=2*bw*length(x) 
 - `dt::T = tt[2]-tt[1]`: sampling rate in time units 
 - `nz::Float64 = 0.0`: zero padding factor
 - `Ftest::Bool = true`: Compute the F-test p-value
 - `jk::Bool = true`: Compute jackknifed confidence intervals
 - `dof::Bool = false`: Compute degrees of freedom for the adaptively weighted 
          spectrum estimate
 - `lambdau::Union{Tuple{Array{Float64,1},Array{Float64,2}},Nothing} = nothing`: Slepians, if precomputed
...

...
# Outputs
 - `pkg::MtSpec` struct containing the spectrum
 - `nu1::Vector{Float64}` optional vector containing the degrees of freedom, given
 if the `dof` kwarg is set to `true`.
...

See also: [`multispec`](@ref), [`mdslepian`](@ref)
"""
function mdmultispec(tt::Union{Vector{Int64},Vector{Float64}}, x::Vector{Float64}; 
                bw=5/length(tt), k=Int64(2*bw*size(x,1)-1), 
                lambdau::Union{Tuple{Array{Float64,1},
                               Array{Float64,2}},Nothing} = nothing,
                dt=tt[2]-tt[1], nz=0, Ftest=true, jk=true,
                dof=false)
  lambda,u = (lambdau == nothing) ? mdslepian(bw, k, tt) : lambdau
  s2    = var(x)
  n, nfft, nfft2 = _pregap(tt, x, nz)
  ak = multispec_coef(tt, x, u, n, nfft, nfft2)
  sxx   = zeros(size(ak,1))
  d = multispec_aweight(ak, lambda, s2)
  Threads.@threads for j in 1:nfft2
        sxx[j] = _dot(d[j,:],abs2.(view(ak,j,:)))
  end
  sxx[2:(nfft2-1)] .*= 2
  d = sxx*sqrt.(lambda')./(sxx*lambda' .+ s2*(1.0 .- lambda'))
  nu1      = 2*d.^2*lambda
  coefswts = Ecoef(ak,d.^2)
  jv       = jk ? jknife(coefswts,nothing,:spec)[2] : nothing
  # F-test
  if Ftest
    dpsw0 = sum(u,dims=1)'
    dpsw02= sum(dpsw0.^2)
    mu    = ak*dpsw0/dpsw02
    num   = (nu1 .- 2).*abs2.(mu)*dpsw02
    denom = 2*sum(abs2.(ak - mu*dpsw0'), dims = 2)
    ft    = vec((num./denom))
    try 
        Fpval = 1.0 .- cdf.(FDist.(2,2*nu1 .- 2), ft)
    catch
        Fpval = nothing
        println("Degrees of freedom get too small to assess F-test p-value.")
    end
  else
    Fpval = nothing
  end
  Tv = nothing
  # Package up the outputs
  pkg = MtSpec((1/dt)*range(0,1,length=nfft)[1:length(sxx)], dt*sxx, nothing, 
              MtParams(bw*n, k, n, tt[2]-tt[1], 2*(nfft2-1), 1, nothing),
              coefswts, Fpval, jv, Tv)
  if dof 
      return pkg, nu1
  else
      return pkg
  end
end

"""
    mdmultispec(tt, x, y; <keyword arguments>)

Multitaper coherence estimation for time series with missing data (gaps)

...
# Arguments
 - `tt::Vector{T} where T<:Float64`: the vector containing the time indices
 - `x::Vector{T}`: data vector 1
 - `y::Vector{T}`: data vector 2
 - `bw::Float64 = 5/length(tt)`: bandwidth of estimate
 - `k::Int64 = 2*bw*length(x)-1`: number of Slepian tapers, must be <=2*bw*length(x) 
 - `dt::T = tt[2]-tt[1]`: sampling rate in time units 
 - `nz::Float64 = 0.0`: zero padding factor
 - `Ftest::Bool = true`: Compute the F-test p-value
 - `jk::Bool = true`: Compute jackknifed confidence intervals
 - `lambdau::Union{Tuple{Array{Float64,1},Array{Float64,2}},Nothing} = nothing`: Slepians, if precomputed
...

...
# Outputs
 - `pkg::MtCoh` struct containing the coherence
...

See also: [`multispec`](@ref), [`mdslepian`](@ref)
"""
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
  (length(x) != length(y)) && error("The two series must have the same lengths.")
  n, nfft, nfft2 = _pregap(t, x, nz)
  lambda,u = (lambdau == nothing) ? mdslepian(bw, k, t) : lambdau
    
  x .-= mean(x)
  y .-= mean(y)
  sx2    = var(x)
  sy2    = var(y)
  sxy   = zeros(nfft2)
    
  cent = (vcat(t, t[end] .+ collect(1:(nfft-n))) .- nfft/2)/nfft
  p   = NFFTPlan(cent, nfft)
  axk = mapreduce(slep -> nfft_adjoint(p, vcat(slep.*x, zeros(nfft-n)) .+ 0.0im)[(end-nfft2+1):end], hcat, eachcol(u))
  ayk = mapreduce(slep -> nfft_adjoint(p, vcat(slep.*y, zeros(nfft-n)) .+ 0.0im)[(end-nfft2+1):end], hcat, eachcol(u))
  outputcoefs = [Ecoef(axk,nothing),Ecoef(ayk,nothing)]
    
  # Jacknife 
  sxy, svar = jknife(outputcoefs..., :coh)
  ph_xy, phvar = jknife_phase(outputcoefs...)
  Tv = nothing

  return MtCoh((1/dt)*range(0,1,length=nfft)[1:nfft2], sxy, ph_xy, 
                MtParams(bw*n, k, n, dt, 2*(nfft2-1), 1, nothing),
                outputcoefs, [svar, phvar], Tv)
end

"""
    mdmultispec(tt, x; <keyword arguments>)

Multitaper coherence estimation for multiple time series with the same missing data (gaps)

...
# Arguments
 - `tt::Vector{T} where T<:Float64`: the vector containing the time indices
 - `x::Matrix{T}`: time series in the columns of a matrix
 - `bw::Float64 = 5/length(tt)`: bandwidth of estimate
 - `k::Int64 = 2*bw*length(x)-1`: number of slepian tapers, must be <=2*bw*length(x) 
 - `dt::T = tt[2]-tt[1]`: sampling rate in time units 
 - `nz::Float64 = 0.0`: zero padding factor
 - `Ftest::Bool = false`: Compute the F-test p-value
 - `jk::Bool = true`: Compute jackknifed confidence intervals
 - `lambdau::Union{Tuple{Array{Float64,1},Array{Float64,2}},Nothing} = nothing`: Slepians, if precomputed
...

...
# Outputs
 - `Tuple{Vector{MtSpec},Matrix{MtCoh},Nothing}` struct containing the spectra, 
coherences, and T^2 test significances (currently set to return nothing)
...

See also: [`multispec`](@ref), [`mdslepian`](@ref)
"""
function mdmultispec(t::Union{Vector{Int64}, Vector{Float64}}, 
                xx::Matrix{Float64};
                bw::Float64 = 5/length(t),
                k::Int64    = Int64(2*bw*size(xx,1) - 1),
                lambdau::Union{Tuple{Array{Float64,1},
                               Array{Float64,2}},Nothing} = nothing,
                dt::Float64 = 1.0,
                nz::Union{Int64,Float64}   = 0, 
                jk::Bool = false, Ftest::Bool = false) 

  n, p = size(xx)
  n, nfft, nfft2 = _pregap(t, xx[:,1], nz)
  fgrid = range(0,1,length=nfft)[1:nfft2]
  if p > 3
    println("You are computing $(Int64((p)*(p+1)/2 - p)) cross-spectra/coherences.")
  end

  # Compute the array of dpss vectors if they aren't given.
  lambdau = (lambdau == nothing) ? mdslepian(bw, k, t) : lambdau
  
  # Get the spectra
  specs     = mapslices(x -> mdmultispec(t, x, bw = bw, k = k, dt = dt, nz = nz, 
                        Ftest = Ftest,
                        lambdau = lambdau, 
                        jk=jk), xx, dims=1)[:]
 
  # Get the coherences
  coherences = Array{MtCoh,2}(undef, p, p)
  for x in CartesianIndex.(filter(x -> x[2]>x[1], 
                Tuple.(eachindex(view(coherences,1:p,1:p)))))
      # Jacknife 
      sxy, svar = jknife(specs[x[2]].coef, specs[x[1]].coef,:coh)
      ph_xy, phvar = jknife_phase(specs[x[2]].coef, specs[x[1]].coef)
      coherences[x] = MtCoh((1/dt)*fgrid, sxy, ph_xy, 
                MtParams(bw*n, k, n, dt, 2*(nfft2-1), 1, nothing),
                nothing, [svar, phvar], nothing)
  end
  Tv = nothing
  return (specs, coherences, Tv)
end

"""
    mdslepian(w, k, t)

Generalized prolate spheroidal sequences for the 1D missing data problem

...
# Arguments
 - `w::Float64`: the bandwidth
 - `k::Int64`: number of Slepian tapers, must be <=2*bw*length(x) 
 - `t::Vector{Int64}`: vector containing the time indices
...

...
# Outputs
 - `lambda,u::Tuple{Vector{Float64}, Vector{Float64}}`: tuple containing the 
 concentrations and the tapers
...

See also: [`mdmultispec`](@ref), [`gpss`](@ref)
"""
function mdslepian(w, k, t)
  n           = length(t)
  a           = 2*w*ones(n,n)
  for i = 1:n
      j       = i .+ 1:n
      a[i,j]  = sin.(2*pi*w*(t[i] .- t[j]))./(pi*(t[i] .- t[j]))
      a[j,i]  = a[i,j]
  end
  lambda,v  = eigs(a, nev = k, which = :LR)
  u         = copy(v)
  for i = 1:2:k
      if mean(real.(u[:,i])) < 0 
        u[:,i] = -u[:,i] 
      end
  end
  for i = 2:2:k-1
      if real(u[2,i] - u[1,i]) < 0
        u[:,i] = -u[:,i] 
      end
  end
  return (lambda, u)
end

"""
    gpss(w, k, t, f; <keyword arguments>)

Generalized prolate spheroidal sequences on an unequal grid

...
# Arguments
 - `w::Float64`: the bandwidth
 - `k::Int64`: number of Slepian tapers, must be <=2*bw*length(x) 
 - `t::Vector{Int64}`: vector containing the time indices
 - `f::Float64`: frequency at which the tapers are to be computed
 - `beta::Float64 = 0.5`: analysis half-bandwidth (similar to Nyquist rate)
...

This function is currently not exported, use `Multitaper.gpss`.

See also: [`mdmultispec`](@ref), [`mdslepian`](@ref)
"""
function gpss(w::Float64, k::Int64, t::Union{Vector{Int64},Vector{Float64}}, 
        f::Float64; beta::Float64 = 0.5)
  n           = length(t)
  a           = 2*w*ones(n,n)
  b = 2.0*beta*ones(n,n) .+ 0.0im
  beta = 0.5
  for i = 1:n
      for j in (i+1):n
          a[i,j]  = sin.(2*pi*w*(t[i] .- t[j]))./(pi*(t[i] .- t[j]))
          b[i,j]  = exp.(-2*pi*1.0im*f*(t[i] .- t[j])).*sin.(2*pi*beta*(t[i] .- t[j]))./(pi*(t[i] .- t[j])) 
      end
  end
  R = cholesky(Hermitian(b))
  V = eigen(Symmetric(real.(inv(R.L)*Symmetric(a)*inv(R.U))))
    lambda = V.values[end:-1:(end-k+1)]
    v = inv(R.U)*V.vectors[:,end:-1:(end-k+1)]
      u         = copy(v)
  for i = 1:2:k
      if mean(real.(u[:,i])) < 0 
        u[:,i] = -u[:,i] 
      end
  end
  for i = 2:2:k-1
      if real(u[2,i] - u[1,i]) < 0
        u[:,i] = -u[:,i] 
      end
  end
  return (lambda, u, R)
end
