
# MT routines by Chave, adapted for Julia from Matlab by C Haley in 2020
# Original Chave license text appended at the end of this doc. 
# Credit is due to the original author, cite: 
#
# Chave, Alan D. "A multitaper spectral estimator for time-series with missing data." Geophysical
# Journal International 218.3 (2019): 2165-2178.
#

using Statistics, FFTW, Arpack, LinearAlgebra, Distributions

""" Multi-taper power spectrum estimation with f-test and reshaping
for time series with missing data
input arguments 
   tt -- real vector of time (required)
   xx -- real vector of data (required)
   bw -- bandwidth of estimate, 5/length(t) default
   k -- number of slepian tapers, must be <=2*bw*length(x), 2*bw*length(x)-1 default
   nz -- zero padding factor, 0 default
   alpha -- probability level for reshaping, 1 if none, 1 default
output arguments 
   sxx -- power spectrum vector of length length(x)/2+1 (required)
   nu1 -- degrees-of-freedom vector for sxx of length length(x)/2+1
   ft -- f-test vector of length length(x)/2+1
   il -- frequency indices of spectral lines
   plin -- spectral line power vector
   srr -- reshaped power spectrum vector of length length(x)/2+1 if alpha<1
   nu2 -- degrees-of-freedom vector for srr of length length(x)/2+1
"""
function MDmwps(tt::Union{Vector{Int64}, Vector{Float64}}, 
                xx::Union{Vector{Float64}, Matrix{Float64}};
                bw::Float64 = 5/length(tt),
                k::Int64    = Int64(2*bw*size(xx,1) - 1),
                nz::Union{Float64,Int64}   = 0, 
                alpha::Float64 = 1.0)
  n, p  = (typeof(tt) == Matrix{Float64}) ? size(tt) : (length(tt),1)
  t     = copy(tt)
  x     = copy(xx)
  nfft  = length(x)
  if mod(nfft,2) != 0 
    nfft = nfft + 1
  end
  nfft  = (nz + 1)*nfft
  nfft2 = Int64(round(nfft/2)) + 1
  s2    = var(x)
  e     = Array{Matrix{ComplexF64},1}(undef, nfft2)
  sxx   = zeros(nfft2)
  ak    = 0.0im*zeros(nfft2,k)
  lambda,u = MDslepian(bw, k, tt)
  d = zeros(nfft2,k)
  for i = 1:nfft2
      f      = (i - 1)/nfft
      e[i]   = transpose(u.*repeat(exp.(0.0 .- 2.0*pi*im*f*t),outer=(1,k)))       
      # eigencoefficients
      ak[i,:] = e[i]*x
      # eigenspectra at frequency index i
      sk     = abs.(ak[i,:]).^2
      # Adaptive weighting 
      (old, new) = vcat(ones(2).*0.5, zeros(k-2)), zeros(k) 
      sxx[i] = aweighted(old, new, k, sk, lambda, s2, 15,0.05)
      d[i,:] = new
  end
  sxx[2:nfft2-1] = 2*sxx[2:nfft2-1]
  nu1   = 2*d*lambda
  # F-test
  dpsw0 = real.(sum(e[1],dims=2))
  mu    = ak*dpsw0/sum(dpsw0.^2)
  num   = (nu1 .- 2).*abs2.(mu).*sum(dpsw0.^2)
  denom = 2*sum(abs2.(ak[1:nfft2,:] - kron(mu,transpose(dpsw0))), dims = 2)
  ft    = (num./denom)[:]
  Falpha= quantile(FDist(2,2*k-2),alpha)
  il1   = findall(ft .>= Falpha)
  # reshape the spectrum after removing lines
  dpsw  = mapreduce(x -> sum(x, dims=2)', vcat, e)
  dpsw  = fftshift(vcat(dpsw, dpsw[nfft2-1:-1:2,:]),1)
  zk    = fftshift(vcat(ak, conj(ak[nfft2-1:-1:2,:])),1)
  il    = []
  plin  = []
  if length(il1) != 0
      n   = Int64((nz + 1)*round(bw*length(x)) + 1)
      i   = 2:length(il1)
      i1  = vcat(findall(il1[i] .!= il1[i.-1] .+ 1), [length(il1)])
      j   = -n .+ 1:n .- 1
     #
     ii = [1]
      for i = 1:length(i1)
        if (i1[i] == ii)
          push!(il,il1[ii[end]]) 
          push!(ii,ii+1)
        else 
          push!(il,findall(ft .== maximum(ft[il1[ii[end]]:il1[i1[i]]]))[1])
          push!(ii,i1[i] + 1)
        end
      end
      for i2 in il
        m         = nfft2 .+ i2 .+ j .- 1
        mm        = m .> size(zk,1)
        m[mm]     = m[mm] .- size(zk,1)
        zk[m,1:k] -= mu[i2]*dpsw[nfft2 .+ j,1:k]
        push!(plin, mean(sum(abs2.(mu[i2]*dpsw[nfft2.+j,1:k]))))
      end
  end
  zk  = ifftshift(zk,1)
  s2  = mean(abs2.(zk[1,:] .+ 2*sum(abs2.(zk[2:nfft2-1,:])) .+ abs2.(zk[nfft2,:])))/nfft
  sk  = abs2.(zk)
  srr = zeros(nfft2)
  for i = 1:nfft2
      # Adaptive weighting again
      (old, new) = vcat(ones(2).*0.5, zeros(k-2)), zeros(k) 
      srr[i] = aweighted(old, new, k, sk[i,:], lambda, s2, 15, 0.05)
      d[i,:] = new
  end
  srr[2:nfft2-1]  .*= 2.0
  if length(il1) != 0
      if il[1] != 1
          plin  .*= 2.0
      else
          plin[2:length(plin)] .*= 2.0
      end
  end
  nu2 = 2*d*lambda
  return (sxx, nu1, ft, il, plin, srr, nu2)
end

""" Helper function MDslepian computes generalized slepian function for 1D missing data problem
  input variables
       w = analysis half bandwidth
       k = number of eigenvalue/vectors to compute
       t = time vector
  output variables
       lambda = eigenvalues
       u = eigenvectors
"""
function MDslepian(w::Float64, k::Int64, t::Union{Vector{Int64},Vector{Float64}})
  # Random.seed!(2147483647)
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

#= 
Copyright (c) 2019, Alan Chave
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution
* Neither the name of Woods Hole Oceanographic Institution nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=#
