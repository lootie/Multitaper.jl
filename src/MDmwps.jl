
# MT routines by Chave, adapted for Julia from Matlab by C Haley in 2020 Original
# Chave license text appended at the end of this doc.  Credit is due to the original
# author: 
#
# Chave, Alan D. "A multitaper spectral estimator for time-series with missing data."
# Geophysical Journal International 218.3 (2019): 2165-2178.
#
# NB: (i) This code removes the sample mean from the data by default, and
# additionally scales the spectrum by dt, when provided.  (ii) Chave used a
# zero-finding routine to compute adaptive weights while here we use the iterative
# method from Thomson, 82. 
#

using Statistics, FFTW, Arpack, LinearAlgebra, Distributions

""" Multi-taper power spectrum estimation with f-test and reshaping
for time series with missing data
t arguments 
   tt -- real vector of time (required)
   x -- real vector of data (required)
   bw -- bandwidth of estimate, kwarg 5/length(t) default
   k -- number of slepian tapers, must be <=2*bw*length(x), kwarg 2*bw*length(x)-1 
        default
   dt -- sampling frequency in time units (kwarg, default tt[2]-tt[1])
   nz -- zero padding factor, kwarg 0 default
   Ftest -- whether to compute the F-test (kwarg, default true)
   alpha -- probability level for reshaping, 1 if none, 1 default
   jk -- whether to compute jackknifed confidence intervals (kwarg, default true)
   Tsq -- lines at which to compute a T^2 test for multiple line components, 
          returns a list of p-values. 
   dof -- whether to return the degrees of freedom for the adaptively weighted 
          spectrum estimate
output arguments 
   output MtSpec struct containing the spectrum and other results.
   nu -- the optional degrees-of-freedom for the adaptively weighted spectrum 
        estimate
"""
function mdmwps(tt, x; bw=5/length(tt), k=Int64(2*bw*size(x,1)-1), 
                dt=tt[2]-tt[1], nz=0, Ftest=true, alpha=1.0, jk=true,
                Tsq=nothing, dof=false)
  x     .-= mean(x)
  n     = length(tt)
  nfft  = length(x)
  (n != nfft) && error("Time vector and data vector must have the same lengths.")
  # Assuming that equal sampling occurs between the first two observations
  if mod(nfft,2) != 0 
    nfft = nfft + 1
  end
  nfft  = Int64(round((nz + 1)*nfft))
  nfft2 = Int64(round(nfft/2)) + 1
  s2    = var(x)
  e     = Array{Matrix{ComplexF64},1}(undef, nfft2)
  sxx   = zeros(nfft2)
  d     = zeros(nfft2,k)
  ak    = 0.0im*zeros(nfft2,k)
  lambda,u = mdslepian(bw, k, tt)
  for i = 1:nfft2
      f      = (i - 1)/nfft
      e[i]   = transpose(u.*repeat(exp.(0.0 .- 2.0*pi*im*f*tt),outer=(1,k)))       
      # eigencoefficients
      ak[i,:] = e[i]*x
      # eigenspectra at frequency index i
      sk     = abs.(ak[i,:]).^2
      # Adaptive weighting 
      (old, new) = vcat(ones(2).*0.5, zeros(k-2)), zeros(k) 
      d[i,:] = aweighted(old, new, k, sk, lambda, s2, 15, 0.05)
      sxx[i] = _dot(d[i,:],sk)
  end
  sxx[2:nfft2-1] .*= 2
  d = sxx*sqrt.(lambda')./(sxx*lambda' .+ s2*(1.0 .- lambda'))
  nu1      = 2*d.^2*lambda
  coefswts = Ecoef(ak,d.^2)
  jv       = jk ? jknife(coefswts,nothing,:spec)[2] : nothing
  # F-test
  if Ftest
    dpsw0 = real.(sum(e[1],dims=2))
    dpsw02= sum(dpsw0.^2)
    mu    = ak*dpsw0/dpsw02
    num   = (nu1 .- 2).*abs2.(mu)*dpsw02
    denom = 2*sum(abs2.(ak - mu*dpsw0'), dims = 2)
    ft    = vec((num./denom))
    Fpval = 1.0 .- cdf.(FDist.(2,2*nu1 .- 2), ft)
  else
    Fpval = nothing
  end
  # T^2 test
  if typeof(Tsq) != Nothing
    Tsq      = (typeof(Tsq) <: Vector{Number}) ? [Tsq] : Tsq
    map!(x -> freq_to_int(Tsq[x], n, dt), Tsq, eachindex(Tsq))
    Tsq = Vector{Vector{Int64}}(Tsq)  
    if (2*k < (true ? 1 : 2)*maximum(length.(Tsq)))
      error("There are too few tapers for the number of Tsq tests.")
    end
    dcs = vec(map(isodd,1:K).*sum(u,dims=1))
    Tv = map(x->testTsq(dcs,Ecoef(coefswts.coef[Tsq[x],:],nothing)),
              eachindex(Tsq)) 
  else
    Tv = nothing
  end
  # Package up the outputs
  pkg = MtSpec((1/dt)*range(0,1,length=nfft)[1:nfft2], dt*sxx, nothing, 
              MtParams(bw*n, k, n, tt[2]-tt[1], 2*(nfft2-1), 1, nothing),
              coefswts, Fpval, jv, Tv)
  if dof 
      return pkg, nu1
  else
      return pkg
  end
end

""" Helper function mdslepian computes generalized slepian function for 1D missing data problem
  input variables
       w = analysis half bandwidth
       k = number of eigenvalue/vectors to compute
       t = time vector
  output variables
       lambda = eigenvalues
       u = eigenvectors
"""
function mdslepian(w, k, t)
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
