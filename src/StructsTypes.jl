
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

#### MT-like structs 

""" 
Multitaper parameters: 
  time bandwidth (NW) as Float, 
  number of tapers (K), 
  number of samples (N), 
  sampling rate (dt) in temporal units (e.g. seconds), 
  padded length (M), 
  number of segments (nsegments), and 
  overlap if nsegments is greater than 1, nothing otherwise.
"""
struct mtparams
  NW      ::Float64
  K       ::Int64
  N       ::Int64
  dt      ::Float64
  M       ::Int64
  nsegments ::Int64
  overlap ::Union{Nothing, Float64}
end

""" 
The ecoef structure holds 
  multitaper eigencoefficients (coef) and, optionally, 
  adaptive weights (wts) 
"""
mutable struct ecoef
  coef    ::Union{Matrix{ComplexF64}}
  wts     ::Union{Matrix{ComplexF64}, Matrix{Float64}, Nothing}
end

""" 
The multitaper spectrum is given as a mtspec structure which holds 
  frequency (f), 
  spectrum (S), 
  phase (optional), 
  chosen values of the multitaper time bandwidth product etc of type mtparams (params)
  eigencoefficients (coef, optional), 
  Ftest values (Fpval, optional), 
  jackknife output (jkvar, optional), and
  Tsquared test results (Tsq_pval, optional). 
"""
mutable struct mtspec
  f       ::Union{Array{Float64, 1}, LinRange{Float64}}
  S       ::Union{Array{Float64, 1}, Array{ComplexF64, 1}}
  phase   ::Union{Array{Float64, 1}, Nothing}
  params  ::mtparams
  coef    ::Union{ecoef, Array{ecoef,1}, Array{Float64,2}, Nothing}
  Fpval   ::Union{Array{Float64, 1}, Nothing}
  jkvar   ::Union{Array{Float64, 1}, Array{Float64, 2}, Array{Array{Float64, 1}, 1}, Nothing}
  Tsq_pval::Union{Float64, Array{Float64, 1}, Nothing}
end

""" 
The multitaper autocovariance function is given in the mtacvf structure, which holds
  lags (lags), 
  autocovariance function (acvf), 
  mtparams (params)
"""
mutable struct mtacvf
  lags    ::Union{Array{Float64, 1}, LinRange, UnitRange}
  acvf    ::Union{Array{Float64, 1}, Array{ComplexF64, 1}} 
  params  ::mtparams
end

""" 
The multitaper autocorrelation function is given in the mtacf structure which holds 
  lags (lags), 
  autocorrelation function (acf), 
  mtparams (params)
"""
mutable struct mtacf
  lags    ::Union{Array{Float64, 1}, LinRange, UnitRange}
  acf     ::Union{Array{Float64, 1}, Array{ComplexF64, 1}} 
  params  ::mtparams
end

""" 
The multitaper cepstrum is given in the mtceps structure which holds 
  lags (lags), 
  cepstrum (ceps), 
  mtparams (params)
"""
mutable struct mtceps
  lags    ::Union{Array{Float64, 1}, LinRange, UnitRange}
  ceps    ::Union{Array{Float64, 1}, Array{ComplexF64, 1}} 
  params  ::mtparams
end

""" 
The multitaper coherence structure, mtcoh, holds 
  frequency (f), 
  coherence (coh), 
  phase (phase), 
  eigencoefficients (ecoeffs, optional), 
  jackknife output (jkvar, optional), and 
  Tsquared test results (Tsq, optional).
"""
mutable struct mtcoh
  f       ::Union{Array{Float64, 1}, LinRange}
  coh     ::Array{Float64, 1}
  phase   ::Union{Array{Float64, 1}, Nothing}
  params  ::mtparams
  coef    ::Union{ecoef, Array{ecoef,1}, Array{Float64,2}, Array{Array{Float64, 2}, 1}, Nothing}
  jkvar   ::Union{Array{Float64, 1}, Array{Float64, 2}, Array{Array{Float64, 1}, 1}, Nothing}
  Tsq_pval::Union{Float64, Array{Float64, 1}, Nothing}
end

""" 
The multitaper transfer function is given as a mttransf structure which holds 
          frequency (f), 
          transfer function (transf), 
          phase (phase), 
          mtparams (params),
          eigencoefficients (ecoeffs, optional), 
          jackknife output (jkvar, optional), and 
          Tsquared test results (Tsq, optional).
"""
mutable struct mttransf
  f       ::Union{Array{Float64, 1}, LinRange}
  transf  ::Union{Array{Float64, 1}, Array{Float64, 2}}
  phase   ::Union{Array{Float64, 1}, Nothing}
  params  ::mtparams
  coef    ::Union{ecoef, Array{ecoef,1}, Array{Float64,2}, Array{Array{Float64, 2}, 1}, Nothing}
  jkvar   ::Union{Array{Float64, 1}, Array{Float64, 2}, Array{Array{Float64, 1}, 1}, Nothing}
end

""" 
The multitaper cross-covariance function is given in the mtccvf structure which holds 
  lags (lags), 
  cross-covariance function (acvf), 
  mtparams (params)
"""
mutable struct mtccvf
  lags    ::Union{Array{Float64, 1}, LinRange, UnitRange}
  ccvf    ::Union{Array{Float64, 1}, Array{ComplexF64, 1}} 
  params  ::mtparams
end

""" 
The multitaper cross-correlation function is given as a mtccf structure which holds 
  lags (lags), 
  cross-correlation function (acf), 
  mtparams (params)
"""
mutable struct mtccf
  lags    ::Union{Array{Float64, 1}, LinRange, UnitRange}
  ccf     ::Union{Array{Float64, 1}, Array{ComplexF64, 1}} 
  params  ::mtparams
end

"""
Multitaper complex demodulates are held in the Demodulate struct which contains
          magnitude (Mag), 
          phase (Phase)
"""
mutable struct Demodulate
  time    ::LinRange{Float64}
  mag     ::Array{Float64, 1}
  phase   ::Array{Float64, 1}
end
