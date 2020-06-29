
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

#### MT-like structs 


const OutputPoints = Union{Vector{Float64}, StepRangeLen{Float64}, UnitRange{Int64}}

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
struct MtParams
  NW      ::Float64
  K       ::Int64
  N       ::Int64
  dt      ::Float64
  M       ::Int64
  nsegments ::Int64
  overlap ::Union{Nothing, Float64} # nothing? fixable?
end

""" 
The Ecoef structure holds 
  multitaper eigencoefficients (coef) and, optionally, 
  adaptive weights (wts) 
"""
mutable struct Ecoef
  coef    ::Matrix{ComplexF64}
  wts     ::Union{Matrix{Float64}, Nothing}
end

""" 
The multitaper spectrum is given as a MtSpec structure which holds 
  frequency (f), 
  spectrum (S), 
  phase (optional), 
  chosen values of the multitaper time bandwidth product etc of type MtParams
    (params)
  eigencoefficients (coef, optional), 
  Ftest values (Fpval, optional), 
  jackknife output (jkvar, optional), and
  Tsquared test results (Tsq_pval, optional). 
"""
mutable struct MtSpec{C,J,P}
  f       ::OutputPoints
  S       ::Vector{Float64} 
  phase   ::Union{Array{Float64, 1}, Nothing}
  params  ::MtParams
  coef    ::C
  Fpval   ::Union{Array{Float64, 1}, Nothing}
  jkvar   ::J
  Tsq_pval::P
end

""" 
The multitaper autocovariance function is given in the MtAcvf structure, which holds
  lags (lags), 
  autocovariance function (acvf), 
  MtParams (params)
"""
mutable struct MtAcvf{T}
  lags    ::OutputPoints
  acvf    ::Vector{T}
  params  ::MtParams
end

""" 
The multitaper autocorrelation function is given in the MtAcf structure which holds 
  lags (lags), 
  autocorrelation function (acf), 
  MtParams (params)
"""
mutable struct MtAcf{T}
  lags    ::OutputPoints
  acf     ::Vector{T}
  params  ::MtParams
end

""" 
The multitaper cepstrum is given in the MtCeps structure which holds 
  lags (lags), 
  cepstrum (ceps), 
  MtParams (params)
"""
mutable struct MtCeps{T}
  lags    ::OutputPoints
  ceps    ::Vector{T}
  params  ::MtParams
end

""" 
The multitaper coherence structure, MtCoh, holds 
  frequency (f), 
  coherence (coh), 
  phase (phase), 
  eigencoefficients (coef, optional), 
  jackknife output (jkvar, optional), and 
  Tsquared test results (Tsq, optional).
"""
mutable struct MtCoh{C,J,P}
  f       ::OutputPoints
  coh     ::Array{Float64, 1}
  phase   ::Union{Array{Float64, 1}, Nothing}
  params  ::MtParams
  coef    ::C
  jkvar   ::J
  Tsq_pval::P
end

""" 
The multitaper transfer function is given as a MtTransf structure which holds 
          frequency (f), 
          transfer function (transf), 
          phase (phase), 
          MtParams (params),
          eigencoefficients (Ecoeffs, optional), 
          jackknife output (jkvar, optional), and 
          Tsquared test results (Tsq, optional).
"""
mutable struct MtTransf{C,J}
  f       ::OutputPoints
  transf  ::Union{Array{Float64, 1}, Array{Float64, 2}}
  phase   ::Union{Array{Float64, 1}, Nothing}
  params  ::MtParams
  coef    ::C
  jkvar   ::J
end

""" 
The multitaper cross-covariance function is given in the MtCcvf structure which holds 
  lags (lags), 
  cross-covariance function (acvf), 
  MtParams (params)
"""
mutable struct MtCcvf{T}
  lags    ::OutputPoints
  ccvf    ::Vector{T}
  params  ::MtParams
end

""" 
The multitaper cross-correlation function is given as a MtCcf structure which holds 
  lags (lags), 
  cross-correlation function (acf), 
  MtParams (params)
"""
mutable struct MtCcf{T}
  lags    ::OutputPoints
  ccf     ::Vector{T} 
  params  ::MtParams
end

"""
Multitaper complex demodulates are held in the Demodulate struct which contains
          magnitude (Mag), 
          phase (Phase)
"""
mutable struct Demodulate
  time    ::OutputPoints
  mag     ::Array{Float64, 1}
  phase   ::Array{Float64, 1}
end
