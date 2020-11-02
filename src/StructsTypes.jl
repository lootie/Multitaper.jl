
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

#### MT-like structs 


const OutputPoints = Union{Vector{Float64}, StepRangeLen{Float64}, UnitRange{Int64}}

""" 
Multitaper parameters `MTParameters` struct contains
  - time bandwidth (NW) as Float, 
  - number of tapers (K), 
  - number of samples (N), 
  - sampling rate (dt) in temporal units (e.g. seconds), 
  - padded length (M), 
  - number of segments (nsegments), and 
  - overlap if nsegments is greater than 1, nothing otherwise.
"""
struct MTParameters
  NW      ::Float64
  K       ::Int64
  N       ::Int64
  dt      ::Float64
  M       ::Int64
  nsegments ::Int64
  overlap ::Union{Nothing, Float64} # nothing? fixable?
end

""" 
The `EigenCoefficient` structure holds 
  - multitaper eigencoefficients (coef) and, optionally, 
  - adaptive weights (wts) 
"""
mutable struct EigenCoefficient
  coef    ::Matrix{ComplexF64}
  wts     ::Union{Matrix{Float64}, Nothing}
end

""" 
The multitaper spectrum is given as a `MTSpectrum` structure which holds 
  - frequency (f), 
  - spectrum (S), 
  - phase (optional), 
  - chosen values of the multitaper time bandwidth product etc of type MTParameters
    (params)
  - eigencoefficients (coef, optional), 
  - Ftest values (Fpval, optional), 
  - jackknife output (jkvar, optional), and
  - Tsquared test results (Tsq_pval, optional). 
"""
mutable struct MTSpectrum{C,J,P}
  f       ::OutputPoints
  S       ::Vector{Float64} 
  phase   ::Union{Array{Float64, 1}, Nothing}
  params  ::MTParameters
  coef    ::C
  Fpval   ::Union{Array{Float64, 1}, Nothing}
  jkvar   ::J
  Tsq_pval::P
end

""" 
The multitaper autocovariance function is given in the `MTAutocovarianceFunction` structure, which holds
  - lags (lags), 
  - autocovariance function (acvf), 
  - MTParameters (params)
"""
mutable struct MTAutocovarianceFunction{T}
  lags    ::OutputPoints
  acvf    ::Vector{T}
  params  ::MTParameters
end

""" 
The multitaper autocorrelation function is given in the `MTAutocorrelationFunction`
structure which holds 
  - lags (lags), 
  - autocorrelation function (acf), 
  - MTParameters (params)
"""
mutable struct MTAutocorrelationFunction{T}
  lags    ::OutputPoints
  acf     ::Vector{T}
  params  ::MTParameters
end

""" 
The multitaper cepstrum is given in the `MTCepstrum` structure which holds 
  - lags (lags), 
  - cepstrum (ceps), 
  - MTParameters (params)
"""
mutable struct MTCepstrum{T}
  lags    ::OutputPoints
  ceps    ::Vector{T}
  params  ::MTParameters
end

""" 
The multitaper coherence structure, `MTCoherence`, holds 
  - frequency (f), 
  - coherence (coh), 
  - phase (phase), 
  - eigencoefficients (coef, optional), 
  - jackknife output (jkvar, optional), and 
  - Tsquared test results (Tsq, optional).
"""
mutable struct MTCoherence{C,J,P}
  f       ::OutputPoints
  coh     ::Array{Float64, 1}
  phase   ::Union{Array{Float64, 1}, Nothing}
  params  ::MTParameters
  coef    ::C
  jkvar   ::J
  Tsq_pval::P
end

""" 
The multitaper transfer function is given as a `MTTransferFunction` structure which holds 
  - frequency (f), 
  - transfer function (transf), 
  - phase (phase), 
  - MTParameters (params),
  - eigencoefficients (EigenCoefficientfs, optional), 
  - jackknife output (jkvar, optional), and 
  - Tsquared test results (Tsq, optional).
"""
mutable struct MTTransferFunction{C,J}
  f       ::OutputPoints
  transf  ::Union{Array{Float64, 1}, Array{Float64, 2}}
  phase   ::Union{Array{Float64, 1}, Nothing}
  params  ::MTParameters
  coef    ::C
  jkvar   ::J
end

""" 
The multitaper cross-covariance function is given in the `MtCrossCovarianceFunction` structure which holds 
  - lags (lags), 
  - cross-covariance function (acvf), 
  - MTParameters (params)
"""
mutable struct MtCrossCovarianceFunction{T}
  lags    ::OutputPoints
  ccvf    ::Vector{T}
  params  ::MTParameters
end

""" 
The multitaper cross-correlation function is given as a `MTCrossCorrelationFunction` structure which holds 
  - lags (lags), 
  - cross-correlation function (acf), 
  - MTParameters (params)
"""
mutable struct MTCrossCorrelationFunction{T}
  lags    ::OutputPoints
  ccf     ::Vector{T} 
  params  ::MTParameters
end

"""
Multitaper complex demodulates are held in the `Demodulate` struct which contains
  - magnitude (Mag), 
  - phase (Phase)
"""
mutable struct Demodulate
  time    ::OutputPoints
  mag     ::Array{Float64, 1}
  phase   ::Array{Float64, 1}
end
