
# Utility functions for Multitaper.jl

# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

""" Unwrap the phase """
function unwrapphase(y::Vector{Float64} , typ::Symbol = :rad)
  lims         = (typ==:rad) ? pi : 180
  y[1], y[end] = (0.0, 0.0)
  x            = copy(y)
  dy           = diff(y)
  for i = (length(y)-1):-1:1
    if abs(dy[i]) >= lims
      x[(i+1):end] .+= -sign(dy[i])*(2*lims)
    end
  end
  return x
end

""" Alternative phase unwrapper which takes a complex argument. """
function unwrapphase(y::Vector{ComplexF64},typ::Symbol=:deg)
  if (typ != :deg) && ( typ != :rad) 
    error("Unwrapping must be done in either degrees (:deg) or radians (:rad)")
  end
  # Difference the phase angles
  x = map(i->angle(y[i]/y[i-1]),2:length(y))
  x = vcat(0.0,x)
  return cumsum(x)*((typ == :deg)*180.0/pi + (typ == :rad))
end

""" If freq is a LinRange corresponding to the Fourier frequencies and f_spot is the frequency of
interest, freqfinder gives the index of the Fourier frequency closest to f_spot """
function freqfinder(f_spot::Float64, freq::LinRange{Float64})::Int64
  out  = 1
  mino = abs(f_spot - freq[1])
  minn = abs(f_spot - freq[1])
  for f in freq
    minn = abs(f_spot - f)
    if minn <= mino
      out += 1
      mino = minn
    else
      break
    end
  end
  out -= 1
  return out
end

""" Frequencies given as floats are transformed into indices, as integers, left unchanged when
integers are given"""
freq_to_int(freq::Union{Vector{Float64},Vector{Int64}}, lengt::Int64, dt::Float64=1.0) = 
(typeof(freq) == Vector{Int64}) ? freq : Vector{Int64}(round.(lengt*freq*dt) .+ 1)

""" Integer valued indices are converted to floats (frequencies) and left unchanged when floats are given"""
int_to_freq(int::Union{Vector{Int64},Vector{Float64}}, lengt::Int64, dt::Float64=1.0) =
(typeof(int) == Vector{Float64}) ? int : Vector{Float64}((int .- 1)/(lengt*dt))

""" Expected Jackknife variance """
function EJN(DoF::Union{Int64,Float64})::Float64
  Kft = Float64(DoF)/2
  return (Kft-1)^3*(Kft-3)/((Kft-0.5)*Kft*(Kft-2)^3)
end

""" Return W when NW is given """
function BW(NW::Float64, N::Int64, dt::Float64=1.0)::Float64
  return (NW/N)*(1.0/dt)
end

""" Logrange - gives values equally spaced on a log10 scale
e.g. LogRange(1,5,5) = [0.9, 0.99, ..., 0.99999]"""
function LogRange(logmin, logmax, N::Int64)
  return 1.0 .- 10.0 .^ (-collect(LinRange(logmin,logmax,N)))
end

""" CDF for jackknifed msc when true coherence is 0. 
Thomson & Chave 1991 Book chapter, integrating Ezq (2.53) when γ^2 = 0 """
function mscsig(csq::Float64, ntf::Union{Float64, Int64})
  if (csq > 1.0)||(csq < 0.0)
    error("Squared coherence must be between 0 and 1.")
  end
  return (1.0 - (1.0 - csq)^(ntf - 1.0))
end

""" Inverse significance for jackknifed msc when true coherence is 0.
Thomson & Chave 1991 Book chapter"""
function invmscsig(sig::Float64, ntf::Union{Float64, Int64})
  if (sig >= 1.0) || (sig <= 0.0)
    error("Significance must be between 0, and 1, noninclusive")
  end
  return 1.0 - (1.0 - sig)^(1/(ntf - 1.0))
end

# Transformed coherences

""" Magnitude squared coherence variance stabilizing transform"""
function atanhtrans(c::Union{Float64,ComplexF64}, ntf::Union{Float64,Int64})
  y = try atanh(abs(c)^2); catch;  atanh(1.0-1e-15); end
  return sqrt(2*ntf-2)*y
end

""" Sort-of-inverse of the MSC variance-stabilizing transform"""
function tanhtrans(trcsq::Union{Float64,ComplexF64}, ntf::Union{Float64,Int64}) 
  return tanh(trcsq/sqrt(2*ntf-2))
end

# For example, if the transformed coherence is 7.3 and 6 dof, the significance is
# invmscsig(tanhtrans(7.3,6),6)
# if the significance is 0.9 and 6 dof the transformed msc is
# atanhtrans(sqrt(mscsig(0.9,6)),6)


