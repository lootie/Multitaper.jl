
# Utility functions for Multitaper.jl

# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

""" 
    unwrapphase(y, typ)

Unwrap the phase 

... 
# Arguments

 - `y::Vector{Float64}`: The vector of phases

 - `typ::Symbol`: Whether to compute the unwrapped phase in degrees or radians
...

...
# Outputs

 - `x::Vector{Float64}`: The vector of unwrapped phases
...
"""
function unwrapphase(y, typ=:rad)
  lims         = (typ==:rad) ? pi : 180
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
function unwrapphase(y::Vector{Complex{T}},typ=:deg) where{T}
  if (typ != :deg) && ( typ != :rad) 
    error("Unwrapping must be done in either degrees (:deg) or radians (:rad)")
  end
  # Difference the phase angles
  x = map(i->angle(y[i]/y[i-1]),1:length(y))
  return cumsum(x)*((typ == :deg)*180.0/pi + (typ == :rad))
end

""" If freq is a StepRangeLen corresponding to the Fourier frequencies and f_spot is the
frequency of interest, freqfinder gives the index of the Fourier frequency closest to
f_spot """
function freqfinder(f_spot, freq)::Int64 
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

""" Frequencies given as floats are transformed into indices, as integers, left
unchanged when integers are given"""
function freq_to_int(freq, lengt, dt=1.0) 
  flag = (typeof(freq) == Vector{Int64}) 
  flag ? freq : Int64.(round.(lengt*freq*dt) .+ 1)
end

""" Integer valued indices are converted to floats (frequencies) and left unchanged
when floats are given""" 
function int_to_freq(int, lengt, dt=1.0) 
  (typeof(int) == Vector{Float64}) ? int : collect((int .- 1)/(lengt*dt))
end

"""
    ejn(DoF)

Expected jackknife variance of an univariate multitaper spectrum estimate with DoF tapers

...
# Arguments
 - `DoF::Int64`: the number of tapers use to compute the spectrum
...

...
# Outputs
 - The output is a `Float64` indicating the expected jackknife variance.
...

See also: [`multispec`](@ref)
"""
function ejn(DoF)
  Kft = Float64(DoF)/2
  return (Kft-1)^3*(Kft-3)/((Kft-0.5)*Kft*(Kft-2)^3)
end

""" Return W when NW is given """
function bw(NW, N, dt=1.0)
  return (NW/N)*(1.0/dt)
end

""" Logrange - gives values equally spaced on a log10 scale
e.g. logRange(1,5,5) = [0.9, 0.99, ..., 0.99999]"""
function logRange(logmin, logmax, N::Int64)
  return 1.0 .- 10.0 .^ (-collect(range(logmin,logmax,length=N)))
end

""" CDF for jackknifed msc when true coherence is 0. 
Thomson & Chave 1991 Book chapter, integrating Ezq (2.53) when γ^2 = 0 """
function mscsig(csq, ntf)
  if (csq > 1.0)||(csq < 0.0)
    error("Squared coherence must be between 0 and 1.")
  end
  return (1.0 - (1.0 - csq)^(ntf - 1.0))
end

""" Inverse significance for jackknifed msc when true coherence is 0.
Thomson & Chave 1991 Book chapter"""
function invmscsig(sig, ntf)
  if (sig >= 1.0) || (sig <= 0.0)
    error("Significance must be between 0, and 1, noninclusive")
  end
  return 1.0 - (1.0 - sig)^(1/(ntf - 1.0))
end

# Transformed coherences

"""
    atanhtrans(c,ntf)

Magnitude squared coherence variance stabilizing transform

...
# Arguments
 - `c::Float64`: the value of the coherence
 - `ntf::Int64`: the number of tapers use to compute the coherence
...

...
# Outputs
 - The output is a `Float64` indicating the transformed MSC
...

See also: [`multispec`](@ref)
"""
function atanhtrans(c, ntf)
  y = try atanh(abs(c)^2); catch;  atanh(1.0-1e-15); end
  return sqrt(2*ntf-2)*y
end

"""
    tanhtrans(trcsq,ntf)

Inverse of the magnitude squared coherence variance stabilizing transform

...
# Arguments
 - `trcsq::Float64`: The transformed squared coherence
 - `ntf::Int64`: the number of tapers use to compute the coherence
...

...
# Outputs
 - The output is a `Float64` indicating the reverse-transformed MSC
...

# Example

If the transformed coherence is 7.3 and 6 dof, the significance is

```julia-repl
julia> invmscsig(tanhtrans(7.3,6),6)
```

If the significance is 0.9 and 6 dof the transformed msc is
```julia-repl    
julia> atanhtrans(sqrt(mscsig(0.9,6)),6)
```

See also: [`multispec`](@ref)
"""
function tanhtrans(trcsq, ntf) 
  return tanh(trcsq/sqrt(2*ntf-2))
end



function get_plan(n)
  work1 = zeros(ComplexF64, n)
  work2 = zeros(ComplexF64, n)
  plan  = plan_fft!(work1)
  (work1, work2, plan)
end

function conv(x,y)
  @assert length(x) == length(y) "This convolution only supports equal length"
  len = nextprod([2,3,5,7], 2*length(x)-1)
  (work1, work2, plan) = get_plan(len)
  fill!(work1, zero(ComplexF64))
  fill!(work2, zero(ComplexF64))
  work1[1:length(x)] .= x
  work2[1:length(y)] .= y
  plan*work1
  plan*work2
  work1.*=work2
  plan\work1
  isreal(x) ? real(work1) : work1
end


