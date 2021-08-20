
# Sometimes it can be helpful to compute a periodogram estimate of the msc
# This code works a lot like R, so for info about spans and the number of segments check there. 
# It ought to be well known that in order for a periodogram estimate of the msc to not be degenerately
# unity, one needs to either compute the cross spectrum over multiple segments or use smoothing.
# This code implements both. 

using FFTW, RecipesBase, Statistics, DSP

# Periodogram-like struct
"""
Periodogram estimator contained in struct pgram contains
          frequency (freq)
          periodogram estimate (Pxx)
          spans parameter for Daniell window (spans)
          number of segments (segments)
          overlap of the segments (overlap)
""" 
struct pgram
  freq    ::Union{Vector{Float64},LinRange{Float64}}
  Pxx     ::Vector{Float64}
  spans   ::Union{Int64,Vector{Int64}}
  segments::Int64
  overlap ::Float64
end
    
# Simple periodogram recipes

@recipe function pgramplot(Pxx::pgram)
  yscale -->  :log10
  xguide --> "Frequency (Hz)"
  yguide --> "Periodogram, spans = $(Pxx.spans), segments = $(Pxx.segments)"
  label --> "Periodogram"
  @series begin
    Pxx.freq[2:end], Pxx.Pxx[2:end]
  end
end

"""
Direct estimate of the spectrum (single taper estimate) contains
          frequency (freq)
          spectrum estimate (Sxx)
          taper used (taper)
          smoother used (smoother)
          number of segments (segments)
          overlap of the segments (overlap)
"""
struct direct
  freq    ::Union{Vector{Float64}, LinRange{Float64}}
  Sxx     ::Vector{Float64}
  taper   ::Union{Nothing, Symbol, Tuple{Symbol, Nothing}, Tuple{Symbol, Vector{Float64}}, Tuple{Symbol, Float64}}
  smoother::Union{Nothing, Symbol, Tuple{Symbol, Nothing}, Tuple{Symbol, Vector{Float64}}, Tuple{Symbol, Float64}, Tuple{Symbol, Vector{Int64}}, Tuple{Symbol, Int64}}
  segments::Int64
  overlap ::Float64
end

@recipe function directplot(Pxx::direct)
  yscale -->  :log10
  xguide --> "Frequency (Hz)"
  if Pxx.segments == 1
    ystr = "Direct estimate, computed using a\n"*taperstring(Pxx.taper)
    ystr *= (Pxx.smoother == (:daniell,1)) ? "." : "\n and smoothed with a " * smootherstring(Pxx.smoother)
  elseif Pxx.segments >= 1
    ystr = "Welch estimate on $(Pxx.segments) segments with\n $(Pxx.overlap*100.0) percent overlap,"
    ystr *= " computed using a\n"*taperstring(Pxx.taper)
    ystr *= (Pxx.smoother == (:daniell,1)) ? "." : "\n and smoothed with a " * smootherstring(Pxx.smoother)
  end
  label --> ystr
  yguide --> "Spectrum Estimate"
  @series begin
    Pxx.freq[2:end], Pxx.Sxx[2:end]
  end
end

""" Helper function to parse the smoother tuple """
function smootherstring(smoother::Union{Symbol, Tuple{Symbol, Vector{Float64}}, Tuple{Symbol, Float64}, Tuple{Symbol, Vector{Int64}}, Tuple{Symbol, Int64}})
  # right now the default is a Daniell smoother with some spans parameter.
  if smoother == (:daniell,1)
    return "no smoothing"
  else
    return "modified\n Daniell smoother with spans: $(smoother[2])"
  end
end

####### Periodogram cross-spectrum (Welch)

""" Demean-er """
function demean(x::Vector{Float64})
  return x .- mean(x)
end

""" Simple smoother.
This one makes sure that a one-sided spectrum gets smoothed with its edges matched circularly. """
function simple_smooth(x::Union{Vector{ComplexF64},Vector{Float64}}, spans::Int64)
  if !isodd(spans)
    error("spans should be odd.")
  elseif spans == 1
    return x
  else
    N = Int64(floor(spans/2))
    # Daniell window is 1/(m-1) except at endpoints where it is 1/(2*(m-1))
    c = (typeof(x) == Vector{ComplexF64}) ? complex.(ones(spans))*(1/(spans-1)) : ones(spans)*(1/(spans-1))
    c[[1 end]] ./= 2
    # The c's sum to 1
    # c ./= sum(abs2.(c))
    # do a periodic continuation of the spectrum before median smoothing
    return conv(vcat(x,x,x),c)[(N+1+length(x)):(end-N-length(x))]
  end
end

""" Smoothed periodogram 
Inputs: Data vector Fx1, 
        spans smooths the periodogram (repeatedly) with a Daniell smoother with parameter spans
        dt sampling rate
        segments number of segments
        overlap - proportion by which to overlap the segments
        taper - the (name (symbol), parameters) of the tapering function
Outputs: a plottable pgram or direct estimate """
function sm_pgram(Fx1::Vector{Float64}; spans::Union{Int64,Vector{Int64}} = 1, pad::Union{Float64,Int64} = 1.0,
                  dt::Float64 = 1.0, segments::Int64 = 1, overlap::Float64 = 0.0, ctr::Bool = true,
                  taper::Union{Symbol, Nothing, Tuple{Symbol, Nothing}, Tuple{Symbol, Float64}, 
                               Tuple{Symbol, Vector{Float64}}} = nothing)
  Fx1 = ctr ? demean(Fx1) : Fx1 
  t,s = blockerr(length(Fx1),segments, overlap = overlap)
  lengt, fftleng, halffreq = Multitaper.output_len(Fx1[1:s],pad)
  u2 = (taper == nothing) ? (1 - (5/8) * 0.1 * 2) : 1.0
  if typeof(taper) == Nothing
    tap = getwindow(s,(:default, 0.1))
    taper = (:default, 0.1)
  elseif typeof(taper) == Symbol
    tap = getwindow(s,taper)
    taper = (taper, nothing)
  elseif ((typeof(taper) == Tuple{Symbol, Vector{Float64}})||(typeof(taper) == Tuple{Symbol,Float64}))
    tap = getwindow(s,taper)
  else
    error("Taper can be left blank in order to get the modified Daniell windows")
  end
  tap ./= sqrt(sum(abs2,tap)) 
  pgram_seq = zeros(fftleng,length(t))
  for i = 1:length(t)
    # Equation (206c) for direct estimates
    pgram_seq[:,i] = abs2.(fft(vcat(tap.*Fx1[t[i]:(t[i]+s-1)], zeros(fftleng -  s))))
    # Repeated application of a modified Daniell window as in Bloomfield.
    for j = 1:length(spans)
      pgram_seq[:,i] = simple_smooth(pgram_seq[:,i],spans[j])
    end
  end 
  Pxx = mean(pgram_seq,dims=2)[1:halffreq,1]
  Pxx[1] = 0.5*(Pxx[2] + Pxx[end])
  freq = (1/dt)*LinRange(0,0.5,halffreq)
  return direct(freq, dt*Pxx/u2, taper, (:daniell, spans), segments, overlap)
end

""" Basic periodogram """
function basic_pgram(Fx1::Vector{Float64}; dt::Float64 = 1.0, pad::Union{Float64, Int64} = 1.0, ctr::Bool = true)
  Fx1 = ctr ? demean(Fx1) : Fx1
  lengt, fftleng, halffreq = Multitaper.output_len(Fx1, pad)
  Pxx = abs2.(fft(vcat(Fx1,zeros(fftleng - length(Fx1)))))[1:halffreq]/length(Fx1)
  return pgram((1/dt)*LinRange(0,0.5,halffreq), dt*Pxx, 0, 1, 0.0)
end
