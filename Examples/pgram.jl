
# Sometimes it can be helpful to compute a periodogram estimate of the msc
# This code works a lot like R, so for info about spans and the number of segments check there. 
# It ought to be well known that in order for a periodogram estimate of the msc to not be degenerately
# unity, one needs to either compute the cross spectrum over multiple segments or use smoothing.
# This code implements both. 

####### Periodogram cross-spectrum (Welsh)

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
  lengt, fftleng, halffreq = output_len(Fx1[1:s],pad)
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
  lengt, fftleng, halffreq = output_len(Fx1, pad)
  Pxx = abs2.(fft(vcat(Fx1,zeros(fftleng - length(Fx1)))))[1:halffreq]/length(Fx1)
  return pgram((1/dt)*LinRange(0,0.5,halffreq), dt*Pxx, 0, 1, 0.0)
end

""" Fejer kernel 
Equation (198b) in Percival and Walden """
function Fejer(N::Int64; dt::Float64 = 1.0)
  freq = LinRange(0,1.0,N+1)[1:N]
  return (sin(N*pi*dt*f)).^2 ./ (N*(sin(pi*f*dt)).^2)
end

""" Guts function for computing a single smoothed cross spectrum """
function smoothedcspec_guts(Fx1::Vector{Float64},
                            Fx2::Vector{Float64};
                            segments::Int64 = 1, 
                            spans::Union{Int64,Vector{Int64}} = 1, pad::Union{Float64,Int64} = 1.0,
                            taper::Union{Symbol, Nothing, Tuple{Symbol,Nothing}, Tuple{Symbol,Float64},
                                         Tuple{Symbol,Vector{Float64}}} = nothing,
                            opt::Symbol = :coh, overlap::Float64 = 0.0
                           )
  # Welsh technique
  t,s = blockerr(length(Fx1),segments,overlap=overlap)
  lengt, fftleng, halffreq = output_len(Fx1[1:s], pad)
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
  u2 = (taper == (:default, 0.1)) ? (1 - (5/8) * 0.1 * 2) : 1.0
  tap ./= sqrt(sum(abs2,tap)) 
  T = zeros(fftleng,length(t))*1.0im
  for i = 1:length(t)
    T[:,i] = fft(vcat(tap.*Fx1[t[i]:(t[i]+s-1)], zeros(fftleng -  s)))
    T[:,i] .*= conj.(fft(vcat(tap.*Fx2[t[i]:(t[i]+s-1)], zeros(fftleng -  s))))
    for j = 1:length(spans)
      T[:,i] = simple_smooth(T[:,i],spans[j])
    end
  end
  S = mean(T,dims=2)[1:halffreq,1]
  S[1] = 0.5*(S[2] + S[end])
  if (opt == :coh)||(opt == :cross)
    return (abs2.(S)/u2, unwrapphase(vcat(0.0,angle.(S[2:(end-1)]),0.0)*180/pi),:deg)
  elseif (opt == :cqspec)
    return (real.(S)/u2, imag.(S)/u2)
  else
    error("Output must be :coh, :cross, or :cqspec")
  end
end

""" Takes multiple data columns and computes the msc's, spectra and phases """
function smoothedmsc(Fx::Matrix{Float64}; spans::Union{Int64,Vector{Int64}} = 1, pad::Union{Float64,Int64} = 1.0, 
                     taper::Union{Symbol, Nothing, Tuple{Symbol,Nothing}, Tuple{Symbol,Float64}, 
                                  Tuple{Symbol,Vector{Float64}}} = nothing, overlap::Float64 = 0.0,
                     segments::Int64 = 1, opt::Symbol = :coh, dt::Float64 = 1.0)
  N,M = size(Fx);
  # Compute the spectra first, discarding the phase
  spectra = Vector{direct}(undef,M)
  for i = 1:M 
    spectra[i] = sm_pgram(Fx[:,i], segments = segments, spans = spans, pad = pad, dt = dt, taper = taper)
  end
  # Nothing fancy for now
  smoother = (:daniell, 0.1)
  taper = (taper == nothing) ? (:default, 0.1) : taper
  u2 = ((taper == (:default, 0.1))||nothing) ? (1 - (5/8) * 0.1 * 2) : 1.0
  # Compute the cross-spectra
  if opt == :coh
    S = Matrix{pmsc}(undef, M, M)
    for sep = 1:M
      for i = 1:(M-sep)
        # Normalizing step
        msc,phase = smoothedcspec_guts(demean(Fx[:,i]), 
                                       demean(Fx[:,i+sep]), pad = pad, segments = segments, spans = spans, 
                                       opt = opt, taper = taper, overlap = overlap)
        msc ./= (u2*(1.0/dt)^2*spectra[i].Sxx.*spectra[i+sep].Sxx)
        S[i,i+sep] = pmsc(LinRange(0,0.5,length(msc))*(1/dt), msc, phase, taper, smoother, spans, segments, 0.0)
      end
    end
  elseif opt == :cross
    S = Matrix{pcross}(undef, M, M)
    for sep = 1:M
      for i = 1:(M-sep)
        # Normalizing step
        Pxy,phase = smoothedcspec_guts(demean(Fx[:,i]), 
                                       demean(Fx[:,i+sep]), pad = pad, segments = segments, spans = spans, 
                                       opt = opt, taper = taper,overlap = overlap)
        S[i,i+sep] = pcross(LinRange(0,0.5,length(Pxy))*(1/dt), dt*Pxy, phase, taper, smoother, 
                            spans, segments, overlap)
      end
    end
  elseif opt == :cqspec
    S = Matrix{pcspecqspec}(undef, M, M)
    for sep = 1:M
      for i = 1:(M-sep)
        cospec, quadspec = smoothedcspec_guts(demean(Fx[:,i]), 
                                              demean(Fx[:,i+sep]),pad = pad, segments = segments, spans = spans, 
                                              taper = taper, opt = opt,overlap = overlap)
        S[i,i+sep] = pcspecqspec(LinRange(0,0.5,length(cospec))*(1/dt), dt*cospec, dt*quadspec, 
                                 taper, smoother, 
                                 spans, segments, overlap)
      end
    end
  else
    error("Output must be one of :coh, :cqspec, :cross")
  end
  return spectra, S
end

""" Implements a Capon spectral estimator. Adapted from matlab code by R Moses.
Input: y (data vector), m (integer length of the Capon filter), L (number of estimated spectral samples)
Output: phi, R, a """
function capon_spec(y::Vector{Float64}, m::Int64, L::Int64)
  #  
  # Form the sample covariance matrix
  N = length(y)
  R = sum(map(i -> y[i:-1:(i-m)]*transpose(y[i:-1:(i-m)]), (m+1):N))
  R ./= (N-m)
  #
  # Invert
  IR = inv(R)
  #
  # Compute the spectrum 
  ak = map(k->exp.(-2.0im*pi*(k-1)/L*collect(0:m)), 1:L)
  phi = map(a->real.(a'*IR*a),ak)
  #
  return phi
end


