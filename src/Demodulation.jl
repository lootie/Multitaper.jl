
""" Complex demodulation: converted by CLH from the R code by Rahim in 2019. 
Inputs: vector x to be demodulated, dt sampling rate (s), 
f0 the frequency to center at, 
nw the time bandpValidth product, 
blockLen, the length of the dpss filter 
wrapphase indicates whether or not to unwrap the phase"""
function demodulate(x::Vector{Float64}, dt::Float64, f0::Float64, NW::Float64, 
                      blockLen::Int64, wrapphase::Bool = true)
  # Generate the zeroth order Slepian taper as a lowpass filter
  dpVal      = dpss_tapers(blockLen, NW, 1, :tap)
  # Complex factor to shift by
  cshift  = exp.(-1.0im * 2.0 * pi * f0 * dt * collect(0:(blockLen-1)))
  cshift  = broadcast(*, broadcast(*, cshift, dpVal), (2.0/sum(dpVal)))
  nResVal = length(x) - blockLen + 1
  # Apply the fiter by convolution
  cdm     = mapreduce(x-> x'*cshift, vcat, x[i:(i+blockLen-1)] for i in 1:nResVal)[:,1]
  # Process the phase information
  phase   = angle.(cdm) * 180/pi
  phase   = wrapphase ? unwrapphase(phase,:deg) : phase
  phase   = phase - 360 * dt * f0 * (1:nResVal)
  return Demodulate(LinRange(0, length(x), nResVal)*dt, abs.(cdm), phase)
end


