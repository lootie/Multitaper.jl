"""
Goals: -multithread the Bronez method for faster computations
       -save the tapers to examine after the Bronez spectrum is computed (use JLD2)
"""

using Multitaper, FINUFFT, LinearAlgebra, Distributions, Statistics, StatsBase, Polynomials
using Base.Threads: nthreads, @spawn, @threads


"""
    gpss_uneven(w, k, t, f)
Generalized, orthogonalized prolate spheroidal sequences on an
unevenly spaced grid
...
# Positional Arguments
 - `w::Float64`: the bandwidth
 - `k::Int64`: number of generalized prolate spheriodal sequence tapers, must be <=2*bw*length(x) 
 - `t::Union{Vector{Int64}, Vector{Float64}}`: vector containing the time indices
 - `f::Float64`: frequency at which the tapers are to be computed
...
...
# Outputs
 - `lambda::Vector{Float64}` the concentrations of the generalized prolate spheroidal sequences
 - `u::Matrix{Float64}` the matrix containing the sequences themselves
...
See also: [`dpss_tapers`](@ref), [`mdmultispec`](@ref), [`mdmultispec_noadapt`](@ref), [`mdslepian`](@ref)
"""
function gpss_uneven(w::Float64, k::Int64, t::Union{Vector{Int64},Vector{Float64}}, f::Float64)
  beta        = 0.5
  n           = length(t)
  a           = 2*w*ones(n,n)
  b = 2.0*beta*ones(n,n) .+ 0.0im
  for i = 1:n
      for j in (i+1):n
          a[i,j]  = sin.(2*pi*w*(t[i] .- t[j]))./(pi*(t[i] .- t[j]))
          b[i,j]  = exp.(-2*pi*1.0im*f*(t[i] .- t[j])).*
                      sin.(2*pi*beta*(t[i] .- t[j]))./(pi*(t[i] .- t[j])) 
      end
  end
  # If the Cholesky factorization fails, add a small number to the diagonal and
  # try again. Then remove that same value from the eigenvalues. 
  R, fudge = try (cholesky(Hermitian(b)), 0.0)
  catch
    (cholesky(Hermitian(b) + Matrix(I, size(b)...)*(1e-10)), 1.0)
  end
  V = eigen(Symmetric(real.(inv(R.L)*Symmetric(a)*inv(R.U))))
  lambda = V.values[end:-1:(end-k+1)] .- fudge*(1e-10)
  u = V.vectors[:,end:-1:(end-k+1)]
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


"""
    rescale_time(times, observations, tscale)
Sort and rescale vector of observation times to get characteristic dt ~ 1
...
# Positional Arguments
 - `times::Vector{T} where T<:Real`: vector of observing times
 - `observations::Vector{P} where P<:Number`: vector of data
 - `tscale::Float64`: characteristic timestep for rescaling
...
...
# Outputs
 - `t::Vector{T} where T<:Real`: sorted, rescaled time stamps
 - `obs::Vector{P} where P<:Number`: sorted data
See also: [`gpss_uneven`](@ref), [`remove_repeats`](@ref), [`bronez`](@ref)
"""
function rescale_time(times::Vector{T}, observations::Vector{P},
                      tscale::Float64) where{T<:Real, P<:Number}
    
    t = copy(times)
    obs = copy(observations)
    
    # Sort so time is monotonically increasing
    tind = sortperm(t)
    t = t[tind]
    t = t .- t[1]
    obs = obs[tind]
    
    # Divide timestamps by characteristic timescale
    t ./= tscale

    # Return sorted and altered timestamps, sorted observed data 
    return t, obs

end


"""
    remove_repeats(t::Union{Vector{Float64}, Vector{Int64}}, x::Vector{Float64})
If there are duplicate time stamps, remove the duplicate by
averaging the two measurements
...
# Positional Arguments
 - `t:: Union{Vector{Float64}, Vector{Int64}}`: vector of time stamps
 - `x::Vector{Float64}`: vector of data
...
Changes to t and x vectors are made inplace
...
...
See also: [`gpss_uneven`](@ref), [`rescale_time`](@ref), [`bronez`](@ref)
"""
function remove_repeats(t::Union{Vector{Float64}, Vector{Int64}}, x::Vector{Float64})
# t = timestamps; x = measurements
    
    if length(t) != length(x)
        println("Vectors of timestamps and measurements must have the same length. Exiting...")
        return
    end

    cm = countmap(t)
    if (maximum(values(cm)) > 1)
        println("Repeats found")
    end
    
    for (k, v) in cm
        if (v > 1)
            repeat_inds = findall(t .== k)
            x[repeat_inds[1]] = mean(x[repeat_inds])
            map(ind -> deleteat!(t, ind), repeat_inds[2:end])
            map(ind -> deleteat!(x, ind), repeat_inds[2:end])
        end 
    end  # All Julia loops must have end statements
    
end 


# Least-squares fit a sinusoid to a tapered time series
function lsq_sinusoid(tobs::Vector{T}, tapered_series::Vector{P}, fr::Float64) where{T<:Number,P<:Number}
    design_sinusoid = [cos.(2*pi*fr .*tobs) sin.(2*pi*fr .*tobs)]
    amplitudes = design_sinusoid \ tapered_series
    return amplitudes[1] - 1im * amplitudes[2]
end


"""
    bspec(time, obs, W, K; <keyword arguments>)
Bronez (1988) multitaper power spectrum estimate
...
# Arguments
## Positional Arguments
 - `time::Union{Vector{Float64}, Vector{Int64}}`: vector of observing times
 - `obs::Vector{Float64}`: vector of observations
 - `W::Float64`: bandwidth
 - `K::Int64`: number of generalized prolate spheroidal sequence tapers, must be <=2*W*length(obs) 
## Keyword Arguments
 - `oversample::Float64 = 1.0`: oversampling factor for frequency grid
 - `Ftest::Bool = false`: whether to perform F-test for harmonic components
 - `return_tapers::Bool = false`: whether to return array containing tapers
...
...
# Outputs
 - `spectrum::MTSpectrum`: MTSpectrum structure containing the spectrum
 - `tapers::Array{Float64}`: optional array containing tapers for each frequency
...
See also: [`gpss_uneven`](@ref), [`mdmultispec`](@ref), [`multispec`](@ref)
"""
function bspec(time::Union{Vector{Float64}, Vector{Int64}},
               obs::Vector{Float64}, W::Float64, K::Int64;
               oversample::Float64 = 1.0, Ftest::Bool = false,
               return_tapers::Bool = false)
    
    if length(time) != length(obs)
        throw(ArgumentError("Time and data vectors must be the same length."))
    end
    
    println("Number of threads: $(nthreads())")
    
    t = copy(time)
    x = copy(obs)

    # Average data with repeated timestamps and sort
    remove_repeats(t, x)
    
    # Number of surviving data points
    N = length(t)
    
    # Get frequency grid
    beta = 0.5 # process bandwidth
    M = Int64(round(N * oversample))
    if mod(M,2) != 0 
        M = M + 1
    end
    M2 = Int(M/2)
    freq = collect(range(-beta, beta, length = M + 1)[1:M])

    # Save multitaper parameters to MTParameters structure
    params = MTParameters(N * W, K, N, 1.0, M, 1, nothing)
    
    # Time shift so t[1] = 0, subtract linear trend
    t = t .- t[1]
    lintrend = Polynomial(Polynomials.fit(t, x, 1))
    x_det = x .- lintrend.(t)
    
    # Make a 3-d array that holds the tapers
    #  (row: time, column: Slepian order, layer: frequency)
    tapers = zeros(Float64, N, K, M2)
        
    # Fill the taper array
    @threads for m = 1:M2
        tapers[1:N, 1:K, m] = gpss_uneven(W, K, t, freq[M2+m])[2]
    end
    println("Tapers computed")
    
    # Function to apply taper for a specific frequency
    taper(m) = mapslices(gpss_k -> gpss_k .* x_det, tapers[1:N, 1:K, m], dims=1)
    
    # Function to calculate eigencoefficient for a specific frequency
    eigenco(m) = mapslices(tapered_k -> lsq_sinusoid(t, tapered_k, freq[M2+m]), taper(m), dims=1)
    
    # Multithreaded calculation of all eigencoefficients
    chunk_size = max(1, M2 ÷ nthreads())
    tasks = map(Iterators.partition(1:M2, chunk_size)) do chunk
        @spawn mapreduce(m -> eigenco(m), vcat, chunk)
    end
    eco = EigenCoefficient(mapreduce(fetch, vcat, tasks), nothing)
    println("Eigencoefficients computed")
    
    # Normalize eigencoefficients
    df = freq[2] - freq[1]
    x_var = var(x_det)
    f_domain_avg = sum(vec(mean(abs2.(eco.coef), dims=2)) .* df) / beta
    eco.coef = eco.coef .* sqrt.(x_var ./ f_domain_avg)
    
    # Jackknife (not so computationally intensive b/c no gpss computation or fitting)
    jknifed = Multitaper.jknife(eco, nothing, :spec)
    println("Jackknifing done")
    
    # F-test with multithreaded spectral window zero-frequency amplitude calculation
    if Ftest
        tasks = map(Iterators.partition(1:M2, chunk_size)) do chunk
            @spawn mapreduce(m -> sum(tapers[1:N, 1:K, m], dims=1), vcat, chunk)
        end
        gpsw0  = mapreduce(fetch, vcat, tasks)
        println("Window amplitude at zero frequency computed")
        gpsw0sq= sum(abs2, gpsw0, dims = 2)
        mu    = sum(broadcast(/, eco.coef .* gpsw0, gpsw0sq), dims = 2)  
        num   = real.((K - 1) * abs2.(mu) .* gpsw0sq)
        denom = real.(sum(abs2, (eco.coef .- broadcast(*, mu, gpsw0)), dims = 2))
        Ft    = vec((num ./ denom))
        Fpval = 1.0 .- cdf.(FDist.(2, 2 * K .- 2), Ft)
    else
        Fpval = nothing
    end
    
    # Package Bronez spectrum into MTSpectrum structure
    spectrum =  MTSpectrum(freq[(M2+1):M], vec(mean(abs2.(eco.coef), dims=2)), 
                           nothing, params, eco, Fpval, jknifed[2], nothing)
    
    # Return Bronez spectrum and (optionally) tapers
    if return_tapers
        return (spectrum, tapers)
    else
        return spectrum
    end
       
end


"""
    rescale_frequency(spectrum, tsc)
If observation times were rescaled, convert frequencies to physical units
...
# Positional Arguments
 - `spectrum::MTSpectrum`: MTSpectrum structure returned by bspec()
 - `tsc::Float64`: characteristic timestep
...
Spectrum is modified inplace
...
...
See also: [`rescale_time`](@ref), [`bspec`](@ref), [`remove_repeats`](@ref)
"""
function rescale_frequency(spectrum::MTSpectrum, tsc::Float64)
    spectrum.f = spectrum.f ./ tsc
    spectrum.S = spectrum.S ./ tsc
end
