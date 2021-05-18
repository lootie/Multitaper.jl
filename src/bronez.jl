
# Here is a manual nfft:
function slow_nfft_adj(t::Union{Vector{Int64}, Vector{Float64}}, 
                x::Union{Vector{ComplexF64}, Vector{Float64}, Matrix{Float64}}, beta::Float64;
                nfft::Union{Int64} = length(t), nfft2::Union{Int64} = Int64(floor(length(t)/2)+1))
  Nyq_ratio = 2*beta # Scale the quantity (i-1)/nfft by the ratio of the process bandwidth to the
                     # ideal Nyquist frequency: Nyq_ratio = beta/0.5
  return map(i -> exp.(0.0 .- 2.0*pi*im*(Nyq_ratio*(i - 1)/nfft)*t)'*x, 1:nfft2)   
end

""" 
    repeated_times(time)

Verify that the time vector is monotone increasing and has no repeated values

...
# Arguments
 - `times::Vector{T} where T<: Number`: the vector containing the times

# Outputs
  - `(good_times, bad_times)::Tuple{Vector{Float64}}`: The tuple of good and bad indices

...

"""
function repeated_times(time::Vector{T}) where T<:Number
    dt = diff(time)
    (sum(dt .< 0.0) > 0) && error("The time vector is not monotone increasing.")
    good_times = vcat(findall(dt .> 1e-15), length(time))
    bad_times  = findall(diff(time) .< 1e-15)
    (length(bad_times) != 0) && println("Warning: There are repeated time stamps. 
                                       Averaging data points with same time index.")
    return (good_times, bad_times)
end

""" 
    ave_repeats(dat, times)

Average data from the repeated indices 

...
# Arguments
  - `dat::Vector{P} where P <: Number`: The vector containing the data
  - `times::Vector{T}` where T <: Number`: The vector containing the times

# Outputs
  - `(new_dat, times)::Tuple{Vector{Float64}`: The tuple containing the new data and the new times
"""
function ave_repeats(dat::Union{Vector, Matrix}, times::Vector)
    good_times, bad_times = repeated_times(time)
    new_dat = (typeof(dat) <: Matrix) ? dat[good_times,:] : dat[good_times] 
    for i = bad_times[end:-1:1]
        temp_idx = findall(times[i] .== times )
        if typeof(dat) <: Vector
            new_dat[i] = mean(dat[temp_idx])
        else typeof(dat) <: Matrix
            new_dat[i,:] = mean(dat[temp_idx,:], dims=1)
        end
    end
    return (new_dat, times[good_times])
end          


"""
    bspec(times, dat, W, K, beta, nz, Ftest)

Computes the Bronez spectrum of an unequally-spaced time series

...
# Positional Arguments
 - `times::Vector{T} where T<: Number`: the vector containing the times
 - `dat::Vector{T} where T<:Number`: the vector containing the time series 
 - `W::Float64`: bandwidth of estimate
 - `K::Int64`: number of slepian tapers, must be <= 2*NW
 - `beta::Float64`: estimated process bandwidth (aka Nyquist rate)
 - `nz::Union{Float64,Int64} = length(t)`: Number of frequencies
 - `Ftest::Bool = false`: Whether to compute the F-test or not
...

...
# Outputs
 - MTSpectrum struct containing the Bronez spectrum
...

...
# Example usage

```<julia>
N = 256
t = collect(1:N).^(1.05)
W = 0.008
K = 5
x = randn(N)
bet = 0.5 / (last(t) / (N-1))
M = 2*N
S = bspec(t, x, W, K, bet, nz, true)
```
...

See also: [`multispec`](@ref), [`mdmultispec`](@ref), [`mdslepian`](@ref), [`gpss`](@ref)
"""
function bspec(times::Vector{T}, dat::Vector{P}, W::Float64, K::Int64, beta::Float64, 
               nz::Float64 = length(times), Ftest::Bool = false) where{T<:Number,P<:Number}
    x, t = ave_repeats(dat, times)
    N, M, M2 = _pregap(t, x, nz)
    freq = collect(range(-1.0, 1.0, length = M + 1) * beta)
    params = MTParameters(N * W, K, N, 1.0, M, 1, nothing)
    eigenc(j, fr, x) = mapslices(slep -> slow_nfft_adj(t, slep .* x, beta, nfft = M, 
                               nfft2 = M2)[j], gpss(W, K, t, fr, beta = beta)[2], 
                               dims = 1)
    eco = EigenCoefficient(mapreduce(j -> eigenc(j, freq[j + Int(M / 2)], x), vcat, 
                                                 1:(Int(M / 2))), nothing)   
    jknifed = jknife(eco,nothing,:spec) 
    if Ftest
        freq = range(-beta, beta, length = M + 1)
        gpsw0  = mapreduce(fr -> sum(gpss(W, K, t, fr, beta = beta)[2],
                                     dims = 1),vcat, freq[1:Int(M / 2)])
        gpsw0sq= sum(abs2, gpsw0, dims = 2)
        mu    = sum(broadcast(/, eco.coef .* gpsw0, gpsw0sq), dims = 2)  
        num   = real.((K - 1) * abs2.(mu) .* gpsw0sq)
        denom = real.(sum(abs2, (eco.coef .- broadcast(*, mu, gpsw0)), dims = 2))
        Ft    = vec((num ./ denom))
        Fpval = 1.0 .- cdf.(FDist.(2, 2 * K .- 2), Ft)
    else
        Fpval = nothing
    end
    return MTSpectrum(freq[(Int(M/2)+1):M], mean(abs2.(eco.coef), dims=2)[:], 
                      nothing, params, eco, Fpval, jknifed[2], nothing)
end

"""
    bspec(time, dat1, dat2, W, K, bet, nz; <keyword arguments>)

Computes the Bronez coherence or cross-spectrum of two unequally-spaced time series

...
# Positional Arguments
 - `time::Vector{T} where T<:Number`: the vector containing the times 
 - `dat1::Union{Vector{P}, EigenCoefficient} where P<:Number`: the vector containing the first time series 
 - `dat2::Union{Vector{P}, EigenCoefficient} where P<:Number`: the vector containing the second time series
 - `W::Float64`: bandwidth of estimate
 - `K::Int64`: number of slepian tapers, must be <= 2*NW
 - `bet::Float64`: Nyquist frequency
 - `nz::Union{Float64,Int64} = length(t)`: Number of frequencies

# Keyword Arguments
 - `outp::Symb = :coh`: Output, either `:cross` for cross spectrum or `:coh` (default)
 - `params::Union{MTParameters,Nothing} = nothing`: parameters struct, important when x,y are EigenCoefficients
 - `Ftest::Bool = false`: Whether to compute the F-test or not
...

...
# Outputs
 - MTCoherence containing the Bronez coherence
...

...
# Example usage

```<julia>
N = 256
t = collect(1:N).^(1.05)
W = 0.008
K = 5
x = randn(N)
y = randn(N) # Incoherent
M = 2*N
beta = 0.5
S = bspec(t, x, y, W, K, M, beta)
```
...

See also: [`multispec`](@ref), [`mdmultispec`](@ref), [`mdslepian`](@ref), [`gpss`](@ref)
"""
function bspec(time::Vector{T}, dat1::Union{Vector{P},EigenCoefficient}, 
        dat2::Union{Vector{P},EigenCoefficient}, W, K, bet, nz = 0.0; 
        outp = :coh, params::Union{MTParameters,Nothing} = nothing,
        Ftest = false) where{T<:Number,P<:Number}
    
    if typeof(dat1) != EigenCoefficient
        x, t = ave_repeats(hcat(dat1,dat2), time)
        y = x[:,2]
        x = x[:,1]
        N, M, M2 = _pregap(t, x, nz)
        freq = range(-0.5, 0.5, length = M + 1)
        params = MTParameters(N * W, K, N, 1.0, M, 1, nothing)
        eigenc(j, fr, x) = mapslices(slep -> slow_nfft_adj(t, slep .* x, nfft = M, nfft2 = M2)[j], 
            gpss(W, K, t, fr, beta = bet)[2], dims=1)
    
        eco_x = EigenCoefficient(mapreduce(j -> eigenc(j, freq[j + Int(M / 2)], x), 
                                           vcat, 1:(Int(M / 2))), nothing)
        eco_y = EigenCoefficient(mapreduce(j -> eigenc(j, freq[j + Int(M / 2)], y), 
                                           vcat, 1:(Int(M / 2))), nothing)
    else
        M = params.M
        freq = range(-0.5, 0.5, length = params.M + 1)
        eco_x = dat1
        eco_y = dat2
    end
    if outp == :coh
        jknifed = jknife(eco_x, eco_y, :coh) 
        jphase = jknife(eco_x, eco_y, :phase) 
        return MTCoherence(freq[(Int(M / 2) + 1):M], jknifed[1], jphase[1], params, 
                           [eco_x, eco_y], [jknifed[2], jphase[2]], nothing) 
    else 
        # returns cross spectrum
        jknifed = jknife(eco_x, eco_y, :spec)
        jphase = jknife(eco_x, eco_y, :phase) 
        return MTSpectrum(freq[(Int(M / 2) + 1):M], 
                          ((eco_x.coef) * conj(eco_y.coef)' / K)[:], nothing, 
                          params, nothing, nothing, jknifed[2], nothing) 
    end
end

"""
    bspec(time, dat, W, K, bet, nz; <keyword arguments>)

Computes the Bronez spectra and coherences of p unequally-spaced time series

...
# Positional Arguments
 - `time::Vector{T} where T<:Number`: the vector containing the times 
 - `dat::Matrix{Vector{P}} where T<:Number`: the matrix containing the time series in its columns
 - `W::Float64`: bandwidth of estimate
 - `K::Int64`: number of slepian tapers, must be <= 2*NW
 - `bet::Float64`: Nyquist frequency
 - `nz::Union{Float64,Int64} = length(t)`: Number of frequencies

# Keyword Arguments
 - `outp::Symb = :coh`: Output, either `:cross` for cross spectrum or `:coh` (default)
 - `Ftest::Bool = false`: Whether to compute the F-test or not
...

...
# Outputs
 - tuple(Vector{MTSpectrum},Matrix{MTCoherence},Nothing) containing the Bronez spectra and coherences
...

...
# Example usage

```<julia>
N = 256
t = collect(1:N).^(1.05)
W = 0.008
K = 5
x = randn(N)
y = randn(N) # Incoherent
M = 2*N
beta = 0.5
S = bspec(t, hcat(x, y), W, K, beta, M, outp = :coh, Ftest = true)
```
...

See also: [`multispec`](@ref), [`mdmultispec`](@ref), [`mdslepian`](@ref), [`gpss`](@ref)
"""
function bspec(t::Vector{T}, x::Matrix{P}, W, K, bet, nz = 0.0; 
        outp = :coh, Ftest = false) where{T<:Number,P<:Number}
    N, p = size(x)
    freq = 1.0 * range(-0.5, 0.5, length = M + 1)
    # Get the spectra
    specs   = map(y -> bspec(t, y, W, K, nz, Ftest), x[:, k] for k in 1:p)
    params  = specs[1].params
    # Get the coherences or cross spectra
    crosspecs = (outp == :cross) ? Array{MTSpectrum, 2}(undef, p, p) : 
                                  Array{MTCoherence, 2}(undef, p, p)
    for i in CartesianIndex.(filter(dat -> dat[2] > dat[1], 
                             Tuple.(eachindex(view(crosspecs, 1:p, 1:p)))))
      crosspecs[i] = bspec(t, specs[i[1]].coef, specs[i[2]].coef, W, K, nz,
                           outp = ((outp == :cross) ? :spec : :coh), 
                           params = params)
    end
    return (specs, crosspecs, nothing)
end
