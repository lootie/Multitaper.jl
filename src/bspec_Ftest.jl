
# Here is a manual nfft:
function slow_nfft_adj(t::Union{Vector{Int64}, Vector{Float64}}, 
                x::Union{Vector{ComplexF64}, Vector{Float64}, Matrix{Float64}}, beta::Float64;
                nfft::Union{Int64} = length(t), nfft2::Union{Int64} = Int64(floor(length(t)/2)+1))
  Nyq_ratio = 2*beta # Scale the quantity (i-1)/nfft by the ratio of the process bandwidth to the
                     # ideal Nyquist frequency: Nyq_ratio = beta/0.5
  return map(i -> exp.(0.0 .- 2.0*pi*im*(Nyq_ratio*(i - 1)/nfft)*t)'*x, 1:nfft2)   
end

"""
    bspec(times, x, W, K, <keyword args>)

Computes the Bronez spectrum of an unequally-spaced time series

...
# Positional Arguments
 - `x::Vector{T} where T<:Number`: the vector containing the time series 
 - `W::Float64`: bandwidth of estimate
 - `K::Int64`: number of slepian tapers, must be <= 2*NW
 - `beta::Float64`: estimated process bandwidth (aka Nyquist rate)
 - `nz::Union{Float64,Int64} = length(t)`: Number of frequencies
 - `Ftest::Bool = true`: Whether to compute the F-test or not
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
function bspec(t::Vector{T}, x::Vector{P}, W::Float64, K::Int64, beta::Float64, 
               nz::Float64, Ftest::Bool) where{T<:Number,P<:Number}
    N, M, M2 = Multitaper._pregap(t, x, nz)
    freq = collect(range(-1,1, length=M+1)*beta)
    # freq = range(-0.5,0.5,length=M+1)
    params = Multitaper.MTParameters(N*W, K, N, 1.0, M, 1, nothing)
    #
    eigenc(j,fr,x) = mapslices(slep -> slow_nfft_adj(t, slep .* x, beta, nfft = M, nfft2 = M2)[j], Multitaper.gpss(W, K, t, fr, beta=beta)[2], dims=1)
    eco = EigenCoefficient(mapreduce(j -> eigenc(j, freq[j+Int(M/2)], x), vcat, 1:(Int(M/2))), nothing)   
    #
    jknifed = Multitaper.jknife(eco,nothing,:spec) 
    #
    if Ftest
        # freq   = range(-0.5,0.5,length=M+1)
        freq = range(-beta, beta, length=M+1)
        gpsw0  = mapreduce(fr -> sum(Multitaper.gpss(W, K, t, fr, beta=beta)[2],dims=1),vcat, freq[1:Int(M/2)])
        gpsw0sq= sum(abs2, gpsw0, dims=2)
        # println("Array dimensions: ", size(eco.coef), " ", size(gpsw0), " ", size(gpsw0sq))
        mu    = sum(broadcast(/, eco.coef.*gpsw0, gpsw0sq), dims = 2) # dims = 2 means sum over k 
        num   = real.((K - 1)*abs2.(mu).*gpsw0sq)
        denom = real.(sum(abs2, (eco.coef .- broadcast(*, mu, gpsw0)), dims = 2))
        Ft    = vec((num./denom))
        Fpval = 1.0 .- cdf.(FDist.(2,2*K .- 2), Ft)
    else
        Fpval = nothing
    end
    #
    return Multitaper.MTSpectrum(freq[(Int(M/2)+1):M], mean(abs2.(eco.coef), dims=2)[:], nothing, params, eco, Fpval, jknifed[2], nothing)
end
