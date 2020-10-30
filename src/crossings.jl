
# Details about these functions can be found in Thomson and Haley Proc. Roy Soc. Ser
# A, 2014 In that paper, numbers of upcrossings of a given significance level of a
# multitaper spectrum are given, as well as the expected dwell bandwidth of large
# excursions.

""" Probability distribution of a Chi^2/(2α) distribution """
p(z, α) = pdf(Gamma(α, 1.0/α), z)

""" Intermediate function to determining the constant psi """
function psi_guts(N, CR, K)
  v   = dpss_tapers(N,CR,K,:tap)
  λ   = dpss_tapers(N,CR,K,:egval)
  C   = [(i - j) for i in 1:N, j in 1:N]
  ups = sum((mapreduce(x->x*transpose(x), +, 
              (sqrt(λ[k])*v[:,k] for k in 1:K)).*C).^2)
  return (2*pi/N)*sqrt(ups/sum(abs2.(λ)))
end

"""
    MT_Upcrossings(z, α, CR, N)

Number of upcrossings of the level z of a multitaper spectrum per Rayleigh resolution
...
# Arguments
 - `z::Float64`: the level of the standardized spectrum
 - `α::Float64`: degrees of freedom of estimate
 - `CR::Float64: time bandwidth product
 - `N::Int64`: sampling rate in time units 
...

...
# Outputs
 - `Float64`: the number of upcrossings per Rayleigh
...

"""
MT_Upcrossings(z, α, CR, N) = psi_guts(1024, CR, α)*sqrt(z/(2*pi*α))*p(z,α)

""" The value of a multitaper estimate of a white noise spectrum over which sig
proportion of the spectrum is larger, when the MT spectrum has α degrees of freedom
"""
z(α, sig) = quantile.(Gamma(α,1/α), sig)

"""
    uctable(α, CR, <keyword arguments>)

Gives an upcrossing table for α = K, CR = NW (Float), optionally num_Ray has
default 1e5, and signficiance levels that can be chosen, or set to default

...
# Arguments

## Positional Arguments

 - `α::Int64`: the number of tapers

 - `CR::Int64`: the time bandwidth product

## Keyword Arguments

 - `num_Ray::Int64 = 10_000`: Number of Rayleigh resolutions

 - `sig::Union{Vector{Float64},Float64,Nothing}=nothing`: Significance levels
...

...
# Outputs

 - `Matrix{Float64}`: the upcrossing table

 - `Array{String,1}`: the header for the table
...

"""

""" Gives an upcrossing table for α = K, CR = NW (Float), optionally num_Ray has
default 1e5, and signficiance levels that can be chosen, or set to default """
function uctable(α, CR; num_Ray=100000, sig=nothing, N_segments=1) 
  s   = (sig == nothing) ? 1.0 .- [0.5, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6] : sig
  Z   = z.(α, s)
  U   = MT_Upcrossings.(Z, α*N_segments, CR, 1024)
  return (hcat(100*s, z(α, s), num_Ray*U, (1.0 .- s)./U), 
          ["P%", "z", "U(z,α)", "D(z,α)"]) 
end

"""
    Pgram_upcrossings(z, N)

Number of upcrossings of the level z of a periodogram spectrum per Rayleigh resolution
...
# Arguments
 - `z::Float64`: the level of the standardized spectrum
 - `N::Int64`: sampling rate in time units 
...

...
# Outputs
 - `Float64`: the number of upcrossings per Rayleigh
...

"""
Pgram_upcrossings(z, N) = sqrt((π*z/3)*(N^2 - 1))*exp(-z)/N

""" Dwell bandwidth for the periodogram, in Rayleighs """
Pgram_dwellband(z) = sqrt(3/(z*π))


