

# GNU GPL v2 licenced to C. Haley and C. Geoga 04/2020

"""
The underlying Kernel function, equation (54c) in Simons & Wang 2011 (SW2011)
Inputs: x - a spatial point
        y - another spatial point
        Kp - the radius of the disk in spectral space
Outputs: D(x,y) - the spacelimited kernel
"""
function dfun(x, y, Kp)
  out   = 0.0
  if x == y
    out = abs2(Kp)/(4.0*pi)
  else
    out = Kp*besselj1(Kp*norm(x-y))/(2.0*pi*norm(x-y))
  end
  return out
end

"""
The matrix with entries equal to the kernel evaluated on each of the 
points pts1, pts2 which are vectors.
"""
function dmatrix(Kp, pts1, pts2)
  Out = zeros(length(pts1), length(pts2))
  @simd for j in eachindex(pts1)
    @simd for k in eachindex(pts2)
      @inbounds Out[j,k] = dfun(pts1[j], pts2[k], Kp)
    end
  end
  return Out
end

"""

    rectsleps(nslep,n,m,Kp,N,M, <keyword arguments>)=)

Slepian functions concentrated in 2 dimensions on a rectangle in 
physical space and a circle in spectral space. 

...
# Arguments
 - `nslep::Int64`: number of output slepians
 - `n::Int64`: number of GL nodes in the x-direction
 - `m::Int64`: number of GL nodes in the y-direction
 - `Kp::Float64`: the radius of the circle in spectral space
 - `N::Int64`: number of Gauss-Legendre nodes in the first dimension
 - `M::Int64`: number of Gauss-Legendre nodes in the second dimension
 - `verbose::Boo`: select true if you would like to see the concentrations
...

...

# Outputs

 - `sleps::Array{Matrix{Float64},1}` - an array of 2D tapers 

...
"""
function rectsleps(nslep, n, m, Kp, N, M; verbose = false)

  # Get the quadrature weights and nodes for each dimensions:
  no1, wt1 = FastGaussQuadrature.gausslegendre(N)
  no2, wt2 = FastGaussQuadrature.gausslegendre(M)
  no       = collect.(vec(collect(Iterators.product(no1, no2))))
  wtv      = prod.(vec(collect(Iterators.product(wt1, wt2))))
  
  # set up the eigenvalue problem and factorize: see eqn (86) of SW2011
  Kf       = dmatrix(Kp, no, no)
  W        = Diagonal(sqrt.(wtv))
  solvme   = Symmetric(W*Kf*W) 
  factd    = eigen(solvme) # formerly eigfact

  # Extract the slepians, show the concentrations if verbose:
  highconcind  = sortperm(factd.values, rev=true)[1:nslep]
  if verbose
    println("The $(nslep) concentrations:")
    for j in 1:nslep
      println(factd.values[highconcind[j]])
    end
  end

  # Get the number of requested sleps: see eqn (88) of SW2011
  points   = collect.(vec(collect(Iterators.product(range(-1.0, 1.0, length=n),
                     range(-1.0, 1.0, length = m)))))
  sleps    = Matrix{Float64}[]
  for l in 1:nslep
    newslep = zeros(Float64, n*m)
    @simd for j in eachindex(points)
      @simd for k in eachindex(no)
        @inbounds newslep[j] += wtv[k]*dfun(no[k], points[j],
                                Kp)*factd.vectors[:,highconcind[l]][k]
      end
    end
    push!(sleps, reshape(newslep, n, m)./factd.values[highconcind[l]])
  end

  return sleps
end

"""
    multispec2_Rectangle(Md, K, ntaper; <keyword arguments>)

Function for computing the 2D multitaper spectrum estimate

...
# Arguments
 - `Md::Matrix{T} where T<:Float64`: data matrix 
 - `K::Float64`: squared radius (bandwidth) of the region of concentration in spectrum
 - `ntaper::Int64`: number of tapers
 - `padn::Int64 = 0`: number of samples to pad to in the x-direction 
 - `padm::Int64 = 0`: number of samples to pad to in the y-direction
 - `nquad::Int64 = 72`: number of quadrature nodes in the x-direction
 - `mquad::Int64 = 72`: number of quadrature nodes in the y-direction
 - `center::Bool = true`: whether to subtract the mean or not
 - `sleps::Union{Matrix{T},Nothing} = nothing`: if you have already computed the tapers, input them here
...

...
# Outputs

 - Matrix containing the 2D multitaper spectrum (simple averaging used)

Note: If this function runs slowly for the size of your input matrix, reducing the number
of quadrature nodes for computing the tapers, or alternatively precomputing the
tapers may help.
...

"""
function multispec2_Rectangle(Md, K, ntaper, padn=0, padm=0, nquad=72, 
                              mquad=72; center=true, sleps=nothing)

  M = complex(Md)
  if center
    M .-= mean(M)
  end
  if sleps == nothing
    sleps = rectsleps(ntaper, size(M)..., K, nquad, mquad)
  end
  szs   = size(M)
  padsz = (max(size(M,1), padn), max(size(M,2), padm))
  Out   = zeros(Float64, padsz...)
  tmp   = zeros(ComplexF64, padsz...)
  pln   = plan_fft!(tmp)
  for j in eachindex(sleps)
    fill!(tmp, zero(ComplexF64))
    setindex!(tmp, sleps[j].*M, 1:szs[1], 1:szs[2])
    pln*tmp
    Out .+= abs2.(tmp)
  end
  Out .*= 1.0/Float64(ntaper)
  return fftshift(Out)
end

