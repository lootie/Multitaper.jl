

# GNU GPL v2 licenced to C. Haley and C. Geoga 04/2020

"""
The underlying Kernel function, equation (54c) in Simons & Wang 2011 (SW2011)
Inputs: x - a spatial point
        y - another spatial point
        Kp - the radius of the disk in spectral space
Outputs: D(x,y) - the spacelimited kernel
"""
function Dfun(x::Vector{Float64}, y::Vector{Float64}, Kp::Float64)
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
function Dmatrix(Kp::Float64, pts1::Vector{Vector{Float64}}, 
                 pts2::Vector{Vector{Float64}})::Matrix{Float64}
  Out = zeros(length(pts1), length(pts2))
  @simd for j in eachindex(pts1)
    @simd for k in eachindex(pts2)
      @inbounds Out[j,k] = Dfun(pts1[j], pts2[k], Kp)
    end
  end
  return Out
end

"""
Slepian functions concentrated in 2 dimensions on a rectangle in 
physical space and a circle in spectral space. 
Inputs: nslep - number of output slepians
        n - number of GL nodes in the x-direction
        m - number of GL nodes in the y-direction
        Kp - the radius of the circle in spectral space
        N - number of Gauss-Legendre nodes in the first dimension
        M - number of Gauss-Legendre nodes in the second dimension
        verbose - select true if you would like to see the concentrations.
Outputs:
        sleps - an array of 2D tapers 
"""
function rectsleps(nslep::Int64, n::Int64, m::Int64, Kp::Float64,
                   N::Int64, M::Int64; verbose::Bool=false)::Vector{Matrix{Float64}}

  # Get the quadrature weights and nodes for each dimensions:
  no1, wt1 = FastGaussQuadrature.gausslegendre(N)
  no2, wt2 = FastGaussQuadrature.gausslegendre(M)
  no       = collect.(collect(Iterators.product(no1, no2))[:])
  wtv      = prod.(collect(Iterators.product(wt1, wt2))[:])
  
  # set up the eigenvalue problem and factorize: see eqn (86) of SW2011
  Kf       = Dmatrix(Kp, no, no)
  W        = Diagonal(sqrt.(wtv))
  solvme   = Symmetric(W*Kf*W) 
  factd    = eigen(solvme) # formerly eigfact

  # Extract the slepians, show the concentrations if verbose:
  goodind  = sortperm(factd.values, rev=true)[1:nslep]
  if verbose
    println("The $(nslep) concentrations:")
    for j in 1:nslep
      println(factd.values[goodind[j]])
    end
  end

  # Get the number of requested sleps: see eqn (88) of SW2011
  points   = collect.(collect(Iterators.product(LinRange(-1.0, 1.0, n), LinRange(-1.0, 1.0, m)))[:])
  sleps    = Matrix{Float64}[]
  for l in 1:nslep
    newslep = zeros(Float64, n*m)
    @simd for j in eachindex(points)
      @simd for k in eachindex(no)
        @inbounds newslep[j] += wtv[k]*Dfun(no[k], points[j], Kp)*factd.vectors[:,goodind[l]][k]
      end
    end
    push!(sleps, reshape(newslep, n, m)./factd.values[goodind[l]])
  end

  return sleps
end

""" 
Function for computing the 2D multitaper spectrum estimate

Inputs: Md - data matrix
        K - squared radius (bandwidth) of the region of concentration in spectrum
        ntaper - number of tapers
        padn - number of samples to pad to in the x-direction
        padm - number of samples to pad to in the y-direction
        nquad - number of quadrature nodes in the x-direction (default 72)
        mquad - number of quadrature nodes in the y-direction (default 72)
        center - boolean whether to subtract the mean or not (default true)
        sleps - if you have already computed the tapers, can use them as input here (default nothing)

Outputs: Matrix containing the 2D multitaper spectrum (simple averaging used)

If this function runs slowly for the size of your input matrix, reducing the number of quadrature
nodes for computing the tapers, or alternatively precomputing the tapers may help.
"""
function multispec2_Rectangle(Md::Matrix{Float64}, K::Float64, ntaper::Int64,
                              padn::Int64=0, padm::Int64=0, nquad::Int64=72, mquad::Int64=72; 
                              center::Bool=true, sleps::Union{Vector{Matrix{Float64}},Nothing}=nothing)
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

