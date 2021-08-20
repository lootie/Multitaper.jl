# All of these windows are reproduced from F.J. Harris except for the R default window #

""" Lag grid """
function lag(N::Int64)
  N2 = Int64(floor(N/2))
  return isodd(N) ? ((-N2):N2) : ((-N2):(N2-1))
end

""" Frequency grids """
function fgrid(N::Int64; dt::T=1.0) where T<:Number
  return dt*LinRange(-0.5,0.5,N+1)[1:N]
end

""" V.A. Rectangle (Dirichlet) window """
function dirichlet(N::Int64)
  ones(N)/sqrt(N)
end

""" V.B. Triangle (Bartlett) window """
function bartlett(N::Int64)
  # N is odd here
  N2 = isodd(N) ? Int64(floor(N/2)) : Int64(N/2)
  bart = 1.0 .- abs.(collect(lag(N))/N2)
end

""" R default is just ones with a 10% cosine taper, and is the default from R
spec.pgram """
function default(N::Int64, p::Union{Float64} = 0.1)
  if ((p .< 0) || (p .> 0.5))
      error("'p' must be between 0 and 0.5")
  end
  m = floor(N * p)
  if (m == 0) 
    # Not recommended cause it usually produces a bad result, but does not throw an error
    y = ones(N)
  else 
    w = 0.5 * (1 .- cos.(pi * collect(1:2:(2*m-1))/(2 * m)))
    y = vcat(w, ones(Int64(N - 2 * m)), w[end:-1:1]) 
  end
  return y
end

""" V.C. Cos^α(x) windows """
function cosalpha(N::Int64,alpha::Float64)
  return cos.(pi*fgrid(N)).^(alpha)
end

""" The cosine squared, or raised cosine, or Hanning"""
function hanning(N::Int64)
  return cosalpha(N,2.0)
end

# This one looks bad
""" V.D. Hamming window """
function hamming(N::Int64, alpha::Float64)
  return alpha .+ (1-alpha)*cos.(fgrid(N)*2*pi)
end

""" V.E. Blackman Window (32) """
function blackman(N::Int64)
  out = zeros(N)
  a = [0.42659071, 0.49656062, 0.07684867]
  for m = 0:(length(a)-1)
    out .+= a[m+1]*cos.(2*pi*fgrid(N)*m)
  end
return out
end

""" V.E. Blackman-Harris window (33) """
function blackmanharris(N::Int64, a0::Float64, a1::Float64, a2::Float64, a3::Float64)
  out = a0 .- a1*cos.(2*pi*fgrid(N)) .+ a2*cos.(2*pi*2*fgrid(N)) .- a3*cos.(2*pi*3*fgrid(N))
  return out
end

""" V.F. Parzen (Riesz) window """
function riesz(N::Int64)
  return 1.0 .- abs.(2.0*fgrid(N)).^2
end

""" V.F. Riemann window  (36) """
function riemann(N::Int64)
  out = sin.(2*pi*fgrid(N))
  out ./= fgrid(N)*2*pi
  # By Calculus, at zero we have 1.0
  N2 = isodd(N) ? Int64((N-1)/2) : Int64(N/2)
  out[N2+1] = 1.0
  return out
end

""" de la Valle Poussin (Jackson, Parzen) window """
function delaValle(N::Int64)
  N4 = Int64(floor(N/4))
  out = 2*(1.0 .- abs.(fgrid(N))*2).^3
  mid = 1.0 .- 6*(2*fgrid(N)).^2 .* (1.0 .- abs.(2*fgrid(N)))
  return vcat(out[1:N4], mid[(N4+1):(end-N4-1)], out[(end-N4):end])
end 

""" Tukey window is a cosine lobe of width (alpha/2)N convolved with a rectangle window of width
(1.0 -alpha/2)N """
function tukey(N::Int64, alpha::Float64)
  N2 = Int64(round(N/2))
  aN2 = Int64(round(alpha*N/2))
  out = ones(N2)
  for i = (aN2+1):N2
    out[i] = 0.5*(1.0+cos(pi*(-i-alpha*N/2)/(2*(1-alpha)*(N/2))))
  end
  return vcat(out[end:-1:1], out)
end

""" Poisson window is a two-sided exponential """
function poisson(N::Int64, alpha::Float64)
  return exp.(-2*alpha*fgrid(N))
end

""" Hanning-Poisson window """
function hanningpoisson(N::Int64, alpha::Float64)
  return 0.5*(1.0 .+ cos.(pi*2*fgrid(N))).*exp.(-alpha*2*fgrid(N))
end

""" Cauchy (abel, poisson) window """
function cauchy(N::Int64,alpha::Float64)
  return 1.0 ./ (1.0 .+ (alpha*fgrid(N)).^2)
end

""" Gaussian or weierstrass window """
function gaussWeier(N::Int64, alpha::Float64)
  return exp.(-0.5*(alpha*fgrid(N)).^2)
end

""" Kaiser-Bessel window """
function kaiserbessel(N::Int64, alpha::Float64)
  out = SpecialFunctions.besseli.(0.0, pi*alpha*sqrt.(1.0 .- fgrid(N).^2))
  out /= SpecialFunctions.besseli(0.0, pi*alpha)
  return out
end

""" Helper function that selects one of the windows for you """
function getwindow(N::Int64, taper::Union{Nothing, Symbol, Tuple{Symbol, Nothing}, Tuple{Symbol, Vector{Float64}}, Tuple{Symbol, Float64}})
  nom, parms = (typeof(taper) == Symbol) ? (taper, nothing) : taper   
  if (parms == nothing)
    if (nom == :cosalpha)||(nom == :hamming)||(nom == :tukey)||(nom == :poisson)||(nom == :hanningpoisson)||(nom == :cauchy)||(nom == :gaussweierstrass)||(nom == :kaiserbessel)
      error("You must supply the parameter for the taper")
    elseif (nom == :blackmanharris)
      error("You must supply a vector of four parameters for the taper.")
    end
  else
    if (nom == :bartlett)||(nom == :hanning)||(nom == :blackman)||(nom == :riesz)||(nom == :riemann)||(nom == :delavalle)
      error("The type of taper selected does not require parameters.")
    end
  end
  if nom == :default
    if parms == nothing
      win = default(N)
    else
      win = default(N,parms[1])
    end
  elseif nom == :dirichlet
    win = dirichlet(N)
  elseif nom == :bartlett
    win = bartlett(N)
  elseif nom == :cosalpha
    win = cosalpha(N, parms[1])
  elseif nom == :hanning
    win = hanning(N)
  elseif nom == :hamming
    win = hamming(N, parms[1])
  elseif nom == :blackman
    win = blackman(N)
  elseif nom == :blackmanharris
    win = blackmanharris(N, parms[1:4]...)
  elseif nom ==  :riesz
    win = riesz(N)
  elseif nom == :riemann
    win = riemann(N)
  elseif nom == :delavalle
    win = delaValle(N)
  elseif nom == :tukey
    win = tukey(N, parms[1])
  elseif nom == :poisson
    win = poisson(N, parms[1])
  elseif nom == :hanningpoisson
    win = hanningpoisson(N, parms[1])
  elseif nom == :cauchy
    win = cauchy(N, parms[1])
  elseif nom == :gaussweierstrass
    win = gaussWeier(N, parms[1])
  elseif nom == :kaiserbessel
    win = kaiserbessel(N, parms[1])
  else
    error("window name must be one of :default, :dirichlet, :bartlett, :cosalpha, :hanning, :hamming,
          :blackman, :blackmanharris, :riesz, :delavalle, :tukey, :poisson, :hanningpoisson, :cauchy,
          :gaussweierstrass, or :kaiserbessel")
  end
  return win
end

""" Helper function for the plotter to get a string describing the taper """
function taperstring(taper::Union{Nothing, Symbol, Tuple{Symbol, Nothing}, Tuple{Symbol, Vector{Float64}}, Tuple{Symbol, Float64}})
  nom, parms = taper
  if (parms == nothing)
    if (nom == :cosalpha)||(nom == :hamming)||(nom == :tukey)||(nom == :poisson)||(nom == :hanningpoisson)||(nom == :cauchy)||(nom == :gaussweierstrass)||(nom == :kaiserbessel)
      error("You must supply the parameter for the taper")
    elseif (nom == :blackmanharris)
      error("You must supply a vector of four parameters for the taper")
    end
  end
  if nom == :default
    win = "10 percent cosine taper" 
  elseif nom == :dirichlet
    win = "Dirichlet taper"
  elseif nom == :bartlett
    win = "Bartlett taper"
  elseif nom == :cosalpha
    win = "Cos^{α} taper with parameter $(parms[1])"
  elseif nom == :hanning
    win = "Hanning taper"
  elseif nom == :hamming
    win = "Hamming taper with parameter $(parms[1])"
  elseif nom == :blackman
    win = "Blackman taper"
  elseif nom == :blackmanharris
    win = "Blackman-Harris taper with parameters "*join(parms[1:4],", ")
  elseif nom ==  :riesz
    win = "Riesz taper"
  elseif nom == :riemann
    win = "Riemann taper"
  elseif nom == :delavalle
    win = "De-la-Valle taper"
  elseif nom == :tukey
    win = "Tukey taper with parameter $(parms[1])"
  elseif nom == :poisson
    win = "Poisson taper with parameter $(parms[1])"
  elseif nom == :hanningpoisson
    win = "Hanning-Poisson taper with parameter $(parms[1])"
  elseif nom == :cauchy
    win = "Cauchy taper with parameter $(parms[1])"
  elseif nom == :gaussweierstrass
    win = "Gauss-Weierstrass taper with parameter $(parms[1])"
  elseif nom == :kaiserbessel
    win = "Kaiser-Bessel taper with parameter $(parms[1])"
  else
    error("window name must be one of :default, :dirichlet, :bartlett, :cosalpha, :sinalpha, :hanning, :hamming,
          :blackman, :blackmanharris, :riesz, :delavalle, :tukey, :poisson, :hanningpoisson, :cauchy,
          :gaussweierstrass, or :kaiserbessel")
  end
  return win
end

