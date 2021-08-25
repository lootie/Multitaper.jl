function maternsdf(f::Number, phi, a, nu, flr)
  return phi*(abs2(a) + abs2(2.0*pi*f))^(-nu - 0.5) + flr
end

""" Simulate a Matern process of length N """
function Matern_sim(p, N)
  R = real(ifftshift(ifft(fftshift(sqrt.(maternsdf.(LinRange(-0.5,0.5,N+1)[1:N], p...))))))*2
  x = Multitaper.conv(randn(N),R)[Int64(N/2-1):Int64(3*N/2-2)]
  sigmasq = var(x)
  return x, sigmasq
end

function dbnll(θ::Vector, MTS::Vector, SpecWin::Vector)
  fmat = maternsdf.(LinRange(-0.5,0.5,2*length(MTS)-1)[1:(2*length(MTS)-2)], θ...)
  fb = circ_conv(SpecWin, fmat/N)[end:-1:1]
  return whitlik(fb, MTS/4)
end

# The objective function in the format NLopt wants.
function dbobjective(θ::Vector, grad::Vector, countr::Vector)
  countr[1] += 1
  if length(grad) > 0
    ForwardDiff.gradient!(grad, db_innerfun, θ) 
  end
  return db_innerfun(θ)
end
