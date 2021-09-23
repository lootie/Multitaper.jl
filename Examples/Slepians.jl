# The underlying Kernel function, equation 54c in Simons & Wang 2011.
function dfn(x, y, Kp)
  x == y && return abs2(Kp[])/(4.0*pi)
  Kp[]*besselj1(Kp[]*sqrt(sum(abs2, x.-y)))/(2.0*pi*sqrt(sum(abs2, x.-y)))
end


function getnodeswts(szs)
    nowtv = [FastGaussQuadrature.gausslegendre(nj) for nj in szs]
end

function givewts(nowtv)
    return sqrt.(prod.(vec(collect(IterTools.product([x[2] for x in nowtv]...)))))
end

"""
**ARGS: 

-M    (number of tapers)

-Kp   (A concentration measure like NW in the 1D case)

-szs  (An AbstractArray of sizes for the QUADRATURE problem)

**KWARGS:

-ex   (EXact flag, as opposed to using H-matrices to approximate)

-p    (Precision for HODLR matrix approximation)

-lvl   (the HODLR LeVel if not computing them exactly)

-maxrank (maximum allowed rank for off-diagonal blocks in HODLR approximation)

So, for example, if I had 4D data and I wanted to use (32, 34, 36, 38) quadrature
nodes in each of those dimensions (for some reason), and I wanted to use H-matrices
because otherwise I'd run out of RAM, I would perhaps set a fixed maximum rank of
256 for my off-diagonal blocks and say

````
sleps = customsleps(5, 10.0, (32, 35, 36, 38), maxrank=256)
````
As a word of warning, though, the actual HODLR structure of the kernel matrix will
really weaken quickly as dimension increases. So these objects will not necessarily
be particularly slepian-like, and at some point may not actually be that much better
than outer-product 1D slepians or something that is far less painful to compute. 
So do apply caution and perform sanity checks.

If you do not supply your own nodes, the region of interest in space is assumed to 
be rectangular. If you want something other than this, you'll provide your own
gauss-legendre nodes in the no and sqwt field.
"""
function customsleps(M, Kp, szs; prec=1.0e-8, exact=false, lvl=6, maxrank=0, no = nothing, sqwt = nothing,
    int = nothing)
  if typeof(Kp) <: Float64 
    return customsleps(M, [Kp], szs, prec=prec, exact=exact, lvl=lvl, maxrank=maxrank, no = no, sqwt = sqwt)
  end
  length(szs) > 1 || @warn("Don't use this function for 1D slepians...")

  if (no == nothing)||(sqwt == nothing)
      # Create quadrature nodes and weights for each dimension:
      nowtv = getnodeswts(szs)

      # Re-format the above nodes and weights under a single index:
      no    = vec(collect(IterTools.product([x[1] for x in nowtv]...)))
      sqwt  = givewts(nowtv)
  end

  # Create the kernel matrix, exactly or via HODLR approximation:
  _K    = KernelMatrix(no, no, Kp, dfn)
  K     = exact ? full(_K) : RHMatrix.rhodlr(_K, lvl, prec, maxrank)

  # Solve the eigenvalue problem to obtain the slepians at the quad nodes:
  s     = eigsolve(z->sqwt.*(K*(sqwt.*z)), prod(szs), M, :LM, issymmetric=true)

  # Prepare even grid points, compute slepians at those points:
  evpts = vec(collect(product([range(-1.0, 1.0, length=s) for s in ((int == nothing) ? szs : int)]...)))
  _K2   = KernelMatrix(evpts, no, Kp, dfn)
  K2    = exact ? full(_K2) : RHMatrix.rhodlr(_K2, lvl, prec, maxrank)
  sleps = [sqwt.*(K2*x[2])./x[1] for x in zip(s[1], s[2])]

  return s[1], [reshape(sp, szs...) for sp in sleps]
end