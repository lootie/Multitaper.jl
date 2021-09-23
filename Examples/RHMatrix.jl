
module RHMatrix

  using  KernelMatrices, LinearAlgebra
  import KernelMatrices: KernelMatrix, full, submatrix, ACA

  export rhodlr

  abstract type LowRankFact end

  struct UVt{T} <: LowRankFact
    U::Matrix{T}
    V::Matrix{T}
  end

  @inline matrix(UV::UVt{T}) where{T} = UV.U*UV.V'

  @inline Base.size(UV::UVt{T}) where{T} = (size(UV.U, 1), size(UV.V, 2))

  @inline Base.size(UV::UVt{T}, j) where{T} = size(UV)[j]

  @inline Base.:*(UV::UVt{T}, v::Vector{T}) where{T} = UV.U*(UV.V'v)

  mutable struct RHODLR{T, LF<:LowRankFact}
    A11::Union{Matrix{T}, RHODLR{T}}
    A22::Union{Matrix{T}, RHODLR{T}}
    A12::LF
    A21::LF
  end

  @inline Base.size(R::RHODLR) = size(R.A11) .+ size(R.A22)
  @inline Base.size(R::RHODLR, j::Int64) = size(R.A11, j) + size(R.A22, j)

  function rhodlr(K::KernelMatrix{T,N,A,Fn}, lvl::Int64, 
                  tol=1.0e-5, maxrank=0) where{T,N,A,Fn}

    # Check sizes:
    iszero(lvl) && return K
    len1, len2 = size(K)
    mpt1, mpt2 = Int64(floor(len1/2)), Int64(floor(len2/2))

    # Compute the sub-blocks:
    K11 = submatrix(K, 1,      mpt1, 1,      mpt2)
    K22 = submatrix(K, mpt1+1, len1, mpt2+1, len2)
    K12 = submatrix(K, 1,      mpt1, mpt2+1, len2)
    K21 = submatrix(K, mpt1+1, len1, 1,      mpt2)

    # Factorize the off-diagonal blocks:
    A12 = UVt{T}(ACA(K12, tol, maxrank)...)
    A21 = UVt{T}(ACA(K21, tol, maxrank)...)

    if isone(lvl)
      return RHODLR{T, UVt{T}}(full(K11), full(K22), A12, A21)
    else
      A11 = rhodlr(K11, lvl-1, tol, maxrank)
      A22 = rhodlr(K22, lvl-1, tol, maxrank)
      return RHODLR{T, UVt{T}}(A11, A22, A12, A21)
    end
  end

  @inline matrix(M::Matrix{T}) where{T} = M

  matrix(R::RHODLR) = [matrix(R.A11) matrix(R.A12) ; matrix(R.A21) matrix(R.A22)]

  function Base.:*(R::RHODLR, V::Vector)
    ix1 = 1:size(R.A11, 2)
    ix2 = (size(R.A11, 2)+1):size(R, 2)
    return vcat(R.A11*V[ix1] + R.A12*V[ix2], R.A21*V[ix1] + R.A22*V[ix2])
  end

end

