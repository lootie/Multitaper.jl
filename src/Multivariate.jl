2
""" T squared test for simultaneous line components """ 
function Tsqtest_pval(dcs::Vector{Float64},tap::ecoef)
  ell, k    = (Float64(size(tap.coef,1)), Float64(size(tap.coef,2)))
  if (ell > k)
    println("Tsq test error: you have K = $k tapers and L = $ell eigencoefficients, but L must be less than K.")
    return NaN
  else
    sum2dc  = sum(abs2, dcs) 
    # Following Thomson (2011) "Some problems ... cyclostationary data"
    Q       = tap.coef*dcs/sum2dc 
    # vector of mean estimates (15), l of them
    Rt      = copy(transpose(tap.coef)) - kron(transpose(Q),dcs) 
    # tap1 is ell by k, Q is 1 by ell, dc is k by 1
    Tsq     = ((k-ell)/ell)*sum2dc*real.(transpose(Q)*pinv(transpose(Rt)*Rt)*Q)
    pval    = StatsFuns.fdistccdf(2*ell,2*(k-ell),Tsq)
    return pval
  end
end

""" Computes multitaper cross-spectrum or coherence with a handful of extra gadgets. """
function multispec(S1::Union{Vector{T},ecoef}, S2::Union{Vector{T},ecoef}; outp::Symbol = :coh, 
                   NW::Real = 4.0, K::Int = 6, offset::Union{Float64,Int64} = 0, 
                   dt::Float64=1.0, ctr::Bool = true, 
                   pad::Union{Int,Float64} = 1.0, 
                   dpVec::Union{Vector{Float64},Matrix{Float64},Nothing} = nothing,
                   guts::Bool = false, 
                   jk::Bool = false, 
                   Tsq::Union{Vector{Float64},Vector{Vector{Float64}},Vector{Int64},Vector{Vector{Int64}},Nothing}=nothing, 
                   alph::Float64 = 0.05) where T<:Number
  
  if (typeof(S1) == ecoef) && (typeof(S2) == ecoef)
    coefswts = [S1, S2]
    halffreq = size(S1.coef,1)
    fftleng  = 2*halffreq - 1 # This is only for real data
    lengt    = Int64(round(fftleng/pad))
  elseif ((typeof(S1) == Vector{Float64}) && (typeof(S2) == Vector{Float64})) || ((typeof(S1) == Vector{ComplexF64}) && 
         (typeof(S2) == Vector{ComplexF64}))
    # Make sure vectors are compatible
    if (length(S1) != length(S2))
      error("Vectors must be the same length")
    end
    # Compute the lengths of the spec obj
    lengt, fftleng, halffreq = output_len(S1,pad)
    # Compute the array of dpss vectors if they aren't given.
    if (dpVec == nothing) 
      dpVec = dpss_tapers(lengt, NW, K)
      dpVec .*= sqrt(dt)
    end 
    coefswts = map(x->mtspec_guts(x,dpVec,fftleng,halffreq,var(x),ctr,nothing),[S1, S2])
  end

  # If offset cross-spectrum or coherence is chosen
  if offset != 0
    offset = freq_to_int([offset], lengt, dt)[1]
    if offset > 0
      ind = vcat(collect(1:halffreq), collect((halffreq-1):-1:2))[offset:(halffreq+offset-1)]
    elseif offset < 0 
      ind = vcat(collect((halffreq):-1:2), collect(1:halffreq))[(halffreq+offset):(2*halffreq+offset-1)]
    end
    # Shift the first set of eigencoefficients forward or backwards accordingly
    coefswts[1] = ecoef(coefswts[1].coef[ind,:], 
                        (coefswts[1].wts == nothing ? nothing : coefswts[1].wts[ind,:]))
    foffset = int_to_freq([offset],lengt,dt)[1]
    freq = (1/dt)*LinRange(foffset, 0.5 + foffset, halffreq)  
  elseif offset == 0
    freq = (1/dt)*LinRange(0, 0.5, halffreq)
  end

  # Output is spectrum, coherence, or transfer function (outp = :spec, :coh or :transf)
  S, jv    = jk ? jknife(coefswts...,outp) : (jknife(coefswts...,outp)[1], nothing)
  S       .*= (outp == :cross) ? dt : 1.0
  # To compute the transfer functions, we need to use the following helper code
  if outp == :transf
    B = mapreduce(f -> transpose(mttransfn(coefswts[1].coef[f,:],coefswts[2].coef[f,:])[1]),vcat, 1:length(S))
  end
  # Now we need to get the phase estimate 
  phase, jphase = jk ? jknife_phase(coefswts...) : (jknife_phase(coefswts...)[1], nothing)

  # Compute the T^2- test if requested:
  dcs      = (Tsq != nothing) ? map(isodd,1:K).*sum(dpVec,dims=1)[:] : nothing

  # Do the Tsquared test
  if typeof(Tsq) != Nothing
    Tsq      = (typeof(Tsq) <: Vector{Number}) ? [Tsq] : Tsq
    map!(x -> freq_to_int(Tsq[x], lengt, dt), Tsq, eachindex(Tsq))
    Tsq = Vector{Vector{Int64}}(Tsq)  
    if (2*K < (true ? 1 : 2)*maximum(length.(Tsq)))
      error("There are too few tapers for the number of Tsq tests.")
    end
    Tv = map(x->Tsqtest_pval(dcs,
             ecoef(vcat(coefswts[1].coef[Tsq[x],:],coefswts[2].coef[Tsq[x],:]),nothing)),eachindex(Tsq)) 
  else
    Tv = nothing
  end

  params = mtparams(NW, K, lengt, dt, fftleng, 1, nothing) # NW, K, N, dt, M, nsegments, overlap
  if outp == :spec
    return mtspec(freq, S, phase, params, (guts ? coefswts : nothing), 
                nothing, (jk ?  [jv, jphase] : nothing), Tv) 
  elseif outp == :coh
    return mtcoh(freq, S, phase, params, (guts ? coefswts : nothing), 
                     (jk ?  [jv,jphase] : nothing), Tv) 
    # frequency, spec/crosspec, coef & weights, jackknife, Tsq test.
  elseif outp == :transf
    return mttransf(freq, abs.(B).^2, unwrapphase(angle.(B[:])*180/pi,:deg), params, (guts ? coefswts : nothing), nothing)
  end

end

""" Computes univariate multitaper cross-covariance/cross-correlation function. Inputs a mtcoh or mtceps struct. """
function mt_ccvf(S::mtspec; typ::Symbol = :ccvf)   
  if typeof(S) == mtcoh
    lags = S.params.dt*S.params.N*LinRange(-1.0, 1.0, length(S.coh))
    if typ == :ccvf
      error("Cannot compute cross covariance from coherence.")
    elseif typ == :ccf
      return mtccf(lags, fftshift(real.(ifft(S.coh))), S.params)
    end
  elseif typeof(S) == mtspec
  lags = S.params.dt*S.params.N*LinRange(-1.0, 1.0, length(S.S))
    if typ == :ccvf
      return mtccvf(lags, fftshift(real.(ifft(S.S))), S.params)
    elseif typ == :ccf
      ccvf = real.(ifft(S.S))[1:length(S.S)]
      return mtccf(lags, fftshift(ccvf)/ccvf[1], S.params)
    end
  else
    error("Select one of :ccvf (cross covariance), :ccf (cross correlation) for output")
  end
end

""" Computes bivariate multitaper cross-covariance/cross-correlation function if you haven't already computed the MT cross-spectrum. """
function mt_ccvf(S1::Vector{T}, S2::Vector{T}; typ::Symbol = :ccvf, NW::Real = 4.0, K::Int = 6, 
                   dt::Float64=1.0, ctr::Bool = true, 
                   pad::Union{Int,Float64} = 1.0, 
                   dpVec::Union{Vector{Float64},Matrix{Float64},Nothing} = nothing,
                   guts::Bool = false, 
                   jk::Bool = false, 
                   Tsq::Union{Vector{Float64},Vector{Vector{Float64}},Vector{Int64},Vector{Vector{Int64}},Nothing}=nothing, 
                   alph::Float64 = 0.05) where T<:Number
  if (typ != :ccvf)&&(typ != :ccf)
    error("Output type must be one of cross-covariance (:ccvf) or cross-correlation (:ccf).")
  end
  S = multispec(S1, S2, outp = :spec, 
                NW = NW, K = K, dt = dt, ctr = ctr, pad = pad, dpVec = dpVec, guts = false, 
                jk = false, Tsq = nothing, alph = alph) 
  return mt_ccvf(S, typ = typ)
end

""" Multivariate version of the multispec call, data are in the columns of a matrix"""
function multispec(S1::Matrix{T}; outp::Symbol=:coh, 
                   NW::Real = 4.0, K::Int = 6, 
                   dt::Float64 = 1.0, ctr::Bool = true,
                   pad::Union{Int,Float64} = 1.0, 
                   dpVec::Union{Vector{Float64},Matrix{Float64},Nothing} = nothing,
                   guts::Bool = false, a_weight::Bool = true, 
                   jk::Bool = false, 
                   Ftest::Bool = false, 
                   Tsq::Union{Array{Int64,1},Array{Array{Int64,1},1},Nothing}=nothing,
                   alph::Float64 = 0.05) where T<:Number

  lengt, p = size(S1)

  if !((outp == :cross)||(outp == :coh) || (outp == :justspecs))
    error("The output must be one of :cross, :coh, or :justspecs")
  end        

  if p > 3
    println("You are computing $(Int64((p)*(p+1)/2 - p)) cross-spectra/coherences.")
  end

  # Compute the array of dpss vectors if they aren't given.
  if (dpVec == nothing) 
    dpVec   = dpss_tapers(lengt, NW, K)
    dpVec .*= sqrt(dt)
  end 
  
  # Get the spectra
  specs     = map(x->multispec(x, NW = NW, K = K, ctr = ctr, dt = dt, pad = pad, guts = true, a_weight = a_weight,
                           jk = jk, alph = alph, dpVec = dpVec, Tsq = Tsq, Ftest = Ftest), 
              (S1[:,k] for k=1:p))
 
  # Get the cross-spectra
  if (outp != :justspecs)
    crosspecs = (outp == :cross) ? Array{mtspec,2}(undef, p, p) : 
                                  Array{mtcoh,2}(undef, p, p)
    for x in CartesianIndex.(filter(x -> x[2]>x[1], Tuple.(eachindex(view(crosspecs,1:p,1:p)))))
      crosspecs[x] = multispec(specs[x[1]].coef, specs[x[2]].coef, outp = ((outp == :cross) ? :spec : :coh), 
                               NW = NW, K = K, 
                               ctr = ctr, dt = dt, pad = pad, guts = false, 
                               jk = jk, alph = alph, dpVec = dpVec, Tsq = Tsq)
    end
  else
    crosspecs = nothing
  end

  # Compute the T^2- test if requested:
  dcs = (Tsq != nothing) ? map(isodd,1:K).*sum(dpVec,dims=1)[:] : nothing
  # Do the Tsquared test
  if typeof(Tsq) != Nothing
    Tsq      = (typeof(Tsq) <: Vector{Number}) ? [Tsq] : Tsq
    map!(x -> freq_to_int(Tsq[x], lengt, dt), Tsq, eachindex(Tsq))
    Tsq = Vector{Vector{Int64}}(Tsq)  
    if (2*K < (true ? 1 : 2)*maximum(length.(Tsq)))
      error("There are too few tapers for the number of Tsq tests.")
    end
    Tv = zeros(length(Tsq)) 
    for x = eachindex(Tsq) 
        temp  = mapreduce(y -> y[Tsq[x],:], vcat, (specs[j].coef.coef for j = 1:p))
        Tv[x] = Tsqtest_pval(dcs, ecoef(temp,nothing))
    end
  else
    Tv = nothing
  end

  # If you wanted a vanilla multispec with no gadgets, you will get it in the simplest form 
  if (Tsq != nothing) || a_weight || jk || Ftest 
    return (specs, crosspecs, Tv)
  else
    if (outp == :cross) || (outp == :coh)
      S = Array{Array{Float64,1},2}(undef,p,p)
      for i = 1:p
        S[i,i]   = specs[i].S
        for j = (i+1):p
          S[i,j] = (outp == :coh) ? crosspecs[i,j].coh : crosspecs[i,j].S
          S[j,i] = crosspecs[i,j].phase 
        end
      end
    else
      S = mapreduce(x,hcat,(specs.S[i] for i = 1:p))
    end
    return S
  end
end


