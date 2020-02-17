
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

""" Eigenvalues/concentrations for the  Slepian sequences.   
They are computed as given in Percival 390, BUT THERE IS A TYPO IN PERCIVAL 390 """
function dpss_eigval(dpVecs::Array{Float64, 2}, n::Int, nw::Real, ntapers::Int)::Vector{Float64}
  sincvec  = [sin(2*pi*nw*j/n)/(pi*j) for j in 1:(n-1)]
  q_tau    = mapreduce(x -> conv(x, x)[n:(2*n-1)], hcat, (dpVecs[:,j] for j = 1:ntapers))
  eigvalss = map(j->2*(nw/n + dot(q_tau[1,j].*q_tau[2:n,j], sincvec)), 1:ntapers)
  return eigvalss
end


""" Simply compute dpss tapers, eigenvalues """
function dpss_tapers(n::Int64, nw::Float64, k::Int64, tap_or_egval::Symbol=:tap) 
  stdm  = SymTridiagonal([cos(2*pi*(nw/n))*abs2(0.5*(n-1)-(j-1)) for j in 1:n],
                          [0.5*j*(n-j) for j in 1:(n-1)])
  # Note that eigvals seems to sort ascending.
  vv = reverse(eigvecs(stdm, reverse(eigvals(stdm, (n-k+1):n))),dims=2)
  if tap_or_egval == :tap
    return vv  
  elseif tap_or_egval == :egval
    return dpss_eigval(vv, n, nw, k)
  elseif tap_or_egval == :both
    return(vv, dpss_eigval(vv, n, nw, k))
  end
end

# a dot function that is faster for very small vectors. for aweight function.
function _dot(v1, v2)
  s = 0.0
  @inbounds @simd for j in eachindex(v1, v2)
    s += v1[j]*v2[j]
  end
  s
end

# a specific weight updating function. for aweight function.
function weightupdate!(new, est, evalues, evar)
  @inbounds @simd for j in eachindex(new, evalues)
    new[j] = sqrt(evalues[j])*est/(evalues[j]*est+evar*(1.0-evalues[j]))^2
  end
  nothing
end

# obtained adaptively weighted estimates. somewhat micro-optimized.
function aweighted(old, new, len, estimates, evalues, evar, maxit, tol)
  est = _dot(estimates, old)/sum(old)
  for _ in 1:maxit
    weightupdate!(new, est, evalues, evar)
    new_est = _dot(estimates, new)/sum(new)
    isapprox(est, new_est, rtol=tol) && return new_est
    est  = new_est
    old .= new
  end
  est
end

""" Guts function for jackknifing spectra """
function jknife_spec(tap::ecoef, tap2::Union{ecoef,Nothing} = nothing,
                     loo_or_psv::Symbol = :loo) 
  # Jackknife the spectra
  if tap2 == nothing
    pseudovals = abs2.(tap.coef)
    # Jackknife the cross-spectra, phase estimates might be a bit wonky
  else
    pseudovals = tap.coef.*conj(tap2.coef)
  end
  allin = sum(pseudovals,dims=2)[:,1]
  pseudovals .-= repeat(allin,1,size(tap.coef,2))
  pseudovals ./= -(size(tap.coef,2)-1)
  map!(log, pseudovals, abs.(pseudovals))
  map!(log, allin, allin./size(tap.coef,2))
  #Scapdot = mean(pseudovals,dims=2)[:,1]
  if loo_or_psv == :loo
    return pseudovals
  else
    return (size(tap.coef,2)*repeat(allin,1,size(tap.coef,2))-(size(tap.coef,2)-1)*(pseudovals))
  end
end

""" Guts function for jackknifing coherences, leave one out values are returned"""
function jknife_coh(tap::ecoef, tap2::ecoef, loo_or_psv::Symbol = :loo)
  loo_coh   = exp.(jknife_spec(tap, tap2, :loo))
  loo_coh ./= sqrt.(exp.(jknife_spec(tap, nothing, :loo)).*exp.(jknife_spec(tap2,nothing,:loo)))
  map!(x->atanhtrans(x,size(tap.coef,2)), loo_coh, loo_coh) 
  if loo_or_psv == :loo
    return loo_coh
  else
    allin = sum(tap.coef.*tap2.coef, dims = 2)[:,1]
    allin ./= sqrt.(sum(abs2.(tap.coef),dims=2)[:,1].*sum(abs2.(tap2.coef),dims=2)[:,1])
    map!(x->atanhtrans(x,size(tap.coef,2)), allin, allin) 
    return (size(tap.coef,2)*repeat(allin,1,size(tap.coef,2)) - (size(tap.coef,2)-1)*loo_coh)
  end
end

""" Guts function for jackknifing the phase, where eigencoefficients are given """
function jknife_phase(tap::ecoef, tap2::ecoef)
  loos    = tap.coef.*conj(tap2.coef)
  allin   = sum(loos,dims=2)[:,1]
  loos  .-= repeat(allin,1,size(tap.coef,2))
  loos  ./= -abs.(loos)
  # From TC91, use Eqn (2.61)
  phase   = mean(loos,dims=2)[:,1]
  # This function gives jackknifed standard deviations for the phase in degrees
  return unwrapphase(angle.(phase)*180/pi,:deg), 
  (180/π)*sqrt.(2.0*(size(tap.coef,2) - 1.0)*(1.0 .- abs.(phase)))
end

""" Guts function for jackknifing the phase, where leave one out coherences are given """
function jknife_phase(loo_coh::Union{Matrix{Float64}, Matrix{ComplexF64}})
  ej     = loo_coh./abs.(loo_coh)
  phase  = mean(ej,dims=2)[:,1]
  return unwrapphase(angle.(phase)*180/pi,:deg),
  (180/π)*sqrt.(angle.(2*(size(tap.coef,2)-1)*(1-abs.(phase))))
end

""" Here's the tdist quantile relevant to jackknife CI's """
Tquant(K,alpha::Float64=0.05) = StatsFuns.tdistinvcdf(K-1, 1.0-alpha)

""" Function that jackknifes spectra and/or coherences """
function jknife(tap::ecoef, tap2::Union{ecoef,Nothing}=nothing, typ::Symbol = :spec)
  K      = size(tap.coef,2)
  # pseudovalues come pre-transformed
  loos   = (typ == :coh) ? jknife_coh(tap, tap2, :loo) : jknife_spec(tap,tap2,:loo)
  allin  = mean(loos, dims=2)[:,1]
  jkvar  = var(loos,dims=2,corrected=false)[:,1]*(K-1)
  jkmean = K*allin .- (K-1)*mean(loos,dims=2)[:,1]
  # Returns squared coherence
  fun(x) = (typ == :spec) ? exp(abs(x)) : abs(x)
  return fun.(jkmean), jkvar
end

""" Function that jackknifes transfer functions from the stacked eigencoefficients of a multivariate time series.
Note that the result is a p by K matrix and represents the coefficients of a regression of one series on the p-1 remaining series.  """
function mttransfn(Y::Vector{ComplexF64}, Ec::Union{Vector{ComplexF64},Matrix{ComplexF64}})
  Ec  = hcat(Y,Ec)
  K,p = size(Ec)
  if K<2
    error("Use more than 2 tapers to get your transfer function.")
  end
  # the right singular vector is a K x p matrix, (K-1) x p in the case of Leave-one-out, so the
  # right singular vector is p x p
  allin  = svd(Ec).V
  # Eigenvector corresponding to the smallest eigenvalue of the matrix
  allin  = -allin[1,p]./allin[2:p,p]
  loos   = map(i -> svd(Ec[vcat(collect(1:(i-1)),collect((i+1):K)),:]).V, 1:K)
  loos   = mapreduce(x -> -x[1,p]./x[2:p,p], hcat, loos)
  jkvar  = var(loos, dims = 2, corrected = false)[:]*(K-1)
  jkmean = K*allin .- (K-1)*mean(loos, dims=2)[:]
  resid  = Ec[:,2:p]*jkmean .- Y
  return jkmean, jkvar, resid
end

""" F-test for line components at every frequency. Can select raw statistic or test p-value. """
function Ftest_pval(dc::Vector{Float64}, tap::Union{Matrix{ComplexF64},ecoef}, 
                    stat_or_p::Symbol = :pval) 
  tap    = (typeof(tap) == ecoef) ? tap.coef : tap
  k      = Float64(size(tap, 2))
  sum2dc = sum(abs2, dc)
  uf     = tap*dc/sum2dc
  ff     = (k-1)*abs2.(uf)*sum2dc
  ff   ./= sum(abs2, tap .- uf.*transpose(dc), dims=2)[:,1]
  if stat_or_p != :stat
    map!(x -> 1.0 - cdf( FDist(2, 2*Int64(k) - 2), x), ff, ff)
  end
  return (ff, uf)
end

""" Compute the output length of the spectrum """
function output_len(S1::Vector{T}, pad::Union{Int,Float64} = 1.0) where T<:Number
  lengt         = length(S1)
  if typeof(pad) == Float64
    fftleng     = (pad >= 1) ? Int64(round(pad*lengt)) : error("Must pad by a factor >= 1.0")
  elseif typeof(pad) == Int64
    fftleng     = (pad >= length(S1)) ? pad : error("Must pad to a number greater than the length of the input vector.")
  end
  # If the series is complex, return all frequencies, if real return only half
  if (typeof(S1)==Vector{ComplexF64})
    halffreq    =  fftleng 
  else 
    halffreq    = isodd(fftleng) ? Int64(ceil((fftleng)/2)) : Int64((fftleng)/2)+1
  end
  return lengt, fftleng, halffreq
end

""" Computes multitaper guts, both eigencoefficients and adaptive weighting """
function mtspec_guts(S1::Vector{T}, 
                     dpVec::Matrix{Float64}, 
                     fftleng::Int64,
                     halffreq::Union{Int64,Nothing},
                     cvar::Float64 = 1.0,
                     ctr::Bool = true, 
                     egval::Union{Nothing,Matrix{Float64}} = nothing,
                     maxit = 15,
                     tol = 0.05
                    ) where T<:Number
  # Compute the eigencoefficients 
  S1      .-= mean(S1).*(ctr==true)
  k         = size(dpVec,2)
  eigcoefs  = mapreduce(slep -> fft(vcat(S1.*slep,zeros(fftleng-length(S1)))), hcat, 
                  (dpVec[:,j] for j in 1:k))
  dsq = (egval == nothing) ? nothing : zeros(halffreq,k)
  if egval != nothing
    (old, new) = vcat(ones(2).*0.5, zeros(k-2)), zeros(k) # is this weird with threads?
    Threads.@threads for j in 1:halffreq
      @inbounds dsq[j,:] = aweighted(old, new, k, abs2.(view(eigcoefs,j,:)), 
                                   egval, cvar, maxit, tol)
    end
  end
  return ecoef(eigcoefs[1:halffreq,:], dsq)
end

""" Computes univariate multitaper spectra with a handful of extra gadgets. """
function multispec(S1::Union{Vector{T}}; NW::Real = 4.0, K::Int = 6, 
                   dt::Float64=1.0, ctr::Bool = true, 
                   pad::Union{Int,Float64} = 1.0, 
                   dpVec::Union{Vector{Float64},Matrix{Float64},Nothing} = nothing,
                   egval::Union{Vector{Float64},Nothing} = nothing,
                   guts::Bool = false, a_weight::Bool = true, 
                   Ftest::Bool = false, highres::Bool = false,
                   reshape::Union{Bool, Float64, Int64, Vector{Float64}, Vector{Int64}} = false,
                   jk::Bool = false, 
                   Tsq::Union{Vector{Float64},Vector{Vector{Float64}},Vector{Int64},Vector{Vector{Int64}},Nothing}=nothing, 
                   alph::Float64 = 0.05) where T<:Number
  # Compute the lengths of the spec obj
  lengt, fftleng, halffreq = output_len(S1,pad)
  # Compute the array of dpss vectors if they aren't given.
  if (dpVec == nothing) || ((egval == nothing)&&(a_weight==true))
    dpVec = dpss_tapers(lengt, NW, K)
    dpVec .*= sqrt(dt)
    egval = a_weight ? dpss_tapers(lengt, NW, K, :vals) : nothing
  end
  cvar     = var(S1)
  coefswts = mtspec_guts(S1,dpVec,fftleng,halffreq,cvar,ctr,egval)

  # Compute the spectrum estimate
  S        = (coefswts.wts == nothing) ? mean(abs2.(coefswts.coef),dims=2)[:,1] : 
                                         sum(abs2.(coefswts.coef).*coefswts.wts, dims=2)[:,1]
  jv       = jk ? jknife(coefswts,nothing,:spec)[2] : nothing

  # Compute the F- and T^2- tests if requested:
  dcs      = (Ftest || Tsq != nothing || highres) ? map(isodd,1:K).*sum(dpVec,dims=1)[:] : nothing

  # Do the F-tests
  fv,uf    = (Ftest || highres) ? Ftest_pval(dcs, coefswts.coef[1:halffreq,:]) : (nothing, nothing)
  S        = highres ? abs2.(uf)/dt : S

  freq = (1/dt)*LinRange(0,1,fftleng+1)[1:halffreq]
  if (reshape != false)
    if (reshape == Float64)||(reshape == Vector{Float64})
      Sigtest = map(x -> freqfinder(x, freq), reshape)
    elseif (reshape == Int64)||(reshape == Vector{Int64})
      Sigtest = reshape
    else
      Sigtest = findall(fv .<= (1.0/length(S1)))
    end
    lastind     = isodd(fftleng) ? halffreq : halffreq-1
    oddev_shift = isodd(fftleng) ? 2 : 3
    Uk = mapreduce(x -> fftshift(fft(vcat(x, zeros(fftleng-length(S1))))), hcat, (dpVec[:,k] for k in 1:K))
    SL = mapreduce(st -> uf[st]*vcat(Uk[(end-st+oddev_shift):end,:],Uk[1:(end-st+oddev_shift-1),:]), +, Sigtest)
    S  = mean(abs2.(coefswts.coef .- SL[lastind:end,:]),dims=2)[:,1]
    for sk = Sigtest
      S[sk] += abs2(uf[sk])/dt
    end
  end

  # Do the Tsquared test
  if typeof(Tsq) != Nothing
    Tsq      = (typeof(Tsq) <: Vector{Number}) ? [Tsq] : Tsq
    map!(x -> freq_to_int(Tsq[x], lengt, dt), Tsq, eachindex(Tsq))
    Tsq = Vector{Vector{Int64}}(Tsq)  
    if (2*K < (true ? 1 : 2)*maximum(length.(Tsq)))
      error("There are too few tapers for the number of Tsq tests.")
    end
    Tv = map(x->Tsqtest_pval(dcs,ecoef(coefswts.coef[Tsq[x],:],nothing)),eachindex(Tsq)) 
  else
    Tv = nothing
  end
  params = mtparams(NW, K, lengt, dt, fftleng, 1, nothing) # NW, K, N, dt, M, nsegments, overlap
  phase = nothing
  return mtspec(freq, S, phase, params, (guts ? coefswts : nothing), fv, jv, Tv) 
  # frequency, spectrum or crosspec, phase, parameters, eigencoefficients & weights, Ftest, jacknife, Tsq test
end

""" Blocker code to divide the data into segments """
function blockerr(lengt::Int64, nsegments::Int64; overlap::Float64=0.0)
  if (overlap < 1.0)&&(overlap >= 0.0)
    seg_len = Int64(ceil((lengt-1)/((nsegments-1)*(1-overlap)+1)))
    step    = Int64(round((1-overlap)*seg_len))
    step    = maximum(vcat(step,1))
    seq     = (nsegments == 1) ? [1] : collect(0:step:lengt-1)[1:nsegments] .+ 1
    if nsegments == 1
      ov = 0.0
    elseif nsegments == 2
      seq[nsegments]  = (lengt-seg_len+1)
      ov = 1 - (seq[2]-seq[1])/seg_len
    else
      seq[nsegments]  = (lengt-seg_len+1)
      ov      = 1 - step/seg_len
    end
  else
    error("Overlap should be between 0 and 1.")
  end
  return seq, seg_len, ov
end

""" Computes univariate multitaper autocovariance/autocorrelation function. Inputs a mtspec struct. """
function mt_acvf(S::mtspec; typ::Symbol = :acvf)   
  lags = S.params.dt*S.params.N*LinRange(0.0, 1.0, S.params.N + 1)[1:length(S.S)]
  spec = mod(S.params.N, 2) == 0 ? vcat(S.S, S.S[end-1:-1:2]) : vcat(S.S, S.S[end:-1:2])
  if typ == :acvf
    return mtacvf(lags, real.(ifft(spec))[1:length(S.S)], S.params)
  elseif typ == :acf
    acvf = real.(ifft(spec))[1:length(S.S)]
    return mtacf(lags, acvf/acvf[1], S.params)
  elseif typ == :ceps
    return mtceps(lags, real.(ifft(log.(spec)))[1:length(S.S)], S.params)
  else
    error("Select one of :acvf (autocovariance), :acf (autocorrelation), :ceps (cepstrum) for output")
  end
end

""" Computes univariate multitaper autocovariance/autocorrelation function if you haven't already computed the MT spectrum. """
function mt_acvf(S1::Union{Vector{T}}; 
                 typ::Symbol = :acvf, 
                 NW::Real = 4.0, 
                 K::Int = 6, 
                 dt::Float64=1.0, 
                 ctr::Bool = true, 
                 pad::Union{Int,Float64} = 1.0, 
                 dpVec::Union{Vector{Float64},Matrix{Float64},Nothing} = nothing,
                 egval::Union{Vector{Float64},Nothing} = nothing,
                 a_weight::Bool = true, 
                 reshape::Union{Bool, Float64, Int64, Vector{Float64}, Vector{Int64}} = false,
                 ) where T<:Number
  S = multispec(S1, NW = NW, K = K, dt = dt, ctr = ctr, pad = pad, dpVec = dpVec, egval = egval, guts = false, 
                a_weight = a_weight, Ftest = false, reshape = reshape, jk = false, Tsq = nothing) 
  return mt_acvf(S, typ = typ)
end

""" Because of the mapreduce function, this tidy piece of code can be used to get either
(a) the MT Welsh spectrum, and (b) the spectrogram """
function welsh(S1::Union{Vector{Float64}, Vector{ComplexF64}, Matrix{Float64}, Matrix{ComplexF64}}, 
               nsegments::Int64, overlap::Float64 = 0.5, ws::Symbol = :welsh; outp::Symbol = :spec,
               NW::Real = 4.0, K::Int = 6, dt::Float64 = 1.0,
               ctr::Bool = true, pad::Union{Int,Float64} = 1.0,
               dpVec::Union{Vector{Float64},Matrix{Float64},Nothing} = nothing, 
               egval::Union{Vector{Float64},Nothing} = nothing, 
               guts::Bool = false, a_weight::Bool = true, Ftest::Bool = false, jk::Bool =
               false, Tsq::Union{Array{Int64,1},Array{Array{Int64,1},1},Nothing}=nothing,
               alph::Float64 = 0.05) 

  multiv    = (typeof(S1) <: Matrix)
  lengt,p   = multiv ? size(S1) : (length(S1),nothing) 
  if multiv && (p > 2)
    error("Can compute cross-spectrograms and cross-welsh spectra for no more than two series at a time.")
  end

  # Get the sizes of the data chunks, note that overlap gets overwritten
  seq,seg_len,overlap = blockerr(lengt,nsegments,overlap=overlap) 

  # Effective bandwidth (2.13 in T&C91)
  bw        = NW*(1 + (lengt - 1)*(mean(diff(seq))/seg_len))/(lengt*dt) 
  if (dpVec == nothing)||((egval == nothing)&&(a_weight == true))
    dpVec   = dpss_tapers(seg_len, NW, K)
    dpVec .*= sqrt(dt)
    egval   = a_weight ? dpss_tapers(seg_len, NW, K, :vals) : nothing
  end

  # Cross-spectra and coherences are a possibility for a future version of this code. That
  # would require the substitution of the jackknife estimate. 
  fun = (ws == :welsh) ? (+) : (hcat)
  if multiv
    if (outp != :spec)&&(outp != :coh)
      error("Output is either cross-spectrum or coherence")
    end
    if ws != :welsh
      # For the spectrogram, we call multispec and set the results on each segment side by side.
      if outp == :spec
        v = mapreduce(x -> multispec(S1[x:(x+seg_len-1),1],S1[x:(x+seg_len-1),2],outp = :spec,
                                     NW = NW, K = K, dt = dt, ctr = ctr, pad = pad,
                                     dpVec = dpVec, guts = false, Tsq = nothing, jk =
                                     false, alph = 0.05).S, hcat, seq) 
      elseif outp == :coh
        v = mapreduce(x -> multispec(S1[x:(x+seg_len-1),1],S1[x:(x+seg_len-1),2],outp = :coh,
                                     NW = NW, K = K, dt = dt, ctr = ctr, pad = pad,
                                     dpVec = dpVec, guts = false, Tsq = nothing, jk = false, alph = 0.05).coh, 
        hcat, seq) 
      end
    else
      # If the welsh spectrum is desired, get the eigencoefficients over each of the segments and
      # jackknife those to get a single estimate. 
      coefswts1 = mapreduce(x->mtspec_guts(S1[x:(x+seg_len-1),1],dpVec,fftleng,halffreq,var(x),ctr,nothing).coef,
                            hcat,seq)
      coefswts2 = mapreduce(x->mtspec_guts(S1[x:(x+seg_len-1),2],dpVec,fftleng,halffreq,var(x),ctr,nothing).coef,
                            hcat,seq)
      v = jknife(ecoef(coefswts1, nothing), ecoef(coefswts2, nothing), outp)[1]
      v .*= (outp == :cross) ? dt : 1.0
    end
  else
    v = mapreduce(x -> multispec(S1[x:(x+seg_len-1)],NW = NW, K = K, dt = dt, ctr = ctr, pad = pad,
                                 dpVec = dpVec, guts = false, Ftest = false, jk =
                                 false, alph = 0.05, egval = egval, a_weight = a_weight).S, fun, seq)
    fftleng, halffreq = output_len(S1[seq[1]:(seg_len+seq[1]-1)],pad)[2:3]
  end   
  params = mtparams(NW, K, lengt, dt, fftleng, nsegments, overlap) # NW, K, N, dt, M, nsegments, overlap
  freq = LinRange(0,1.0,fftleng+1)[1:halffreq]*(1.0/dt)

  if ws == :welsh
    return (mtspec(freq, v/nsegments, nothing, params, nothing, nothing, nothing), bw)
  elseif ws == :sgram
    return SpectroGram(dt*seq, freq, v, params, bw) 
  end

end

