
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

# Recipes for easy plots of the output variables made with the Multitaper package

##### Plots of multitaper spectra #####

@recipe function plot(S::mtspec; cross = true,
                      phase = false) 
  #
  yscale --> :log10
  dt = S.params.dt
  NW = S.params.NW
  K = S.params.K
  N = S.params.N
  xlab = (S.phase == nothing) ? "Frequency" : " "
  xlabel --> xlab
  ylab = (S.phase == nothing) ? "Spectrum" : "Cross-Spectrum"
  if (S.params.nsegments > 1); ylab *= "Welsh "; end
  ylabel --> ylab
  label --> ylab
  if S.jkvar != nothing
    @series begin
      label --> "Spectrum & 95% CI"
      fill := 1
      fillalpha --> 0.25
      linealpha --> 0.25
      z = norminvcdf(0,1,0.975)
      vcat(S.f[2:end],S.f[end:-1:2]), vcat(S.S[2:end].*exp.(z*sqrt.(S.jkvar[2:end])), S.S[end:-1:2].*exp.(-z*sqrt.(S.jkvar[end:-1:2])))
    end
  end
  @series begin
    primary --> false
    S.f[2:end], S.S[2:end]
  end
  if cross
    ejn = EJN((2*K)*S.params.nsegments)
    cent = [0.10*S.f[end],maximum(S.S)*0.8]
    @series begin 
      primary --> false
      cent[1] .+ BW(NW, N, dt)*[-1.0, 1.0], cent[2] .* [1.0, 1.0]
    end
    @series begin
      primary --> false
      cent[1] .* [1.0,1.0] , cent[2].*exp.(ejn*[-1.0, 1.0])
    end
  end
end

@recipe function plot(S::Vector{mtspec}; cross = true,
                      phase = false, xscaling = :lin) 
  #
  yscale --> :log10
  xlab = (S[1].phase == nothing) ? "Frequency" : " "
  xlabel --> xlab
  ylab = (S[1].phase == nothing) ? "Spectrum" : "Cross-Spectrum"
  ylabel --> ylab
  J = length(S)
  #
  for j = 1:J  
    if S[j].jkvar != nothing
      @series begin
        label --> "Spectrum & 95% CI"
        fill := 1
        fillalpha --> 0.25
        linealpha --> 0.25
        z = norminvcdf(0,1,0.975)
        vcat(S[j].f[2:end],S[j].f[end:-1:2]), vcat(S[j].S[2:end].*exp.(z*sqrt.(S[j].jkvar[2:end])), S[j].S[end:-1:2].*exp.(-z*sqrt.(S[j].jkvar[end:-1:2])))
      end
    end
    @series begin
      primary --> false
      S[j].f[2:end], S[j].S[2:end]
    end  
    # If the cross is desired
    if cross
      dt = S[j].params.dt
      K = S[j].params.K
      NW = S[j].params.NW
      N = S[j].params.N
      ejn = EJN(2*K*S[j].params.nsegments)
      cent = [0.15*S[j].f[end],maximum(S[j].S)*0.8]
      #
      @series begin 
        primary --> false
        cent[1] .+ BW(NW, N, dt)*[-1.0, 1.0], cent[2] .* [1.0, 1.0]
      end
      #
      @series begin
        primary -->  false
        cent[1] .* [1.0,1.0] , cent[2].*exp.(ejn*[-1.0, 1.0])
      end
    end
  end
end

### Plots of multitaper autocovariances

""" Simply plot the multitaper autocorrelation """
@recipe function mtacfplt(A::Multitaper.mtacf)
  label --> "MT acf"
  title --> "Autocorrelation Function"
  ylabel --> "Autocorrelation"
  xlabel --> "Lag"
  xlims --> [0.0, A.lags[end]]
  l --> :stem
  @series begin
    A.lags, A.acf
  end
end

""" Simply plot the multitaper autocovariance function """
@recipe function mtacvfplt(A::Multitaper.mtacvf)
  label --> "MT acvf"
  title --> "Autocovariance Function"
  ylabel --> "Autocovariance"
  xlabel --> "Lag"
  xlims --> [0.0, A.lags[end]]
  l --> :stem
  @series begin
    A.lags, A.acvf
  end
end

""" Simply plot the multitaper cepstrum coefficients """
@recipe function mtcepsplt(A::Multitaper.mtceps)
  label --> "MT Cepstrum"
  title --> "Cepstrum"
  ylabel --> "Cepstrum coefficient"
  xlabel --> "Lag"
  xlims --> [0.0, A.lags[end]]
  l --> :stem
  @series begin
    A.lags, A.ceps
  end
end

