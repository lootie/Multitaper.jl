
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

# Recipes for easy plots of the output variables made with the Multitaper package

##### Plots of multitaper spectra #####

@recipe function plot(S::MtSpec; cross = true,
                      phase = false) 
  #
  yscale --> :log10
  dt   = S.params.dt
  NW   = S.params.NW
  K    = S.params.K
  N    = S.params.N
  xlab = (S.phase == nothing) ? "Frequency" : " "
  xguide --> xlab
  ylab = (S.phase == nothing) ? "Spectrum" : "Cross-Spectrum"
  if (S.params.nsegments > 1); ylab *= "Welch "; end
  yguide --> ylab
  @series begin
    S.f[2:end], S.S[2:end]
  end 
  if S.jkvar != nothing
    @series begin
      primary --> false
      label --> "Spectrum & 95% CI"
      fill := 1
      fillalpha --> 0.25
      linealpha --> 0.25
      z = norminvcdf(0,1,0.975)
      vcat(S.f[2:end],S.f[end:-1:2]), vcat(S.S[2:end].*exp.(z*sqrt.(S.jkvar[2:end])), S.S[end:-1:2].*exp.(-z*sqrt.(S.jkvar[end:-1:2])))
    end
  end
  if cross
    EJN = ejn((2*K)*S.params.nsegments)
    cent = [0.10*S.f[end],maximum(S.S)*0.8]
    @series begin 
      primary --> false
      cent[1] .+ (NW/(N*dt))*[-1.0, 1.0], cent[2] .* [1.0, 1.0]
    end
    @series begin
      primary --> false
      cent[1] .* [1.0,1.0] , cent[2].*exp.(EJN*[-1.0, 1.0])
    end
  end
end

### Plots of multitaper autocovariances

""" Simply plot the multitaper autocorrelation """
@recipe function acfplt(A::MtAcf)
  label --> "MT acf"
  title --> "Autocorrelation Function"
  yguide --> "Autocorrelation"
  xguide --> "Lag"
  xlims --> [0.0, A.lags[end]]
  l --> :stem
  @series begin
    A.lags, A.acf
  end
end

""" Simply plot the multitaper autocovariance function """
@recipe function acvfplt(A::MtAcvf)
  label --> "MT acvf"
  title --> "Autocovariance Function"
  yguide --> "Autocovariance"
  xguide --> "Lag"
  xlims --> [0.0, A.lags[end]]
  l --> :stem
  @series begin
    A.lags, A.acvf
  end
end

""" Simply plot the multitaper cepstrum coefficients """
@recipe function cepsplt(A::MtCeps)
  label --> "MT Cepstrum"
  title --> "Cepstrum"
  yguide --> "Cepstrum coefficient"
  xguide --> "Lag"
  xlims --> [0.0, A.lags[end]]
  l --> :stem
  @series begin
    A.lags, A.ceps
  end
end

""" Demodulate recipe """
@recipe function demodplt(cdm::Demodulate)
  layout --> (2,1)
  @series begin
    subplot := 1
    yguide --> "Magnitude"
    label --> "Magnitude"
    cdm.time, cdm.mag
  end
  @series begin
    subplot := 2
    yguide --> "Phase"
    label --> "Phase"
    xguide --> "Time"
    cdm.time, cdm.phase
  end
end

