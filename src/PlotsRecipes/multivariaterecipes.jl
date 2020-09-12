
# GNU GPL v2 licenced to C. Haley and C. Geoga 12/2019

# Recipes for multivariate spectrum analysis, mostly coherences at this point.

@recipe function plot(C::MtCoh; phase = false, sigMax = 0)
  jk = (C.jkvar != nothing) 
  link := :xaxis 
  NW = C.params.NW
  K = C.params.K
  dt = (1/(2*C.f[end]))
  z = norminvcdf(0,1,0.975)
  xguide --> "Frequency"
  if jk
    @series begin
      subplot := 1
      fillalpha --> 0.25
      label --> "MSC & 95% CI"
      fill := 1
      linealpha --> 0.25
      label := "Sq. Coherence"
      vcat(C.f,C.f[end:-1:1]), vcat(C.coh .+ z*sqrt.(C.jkvar[1]), (C.coh .- z*sqrt.(C.jkvar[1]))[end:-1:1])
    end
    @series begin
      yguide --> "Transformed Squared Coherence"
      subplot := 1
      primary --> false
      C.f, C.coh
    end
  else
    @series begin
      yguide --> "Squared Coherence"
      subplot := 1
      label --> "MSC"
      C.f, tanhtrans.(C.coh, K)
    end
  end
  # Label with significance levels
  sigs = logRange(1,sigMax,sigMax)
  labs = jk ? map(x -> atanhtrans(sqrt(invmscsig(x, K)), K), sigs) : map(x -> invmscsig(x,K),sigs)
  #
  @series begin 
    subplot := 1
    color := :red
    label --> 100*sigs'
    [C.f[1], C.f[end]], ones(2,1)*labs'
  end
  #
  if phase
    layout := (2,1)
    if C.jkvar != nothing
      @series begin
        subplot := 2
        fill := 1
        fillalpha --> 0.25
        yguide --> "Phase (deg)" 
        linealpha --> 0.25
        label --> "Phase"
        vcat(C.f,C.f[end:-1:1]), vcat(C.phase .+ z*C.jkvar[2], (C.phase .- z*C.jkvar[2])[end:-1:1])
      end
    end
    #
    @series begin
      xguide --> "Frequency"
      yguide --> "Phase (deg)"  
      primary --> false
      subplot := 2
      C.f, C.phase
    end
  end
end

@recipe function plot(C::MtTransf; phase = false)
  link := :xaxis 
  NW = C.params.NW
  K = C.params.K
  dt = C.params.dt
  @series begin
    yguide --> "Transfer function"
    subplot := 1
    yscale := :log10
    label --> "MSC"
    xlab = phase ? " " : "Frequency"
    xguide --> xlab
    C.f, C.transf
  end
  #
  if phase
    layout := (2,1)
    #
    @series begin
      xguide --> "Frequency"
      yguide --> "Phase (deg)"  
      primary --> false
      subplot := 2
      C.f, C.phase
    end
  end
end

@recipe function plot(C::Matrix{MtCoh}; phase = false, jk = false, sigMax = 0)
  link := :xaxis 
  z = norminvcdf(0,1,0.975)    
  # Label with significance levels
  sigs = logRange(1,sigMax,sigMax)
  labs = jk ? map(x -> atanhtrans(sqrt(invmscsig(x, C[1,2].params.K)), C[1,2].params.K), sigs) : 
              map(x -> invmscsig(x,C[1,2].params.K),sigs)
  if (C[1,2].jkvar != nothing) && jk
    for j in CartesianIndices(C)
      if (j[1] < j[2])*(j[1] != j[2])
        @series begin
          subplot := 1
          fillalpha --> 0.25
          fill := 1
          linealpha --> 0.25
          label := "Coherence & 95% CI"
          vcat(C[j].f,C[j].f[end:-1:1]), vcat(C[j].coh .+ z*sqrt.(C[j].jkvar[1]), (C[j].coh .- z*sqrt.(C[j].jkvar[1]))[end:-1:1])
        end
        # 
        @series begin
          yguide --> "Transformed Squared Coherence"
          subplot := 1
          #  seriescolor --> colour
          label := "Sq. Coherence"
          C[j].f, C[j].coh
        end
        @series begin 
          subplot := 1
          seriescolor := :red
          label --> 100*sigs'
          [C[j].f[1], C[j].f[end]], ones(2,1)*labs'
        end
      end
    end
  else
    for j in CartesianIndices(C)
      if (j[1] < j[2])*(j[1] != j[2])
        @series begin
          yguide --> "Squared Coherence"
          subplot := 1
          label --> "Sq. Coherence"
          C[j].f, tanhtrans.(C[j].coh, C[j].params.K)
        end    
        @series begin 
          subplot := 1
          seriescolor := :red
          label --> 100*sigs'
          [C[j].f[1], C[j].f[end]], ones(2,1)*labs'
        end
      end
    end
  end
  #
  if phase
    layout := (2,1)
    if jk
      for j in CartesianIndices(C)
        if (j[1] < j[2])*(j[1] != j[2])
          @series begin
            label --> "Phase (deg) & 95% CI"
            xguide --> "Frequency"
            subplot := 2
            fill := 1
            fillalpha --> 0.25
            linealpha --> 0.25
            vcat(C[j].f,C[j].f[end:-1:1]), vcat(C[j].phase .+ z*C[j].jkvar[2], (C[j].phase .- z*C[j].jkvar[2])[end:-1:1])
          end
        end
      end
    end
    #
    for j in CartesianIndices(C)
      if (j[1]<j[2])*(j[1] != j[2])
        @series begin
          yguide --> "Phase (deg)"
          label --> "Phase"
          xguide --> "Frequency"
          #  seriescolor --> colour
          subplot := 2
          C[j].f, C[j].phase
        end
      end
    end
  end
end

# Finally, suppose we want to plot the output of a coherence calculation. Currently the output is
# spectra in the first index and coherences in the second. This is not fancy, and does not provide
# CIs, you have to plot those separately to see detail. 
@recipe function specCohMatrix(Spec::Tuple{Array{MtSpec{C,T,P},1},Array{MtCoh,2},Nothing}; 
                                sigMax = 0) where C where T where P
  J = size(Spec[2],1)
  layout := (J,J)
  # Label with significance levels
  sigs = logRange(1,sigMax,sigMax)
  labs = map(x -> invmscsig(x,Spec[1][1].params.K),sigs)
  ser = 1
  for j in CartesianIndices(Spec[2])
    if j[2] > j[1]
      # Phase case - subdiagonal
      @series begin
        seriescolor --> 1
        label --> "Phase ($(j[1]),$(j[2]))"
        subplot := ser
        ser += 1
        Spec[2][j].f, Spec[2][j].phase
      end
    elseif j[2] == j[1]
      # Spectrum case - diagonal
      @series begin
        yscale --> :log10
        seriescolor --> 2
        label --> "Spectrum $(j[1])"
        subplot := ser
        ser += 1
        Spec[1][j[1]].f[2:end], Spec[1][j[1]].S[2:end]
      end
    else
      # MSC case - superdiagonal
      @series begin
        seriescolor --> :black
        ylims --> [0,1.0]
        label --> "MSC ($(j[1]),$(j[2]))"
        subplot := ser
        Spec[2][j[2],j[1]].f, tanhtrans.(Spec[2][j[2],j[1]].coh,Spec[2][j[2],j[1]].params.K)
      end
      @series begin 
        seriescolor := 2
        label --> 100*sigs'
        subplot := ser
        ser += 1
        [Spec[2][j[2],j[1]].f[1], Spec[2][j[2],j[1]].f[end]], ones(2,1)*labs'
      end
    end
  end
end

#### Recipes related to cross-covariance and cross-correlation

""" Plots the multitaper cross covariance function """
@recipe function ccvfplot(A::MtCcvf)
  label --> "MT Cross-Covariance"
  title --> "Cross-Covariance Function"
  yguide --> "Cross Covariance"
  xguide --> "Lag"
  l := :stem
  @series begin
    A.lags, A.ccvf
  end
end

""" Plots the multitaper cross correlation function """
@recipe function ccfplot(A::MtCcf)
  label --> "MT Cross-Correlation"
  title --> "Cross-Correlation Function"
  yguide --> "Cross Correlation"
  xguide --> "Lag"
  l := :stem
  @series begin
    A.lags, A.ccf
  end
end

