using StatsFuns, Multitaper, Plots, RecipesBase

@userplot mtcoh

@recipe function f(h::MTCoherence; siglines = true, msclines = true, sigMax = 4, legtext = false, 
        force_xlims = nothing, force_ylims = nothing, mscaxis = true, sigaxis = true, jk = true,
        vlines = nothing, bwmark = nothing)
    
    # Layout, ylabels, etc
    if mscaxis && sigaxis
        layout := @layout [a{0.01w} b{0.95w} c{0.01w}]
        yguide --> [L"\hat{C}^2_{xy}(f)" L"z(f)" "Significance"]
        grid --> [false true false]
        legend --> [false legtext false]
        foreground_color_guide --> [:blue :black :red]
        xguide --> ["" "Frequency" ""]
        sp = [1,2,3]
    elseif jk && !mscaxis && sigaxis
        layout := @layout [b{0.95w} c{0.01w}]
        yguide --> [L"z(f)" "Significance"]
        grid --> [true false]
        legend --> [legtext false]
        xguide --> ["Frequency" ""]
        foreground_color_guide --> [:black :red]
        sp = [3, 1, 2]
    elseif mscaxis && !sigaxis
        layout := @layout [a{0.01w} b{0.95w}]
        yguide --> [L"\hat{C}^2_{xy}(f)" "Transformed MSC"]
        grid --> [false true]
        legend --> [false legtext]
        foreground_color_guide --> [:blue :black]
        xguide --> ["" "Frequency"]
        sp = [1,2,3]       
    elseif jk && !mscaxis && !sigaxis
        layout = (1,1)
        grid --> true
        legend --> legtext
        yguide --> L"z(f)"
        xguide --> "Frequency"
        sp = [1,1,1]
    elseif !jk && !mscaxis && !sigaxis
        layout = (1,1)
        grid --> true
        legend --> legtext
        yguide --> L"\hat{C}^2_{xy}(f)"
        xguide --> "Frequency"
        sp = [1,1,1]
    elseif !jk && !mscaxis && sigaxis
        layout := @layout [b{0.95w} c{0.01w}]
        yguide --> [L"\hat{C}^2_{xy}(f)" "Significance"]
        grid --> [true false]
        legend --> [legtext false]
        foreground_color_guide --> [:black :red]
        xguide --> ["Frequency" ""]
        sp = [1,1,2]
    end
    
    # Compute ticks
    z = norminvcdf(0,1,0.975)
    xl = (force_xlims == nothing) ? [0, h.f[end]] : force_xlims
    yl = (force_ylims == nothing) ? [0, maximum(h.coh .+ jk*z*sqrt.(h.jkvar[1]))] : force_ylims
    yl_transf = Multitaper.tanhtrans.(yl, h.params.K)
    ymsc = [0.4, 0.6, 0.8, 0.9, 0.99, 0.999]
    ymsc_transf = Multitaper.atanhtrans.(ymsc, h.params.K)
    sigs = Multitaper.logRange(1,sigMax,sigMax)
    labs = map(x -> Multitaper.atanhtrans(sqrt(Multitaper.invmscsig(x, h.params.K)), h.params.K), sigs) 

    if !jk && !mscaxis 
    yl = (force_ylims == nothing) ? Multitaper.tanhtrans.([0, maximum(h.coh)], h.params.K) : force_ylims
    labs = map(x -> Multitaper.invmscsig(x, h.params.K), sigs)
    @series begin
        subplot := sp[1]
        xlims --> xl
        ylims --> yl
        h.f[2:end], Multitaper.tanhtrans.(h.coh[2:end], h.params.K)
    end  
    # If needed to guide a visual comparison
    if siglines
      @series begin
        subplot := sp[1]
        ylims --> yl
        seriescolor --> :red
        line --> :dash
        label --> "" # sigs'
        h.f[[1,end]], vcat(labs',labs')
      end
    end
    # Significance if true coherence were zero
    if sigaxis
    @series begin
        subplot := sp[3]
        linealpha := 0.0
        ymirror := true
        yticks := true
        ylims := yl
        xticks := []
        yticks := (labs, sigs)
        tick_direction := :out
        xshowaxis := false
        h.f[2:end], h.coh[2:end]
    end
    end
    else
    # Un-transformed scale
    if mscaxis || !jk
    @series begin
        subplot := sp[1]
        linealpha := 0.0
        yticks --> (ymsc_transf, ymsc)
        xticks --> []
        ylims --> yl
        linewidth --> 2
        xshowaxis := false
        h.f[2:end], h.coh[2:end]
    end 
    end

    # Transformed scale
    if jk
    @series begin
      xlims := xl
      subplot := sp[2]
      seriescolor --> :mediumblue
      fillalpha --> 0.25
      fill := 1
      linealpha --> 0.25
      vcat(h.f[2:end],h.f[end:-1:2]), vcat((h.coh .+ z*sqrt.(h.jkvar[1]))[2:end], (h.coh .- z*sqrt.(h.jkvar[1]))[end:-1:2])
    end
    end
    if vlines != nothing
      @series begin
        subplot := sp[2]
        seriescolor --> :black
        line --> :dot
        linewidth --> 2
        vcat(vlines', vlines'), yl
      end
    end
    if bwmark != nothing
      @series begin
        subplot := sp[2]
        seriestype  :=  :path
        markershape := :vline
        seriescolor --> :black
        linealpha --> 0.6
        linewidth --> 2
        [bwmark[1]-h.params.NW/h.params.N/h.params.dt, bwmark[1]+h.params.NW/h.params.N/h.params.dt], [bwmark[2], bwmark[2]]
      end
    end
    @series begin
      subplot := sp[2]
      xlims := xl
      ylims := yl
      seriescolor --> :mediumblue
      linewidth --> 2
      if jk
        primary --> false
      else
        primary --> true
      end
      h.f[2:end], h.coh[2:end]
    end
    # If needed to guide a visual comparison
    if siglines
      @series begin
        subplot := sp[2]
        seriescolor --> :red
        line --> :dash
        label --> "" # sigs'
        xlims --> xl
        ylims --> yl
        h.f[[2,end]], vcat(labs',labs')
      end
    end
    if msclines
      @series begin
        subplot := sp[2]
        seriescolor --> :blue
        line --> :dash
        label --> "" # ymsc'
        xlims --> xl
        ylims --> yl
        h.f[[2,end]], vcat(ymsc_transf',ymsc_transf')
      end
    end
    # Significance if true coherence were zero
    if sigaxis
    @series begin
        subplot := sp[3]
        linealpha := 0.0
        xticks := []
        ymirror := true
        yticks := true
        ylims := yl
        yticks := (labs, sigs)
        tick_direction := :out
        xshowaxis := false
        h.f[2:end], h.coh[2:end]
    end
    end
    end
end
