
using StatsFuns, Multitaper, Plots, RecipesBase

@userplot mtcoh

@recipe function f(h::MTCoherence; siglines = true, msclines = true, sigMax = 4, legtext = false, 
        force_ylims = nothing, mscaxis = true, sigaxis = true, jk = true)
    
    # Layout, ylabels, etc
    if mscaxis && sigaxis
        layout := @layout [a{0.01w} b{0.95w} c{0.01w}]
        yguide --> ["MSC" "Transformed MSC" "Significance"]
        grid --> [false true false]
        legend --> [false legtext false]
        foreground_color_guide --> [:blue :black :red]
        xguide --> ["" "Frequency" ""]
        sp = [1,2,3]
    elseif jk && !mscaxis && sigaxis
        layout := @layout [b{0.95w} c{0.01w}]
        yguide --> ["Transformed MSC" "Significance"]
        grid --> [true false]
        legend --> [legtext false]
        xguide --> ["Frequency" ""]
        foreground_color_guide --> [:black :red]
        sp = [3, 1, 2]
    elseif mscaxis && !sigaxis
        layout := @layout [a{0.01w} b{0.95w}]
        yguide --> ["MSC" "Transformed MSC"]
        grid --> [false true]
        legend --> [false legtext]
        foreground_color_guide --> [:blue :black]
        xguide --> ["" "Frequency"]
        sp = [1,2,3]       
    elseif jk && !mscaxis && !sigaxis
        layout = (1,1)
        grid --> true
        legend --> legtext
        yguide --> "Transformed MSC"
        xguide --> "Frequency"
        sp = [1,1,1]
    elseif !jk && !mscaxis && !sigaxis
        layout = (1,1)
        grid --> true
        legend --> legtext
        yguide --> "MSC"
        xguide --> "Frequency"
        sp = [1,1,1]
    elseif !jk && !mscaxis && sigaxis
        layout := @layout [b{0.95w} c{0.01w}]
        yguide --> ["MSC" "Significance"]
        grid --> [true false]
        legend --> [legtext false]
        foreground_color_guide --> [:black :red]
        xguide --> ["Frequency" ""]
        sp = [1,1,2]
    end
    
    # Compute ticks
    z = norminvcdf(0,1,0.975)
    yl = (force_ylims == nothing) ? [0, maximum(h.coh .+ jk*z*sqrt.(h.jkvar[1]))] : force_ylims
    yl_transf = Multitaper.tanhtrans.(yl, h.params.K)
    ymsc = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.99, 0.999]
    ymsc_transf = Multitaper.atanhtrans.(ymsc, h.params.K)
    sigs = Multitaper.logRange(1,sigMax,sigMax)
    labs = map(x -> Multitaper.atanhtrans(sqrt(Multitaper.invmscsig(x, h.params.K)), h.params.K), sigs) 

    if !jk && !mscaxis 
    yl = (force_ylims == nothing) ? Multitaper.tanhtrans.([0, maximum(h.coh)], h.params.K) : force_ylims
    labs = map(x -> Multitaper.invmscsig(x, h.params.K), sigs)
    @series begin
        subplot := sp[1]
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
        ylims --> yl
        xshowaxis := false
        h.f[2:end], h.coh[2:end]
    end 
    end

    # Transformed scale
    if jk
    @series begin
      subplot := sp[2]
      fillalpha --> 0.25
      label --> "MSC & 95% CI"
      fill := 1
      linealpha --> 0.25
                vcat(h.f[2:end],h.f[end:-1:2]), vcat((h.coh .+ z*sqrt.(h.jkvar[1]))[2:end], (h.coh .- z*sqrt.(h.jkvar[1]))[end:-1:2])
    end
    end
    @series begin
      subplot := sp[2]
      ylims := yl
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
        ylims --> yl
        h.f[[2,end]], vcat(ymsc_transf',ymsc_transf')
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
        yticks := (labs, sigs)
        tick_direction := :out
        xshowaxis := false
        h.f[2:end], h.coh[2:end]
    end
    end
    end
end
