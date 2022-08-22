## scripts/groupplotting.jl --- plotting with group information
## statistical happiness and hardships

using CairoMakie
import CytCKinetics, Distributions
import TOML

# get the groups you are interested in
arglist = ["WTC:GroupC", "WTC:GroupD"]
split_arg(arg) = split(arg, ':') .|> string

# const resultsmetaloc = joinpath(@__DIR__, "..", "data", "results.toml")
const volloc = joinpath("/Volumes/ESSENTIALS", "CytCKinetics")
const dataloc = joinpath(volloc, "data")
const resultsmetaloc = joinpath(dataloc, "results.toml")

resultsmeta = open(resultsmetaloc) do io
    TOML.parse(io)
end

# get the groupsetresults
results = map(arglist) do arg
    strain, group = split_arg(arg)
    return resultsmeta[strain][group] |> q->CytCKinetics.deserialize(CytCKinetics.SetGroupResults, q)
end
@info "Obtained results... Proceeding to plotting routine"

# plot all of the objects in the arglist i guess

f = Figure()
ax = Axis(f[1,1],
    title = "Reaction rate curves for each group",
    xlabel = "Initial [Cyt c²⁺]₀ (μM)",
    ylabel = "Initial reaction rate (nM/s)"
)
ax2 = Axis(f[2,1],
    title = "95% confidence intervals for kinetic parameters",
    ylabel = "Km (μM)", xlabel = "kcat (s⁻¹)"
)

paired = zip(arglist, results)

for (i, (name, result)) in paired |> enumerate
    (; concentrations, initrates) = result
    scheme_cycle = [ 1//5, 3//5, 2//5, 4//5 ]
    gcolor = get(Makie.ColorSchemes.viridis, scheme_cycle[i])
    scatter!(ax, concentrations, initrates;
        label=name, color=gcolor, markersize=5)

    # declaration of distribution parameters
    α = 0.05; dofN = length(result.concentrations) - 2

    # the fitting curve and its confidence region
    let fparams = result.fitparams
        xvalues = 0:0.01:7.5
        yvalues = CytCKinetics.menten(xvalues, fparams)
        lines!(ax, xvalues, yvalues; color=gcolor)

        # the confidence band
        # declare a significance value of 0.05 for all experiments
        q = Distributions.quantile(Distributions.TDist(dofN), 1-α/2)
        
        σ_fse = map(xvalues) do x CytCKinetics.fitstderror(x, result) * q end
        #band!(ax, xvalues, yvalues .- σ_fse, yvalues .+ σ_fse;
        #    color=(gcolor, 0.5))
    end

    # the confidence regions
    let centers = result.fitparams, widths = 4 .* result.stderrors
        θ1range = LinRange(centers[1]-widths[1], centers[1]+widths[1], 150)
        θ2range = LinRange(centers[2]-widths[2], centers[2]+widths[2], 150)
        
        q = Distributions.quantile(Distributions.TDist(dofN), 1-α/2)
        scatter!(ax2, centers; color=gcolor, markersize=5, label=name)

	# error bars representing 95% confidence regions from each group
	errs = result.stderrors .* q
	x, y = centers
	errorbars!(ax2, [x], [y], errs[1]; direction=:x, whiskerwidth=10, color=gcolor)
	errorbars!(ax2, [x], [y], errs[2]; direction=:y, whiskerwidth=10, color=gcolor)
    end

end

# legends
leg1 = Legend(f[1,2], ax; framevisible=false)
leg2 = Legend(f[2,2], ax2; framevisible=false)

const plotdump = joinpath(volloc, "plotdump/groups", "GroupCD_.pdf")
save(plotdump, f)
@info "Saved image at $plotdump"
