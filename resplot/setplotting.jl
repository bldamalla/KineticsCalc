## scripts/setplotting.jl --- for plotting stuff relating to set calculations
## groups are handled quite differently in another file (much later, sorry)

using CairoMakie
import CytCKinetics, Distributions
import TOML

# read the results file

# get the sets you are interested in
arglist = ["WTC:Set7", "WTC:Set8"]
split_arg(arg) = split(arg, ':') .|> string

# const resultsmetaloc = joinpath(@__DIR__, "..", "data/results.toml")
const volloc = joinpath("/Volumes/ESSENTIALS/", "CytCKinetics")
const resultsmetaloc = joinpath(volloc, "data/results.toml")

resultsmeta = open(resultsmetaloc) do io
    TOML.parse(io)
end

# get the seriessetresults
results = map(arglist) do arg
    strain, set = split_arg(arg)
    return resultsmeta[strain][set] |> q->CytCKinetics.deserialize(CytCKinetics.SeriesSetResults, q)
end
@info "Obtained results... Proceeding to plot routine"

# plot all of the objects in the arglist i guess

f = Figure()
ax = Axis(f[1,1],
    title="Comparison of activities in different measurement sets",
    xlabel="Initial [Cyt c²⁺]₀ (μM)",
    ylabel="Initial reaction rate (nM/s)"
)

paired = zip(arglist, results)
N = length(arglist)

for (i, (name, result)) in paired |> enumerate
    (; concentrations, initrates) = result
    gcolor = get(Makie.ColorSchemes.viridis, (i)/(N+1))
    scatter!(ax, result.concentrations, result.initrates; label=name, color=gcolor)

    # do that thing where !isnothing(results.fitparams)
    if !isnothing(result.fitparams)
        fparams = result.fitparams
        xvalues = 0:0.01:7.5
        yvalues = CytCKinetics.menten(xvalues, fparams)
        lines!(ax, xvalues, yvalues; color=gcolor)

        ## add confidence bands to see if the respective curves are within others' bands
        # declaration of distribution parameters
        α = 0.05; dofN = length(result.concentrations) - 2
        q = Distributions.quantile(Distributions.TDist(dofN), 1-α/2)
        σ_fse = map(xvalues) do x CytCKinetics.fitstderror(x, result) * q end

        # band!(ax, xvalues, yvalues .- σ_fse, yvalues .+ σ_fse; color=(gcolor, 0.5))
    end
end

leg = Legend(f[1,2], ax; framevisible=false)

# save the plot as pdf
const plotdump = joinpath(volloc, "plotdump/sets", "thing.pdf")
save(plotdump, f)
@info "Saved image at $plotdump"
