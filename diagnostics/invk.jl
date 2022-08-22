## invk.jl -- checking the linearity of 1/k as a function of total cytochrome concentration

using CairoMakie, Statistics
import TOML, CytCKinetics

import Distributions

const volloc = joinpath("/Volumes/ESSENTIALS/CytCKinetics")
const dataloc = joinpath(volloc, "data")
const resultsloc = joinpath(dataloc, "results.toml")

const arglist = ["WTC:GroupA", "WTC:GroupC", "WTC:GroupD"]
const printname = ["May 16", "July 21", "July 25"]
parsearg(str) = split(str, ':') .|> string

function ssr(xs, ys, f)
	return sum(zip(xs, ys)) do (x, y)
		abs2(y - f(x))
	end
end

function r2fromvecs(xs, ys, f)
	mv = mean(ys)
	sst = sum(ys) do y
		abs2(mv - y)
	end
	ssr_ = ssr(xs, ys, f)

	return 1 - (ssr_/sst)
end
function quadfit(xs, ys)
	@assert length(xs) == length(ys)
	len = length(xs); T = eltype(xs)
	# build coef matrix
	coefmat = Matrix{T}(undef, len, 3)
	coefmat[:,1] .= one(T)
	@inbounds for i in (1,2)
        @. coefmat[:,i+1] = coefmat[:,i] * xs
	end
	return coefmat \ ys		# is there something faster than this?
end

# load the results file for each group
results = open(resultsloc) do io
    TOML.parse(io)
end
# calculate the 1/k for each point 

pairs = map(arglist) do arg
    strain, set = parsearg(arg)

    result = results[strain][set]
    casted = CytCKinetics.deserialize(CytCKinetics.SetGroupResults, result)

    concs = casted.concentrations
    invks = casted.concentrations ./ casted.initrates

    return concs, invks
end

@info "Proceeding to plotting routine"

f = Figure(textsize=25, font="CMU Sans Serif", resolution=(600,600))
ax = Axis(f[1,1],
    xlabel = "[Cyt c²⁺]₀ (μM)",
    ylabel = "Apparent 1/k (×10³ s)",
	title = "Linearity of 1/k against [Cyt c²⁺]₀"
)

for (idx, pair) in enumerate(pairs)
    xvals, yvals = pair
    scatter!(ax, xvals, yvals; color=Cycled(idx), markersize=5)

    ## simple linear (unweighted) regression
    m, b = CytCKinetics.weightedlinfit(xvals, yvals)
    ## weighted regression
    # m, b = CytCKinetics.weightedlinfit(xvals, yvals, 1 ./(yvals .^ 2))

	## quadratic fitting
	c, d, a = quadfit(xvals, yvals)

	r2 = r2fromvecs(xvals, yvals, x->m*x+b)
	r2_ = r2fromvecs(xvals, yvals, x->a*x^2+d*x+c)

	ablines!(ax, b, m; color=Cycled(idx), label=printname[idx])
	@info printname[idx] m b r2
	@info printname[idx] a d c r2_

	# calculate the pvalue
	ssr1 = ssr(xvals, yvals, x->m*x+b); dof1 = length(xvals) - 2
	ssr2 = ssr(xvals, yvals, x->a*x^2+d*x+c); dof2 = length(xvals) - 3

	score = ((ssr1 - ssr2)/(dof1 - dof2))/(ssr2 / dof2)
	@info "F-test pvalue" 1 - Distributions.cdf(Distributions.FDist(dof1-dof2, dof2), score)
end

leg = Legend(f[2,1], ax; 
			 orientation=:horizontal, tellheight=true, tellwidth=false, framevisible=false)

const plotdump = joinpath(volloc, "plotdump/diagnostics/invk")
const fname = "ACD.pdf"
const dump = joinpath(plotdump, fname)
save(dump, f)
@info "Saved image at $dump"

