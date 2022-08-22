## diagnostics/order.jl -- check the range for which assuming steady-state solution is valid

using CytCKinetics, Statistics
using CairoMakie

import TOML

# work on sets and not groups
const arglist = ["WTC:Set2", "WTC:Set6", "WTC:Set7"]
const printnames = ["May 16", "July 21", "July 25"]
const argnum = 2		# work on concentration = 0.5μM
split_arg(arg) = split(arg, ':') .|> string

# import the appropriate series for each element in the arglist
# blank it; do an exponential fitting
const testrange = 20:2:40
const offset = 0.0

# threshold for plotting purposes
const r2threshold = 0.975

const volloc = joinpath("/Volumes/ESSENTIALS", "CytCKinetics")
const dataloc = joinpath(volloc, "data")

function seriesfrominfo(seriesmeta, root, sub)
	ardata = open(joinpath(dataloc, root, sub, seriesmeta["filename"])) do io
		read(io, AbsorbanceRaw)
	end

	framestart, _ = seriesmeta["viewrange"]
	blanktime = seriesmeta["blankat"]
	return AbsorbanceSeries(ardata, framestart, 
							framestart+maximum(testrange), blanktime)
end

function r2exp(x, y, p)
	my = mean(y)
	sst = sum(y) do val
		abs2(my - val)
	end
	ssr = sum(zip(x, y)) do (v1, v2)
		py = p[1] * exp(-v1 * p[2])
		abs2(v2 - py)
	end

	return 1 - (ssr/sst)
end

function vecs2r2(times, values)
	ftmodel = CytCKinetics.expfit(times, values)
	p = ftmodel.param

	return p, r2exp(times, values, p)
end

# load information from metadata
series = map(arglist) do arg
	strain, set = split_arg(arg)
	
	metaall = open(joinpath(dataloc, "$strain.toml")) do io
		TOML.parse(io)
	end
	root = metaall["rootfolder"]
	setmeta = metaall[set]

	return seriesfrominfo(setmeta["Series"][argnum], root, setmeta["subfolder"])
end

@info "Proceeding to plotting routine..."

f = Figure(textsize=25,font="CMU Sans Serif",resolution=(1200,600))

# make grid layout
ga = f[1,1] = GridLayout()
gb = f[1,2] = GridLayout()

ax = Axis(gb[1,1],
	title = "Apparent R² against frame width (τ)",
	xlabel = "Frame width, τ (s)", ylabel = "Apparent R²"
)
ax2 = Axis(ga[1,1],
	title = "Exponential fitting, [Cyt c²⁺]₀ = 0.5μM",
	xlabel = "Time from start of reaction (s)", ylabel = "100× ΔAbs(t) (a.u.)"
)

# for each of the absorbance series blank it and get the appropriate range
for (i, ser) in enumerate(series)
	blanked = blankedframe(ser)
	zerotime = blanked.times .- minimum(blanked.times) .+ offset
	ΔAbs = blanked.values

	# get the testing range for the vectors
	testlen = map(testrange) do tval
		ceil(Int, tval / CytCKinetics.timestep(ser.ardata)) + 1
	end

	# do the fitting here and r2 calculation next to each other
	pairs = map(testlen) do tval
		ttime = zerotime[1:tval]; tvalues = ΔAbs[1:tval]
		return vecs2r2(ttime, tvalues)
	end

	r2vals = map(pairs) do (_, y)
		return y
	end
	kvals = map(pairs) do (x, _)
		return x[2]
	end
	params = map(pairs) do (x, _)
		return x
	end
	# add these to the axis
	lines!(ax, testrange, r2vals; label=printnames[i])
	scatter!(ax2, zerotime, ΔAbs; label=printnames[i], markersize=3)
	lines!(ax2, zerotime, params[end][1] .* exp.(-1 .* zerotime .* params[end][2]))

	# get statistics on the k values
	@info arglist[i] mean(kvals) std(kvals)
end

# draw the hline for the r2 threshold
hlines!(ax, r2threshold; color=:red, label="Threshold")

# draw the legend
leg = Legend(f[2,:], ax; 
			 framevisible=false, tellheight=true, orientation=:horizontal)

# letter labels
Label(ga[1,1,TopLeft()], "A", textsize=32)
Label(gb[1,1,TopLeft()], "B", textsize=32)

# save plot
const plotdump = joinpath(volloc, "plotdump/diagnostics/order/")
const sstring = joinpath(plotdump, "Set267.pdf")
save(sstring, f)
@info "Saved figure at $sstring"

