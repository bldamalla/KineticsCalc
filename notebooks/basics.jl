### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ d42324a2-279e-4f00-889f-7539efb6d41a
begin
	import Pkg;
	Pkg.activate("..");
	Pkg.instantiate();
end

# ╔═╡ 97d1ab2a-4a56-49b4-b67a-8633ab39c828
# we will need the following packages for reading
using CytCKinetics; import TOML

# ╔═╡ 6eeafb5b-f669-4f24-a6b3-be64ce692853
using CairoMakie

# ╔═╡ 3e28b968-22ef-11ed-2861-35bfb94e9fb2
md"# CytCKinetics.jl basics

This notebook presents the basic functions included in the package including:
+ Reading Hitachi spectrometer output files
+ Setting up `SeriesSet` and `SetGroup` objects
+ Reaction rate calculations and nonlinear regression with Michaelis-Menten model

Results will be plotted for demonstration purposes only.
"

# ╔═╡ 9abd23ff-957a-4749-8686-062702f0ac03
# define location from which to get metadata and raw absorbance input
const dataloc = joinpath(@__DIR__, "..", "data")

# ╔═╡ 6c9c3f82-7fd6-45d9-a707-718e919bf41e
md"
## Reading input and `SeriesSet` construction

Naming of separate kinetic assays (`SeriesSet`) has the format `strain:set` due
to historical reasons. Mostly because the package may be used if the author will
study effects of mutations on the oxidation activity. With this, I have chosen to
contain information from one strain into a single TOML file (except when it gets
_very_ large).

The mechanism for reading input will be explored in this section.
"

# ╔═╡ bfa22544-1ccd-481b-b561-b0fc05507dcb
"""
	argparse(str::AbstractString)

Split `strain:set` string into a tuple `(strain, set)` of strings.
"""
argparse(str) = split(str, ':') .|> string

# ╔═╡ d3cd92ba-f634-4e5c-861b-09c6ff6433cf
"""
	getsetmeta(strain, set)

Get metadata for the set of absorbance series data by cytochrome _c_ strain `strain`
and denoted by set label `set`.
"""
function getsetmeta(strain, set)
	# read strain metadata in data location
	strainmeta = open(joinpath(dataloc, "$strain.toml")) do io
		TOML.parse(io)
	end

	# obtain set information
	return strainmeta[set]
end

# ╔═╡ 4e929994-4670-4390-b64f-e3cbec934442
md"
The key function interfacing between this notebook and the package will be the
following function. Comments in the function definition will explain what the
lines do. In addition, paragraphs later will also discuss these more in detail.
"

# ╔═╡ cdaf0f8c-94a4-44ec-aba8-9fa5f755b2ab
"""
	meta2series(meta, rawsfolder)

From the metadata of the set, construct a `SeriesSet` containing `AbsorbanceSeries`
objects for downstream calculation of rate constants and kinetic parameters. The
argument `rawsfolder` is a path to the raw absorbance data.
"""
function meta2series(meta, rawsfolder; method=:linear)
	# series information is contained in "Series" key

	# get the initial substrate concentrations at which the measurements are done
	# units are (μM)
	concentrations = map(meta["Series"]) do ser
		ser["concentration"]
	end

	## get the raw absorbance and wrap with `AbsorbanceRaw`
	## other information for manipulation are in keys "viewrange" and "blankat"
	absraws = map(meta["Series"]) do ser
		fpath = joinpath(rawsfolder, ser["filename"])
		open(fpath) do io
			read(io, AbsorbanceRaw)
		end
	end
	manipdata = map(meta["Series"]) do ser
		starttime, stoptime = ser["viewrange"]
		blanktime = ser["blankat"]
		method === :linear && return (starttime, stoptime, blanktime)
		method === :exp && return (starttime, starttime+30, blanktime)
		throw(ArgumentError("Method $method not implemented here"))
	end

	# create the absorbance series objects
	absseries = map(zip(absraws, manipdata)) do (raws, manips)
		start, stop, blank = manips
		return AbsorbanceSeries(raws, start, stop, blank)
	end

	# create the series set
	return SeriesSet(absseries, concentrations, true, Dict{String,Any}())
end

# ╔═╡ 8bcc1551-fed3-4aa9-8b3f-d570f4a7f69a
md"
In this notebook, we will use `SeriesSets` obtained from `WT:Set1` and `WT:Set2` for
demonstration. We first define the list of assays we are working on, and perform a
series of calculations to get `AbsorbanceSeries` in `SeriesSet` objects.
"

# ╔═╡ d422d5d2-a177-4ff6-94f8-e114ad9008da
const arglist = ["WT:Set1", "WT:Set2"]

# ╔═╡ 3b6b19fd-32a8-4278-ba36-2c73b7af49be
set1series = let arg = arglist[1]
	strain, set = argparse(arg)

	# get metadata and folder for raw data
	meta = getsetmeta(strain, set)
	rawsfolder = joinpath(dataloc, strain, set)

	# construct the SeriesSet from metadata
	meta2series(meta, rawsfolder)
end

# ╔═╡ 0d515d07-d617-45af-918d-41ec5e82c1a9
md"
If the notebook is run, it can be seen that the type of the above created object is
`SeriesSet`. Its fields are:
+ `series`: a `Vector` of `AbsorbanceSeries` objects
+ `concentrations`: a `Vector` of `Float64`s
+ `menten`: indicating whether to perform Michaelis-Menten fitting for set points
+ `metadata`: a `Dict{String,Any}` containing other info (empty in this case)

Also `AbsorbanceSeries` objects have the following fields:
+ `ardata`: containing absorbance information (time of measurement and intensities)
+ `framestart` to `framestop`: time frame in which a rate constant is calculated
+ `blankat`: indicating time containing absorbance serving as blank

For this set, since `menten=true`, we can proceed to fitting.
"

# ╔═╡ 975a0fcf-ba2d-4f02-8237-c52ced1a67d5
md"
## Fitting and plotting

The data we have obtained previously can be fit and plotted as follows. The
packate provides a `fit` function that accepts a target results type, mostly
for clarity upon reading by the user. See the documentation for the keyword
arguments accepted.
"

# ╔═╡ 8fbe21d0-bc2c-4058-b939-0585ef165cee
set1fitresults = CytCKinetics.fit(SeriesSetResults, set1series; method=:exp)

# ╔═╡ bcc900d8-920d-47aa-ba4c-4ca6acdcf177
md"
As can be seen for the above fitting, the function complains that the ``R^2`` value
for the calculation of a rate constant is less than the threshold (that I set) of
``0.96``. A value slightly less than this threshold usually (by experience) that the
data is _relatively noisy_ (low signal-to-noise). In this case, this is given by
measurements at ``0.3\text{μM}``.

Cases where ``R^2 < 0.9`` usually occur when the manipulation parameters in a
`AbsorbanceSeries` object is not properly chosen. It has been shown by previous
reports that cytochrome _c_ oxidation by its oxidase is first-order over a wide
range of concentrations in _in vitro_ assays, so it is unlikely that the fitting
model is at fault.
"

# ╔═╡ 9e219bdd-7401-45f5-b072-2cec1ae95dce
md"
The same procedures can be done for the other series set.
"

# ╔═╡ 4d544fc3-14d7-458c-850b-ba2d812ea842
set2series = let arg = arglist[2]
	strain, set = argparse(arg)

	meta = getsetmeta(strain, set)
	rawsfolder = joinpath(dataloc, strain, set)

	meta2series(meta, rawsfolder)
end

# ╔═╡ 417e6086-095b-43b6-9ec8-364afe39cffb
set2fitresults = CytCKinetics.fit(SeriesSetResults, set2series; method=:exp)

# ╔═╡ a3b59dc9-e84a-4eb7-bc38-e7269f56a7e0
md"
Finally we can plot the fitted results!
"

# ╔═╡ cce345e7-cc71-4232-af9c-aec31b8c9fd8
plottheme = Theme(
	resolution = (900, 600),
	textsize = 35,

	Axis = (
		title = "Reaction rate and fitted curves", 
		xlabel = "[Cyt c²⁺]₀ (μM)", ylabel = "Initial rate (nM/s)",
		xticks = 0.0:1.0:8.0, yticks = 0:15:75,
		xminorticks = IntervalsBetween(2), xminorticksvisible = true,
		yminorticks = IntervalsBetween(3), yminorticksvisible = true
	),

	Scatter = (
		markersize = 6, colormap = :seaborn_colorblind
	),

	Lines = (
		strokewidth = 2, colormap = :seaborn_colorblind
	)
)

# ╔═╡ 2c039d5a-ecaf-46d2-b4c1-6bf2b429addb
with_theme(plottheme) do
	(; concentrations, initrates, fitparams) = set1fitresults
	fig1, _, _ = scatter(concentrations, initrates; label="Set1")
	lines!(0.0:0.01:7.5, x->CytCKinetics.menten(x, fitparams); label="Set1")

	(; concentrations, initrates, fitparams) = set2fitresults
	scatter!(concentrations, initrates; label="Set2")
	lines!(0.0:0.01:7.5, x->CytCKinetics.menten(x, fitparams); label="Set2")

	axislegend("Sets"; merge=true, position=:rc, orientation=:horizontal)
	return fig1
end

# ╔═╡ Cell order:
# ╟─d42324a2-279e-4f00-889f-7539efb6d41a
# ╟─3e28b968-22ef-11ed-2861-35bfb94e9fb2
# ╟─9abd23ff-957a-4749-8686-062702f0ac03
# ╟─6c9c3f82-7fd6-45d9-a707-718e919bf41e
# ╠═97d1ab2a-4a56-49b4-b67a-8633ab39c828
# ╟─bfa22544-1ccd-481b-b561-b0fc05507dcb
# ╟─d3cd92ba-f634-4e5c-861b-09c6ff6433cf
# ╟─4e929994-4670-4390-b64f-e3cbec934442
# ╟─cdaf0f8c-94a4-44ec-aba8-9fa5f755b2ab
# ╟─8bcc1551-fed3-4aa9-8b3f-d570f4a7f69a
# ╠═d422d5d2-a177-4ff6-94f8-e114ad9008da
# ╠═3b6b19fd-32a8-4278-ba36-2c73b7af49be
# ╟─0d515d07-d617-45af-918d-41ec5e82c1a9
# ╟─975a0fcf-ba2d-4f02-8237-c52ced1a67d5
# ╠═8fbe21d0-bc2c-4058-b939-0585ef165cee
# ╟─bcc900d8-920d-47aa-ba4c-4ca6acdcf177
# ╟─9e219bdd-7401-45f5-b072-2cec1ae95dce
# ╠═4d544fc3-14d7-458c-850b-ba2d812ea842
# ╠═417e6086-095b-43b6-9ec8-364afe39cffb
# ╟─a3b59dc9-e84a-4eb7-bc38-e7269f56a7e0
# ╠═6eeafb5b-f669-4f24-a6b3-be64ce692853
# ╟─cce345e7-cc71-4232-af9c-aec31b8c9fd8
# ╠═2c039d5a-ecaf-46d2-b4c1-6bf2b429addb
