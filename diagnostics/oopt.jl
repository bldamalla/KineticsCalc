# diagnostics/oopt.jl --- optimize offsetting for each absorbance series

# working input will be the set names

import TOML
using CytCKinetics, Statistics

const volloc = joinpath("/Volumes/ESSENTIALS/CytCKinetics")
const datapath = joinpath(volloc, "data")

# test range for possible values of the offset
offsettrange = 0.0:0.1:5.0
function modelr2(f, x, y)
    mval = mean(y)

    sst = sum(y) do val
		abs2(val - mval)
    end
    ssr = sum(zip(x, y)) do (vx, vy)
		abs2(vy - f(vx))
    end

    return 1 - (ssr/sst)
end
function blankedtor2(ser::CytCKinetics.ARAny, offset)
    zeroed_times = ser.times .- minimum(ser.times) .+ offset

    ftmodel = CytCKinetics.expfit(zeroed_times, ser.values)
    p = ftmodel.param

    return modelr2(zeroed_times, ser.values) do x
		p[1] * exp(-x * p[2])
    end
end

arglist = [ "WTC:Set5", "WTC:Set6", "WTC:Set7", "WTC:Set8" ]
parsearg(arg) = split(arg, ':') .|> string

# obtain the series sets from reading
seriessets = map(arglist) do arg
    strain, set = parsearg(arg)

    metaloc = joinpath(datapath, string(strain, ".toml"))
    setmeta = open(metaloc) do io   # basically read the metadata for the strain
		TOML.parse(io)		    # and store as a dictionary
    end |> q->q[set]		    # obtain the one specifically for the sets

    # create the absorbance series objects
    absseries = map(setmeta["Series"]) do series
	start, stop = series["viewrange"]
	blank = series["blankat"]
	
	# this is a bit cheating because I know the format of the folders
	absdatapath = joinpath(datapath, strain, set, series["filename"])
	ardata = open(absdatapath) do io
	    read(io, AbsorbanceRaw)
	end

	# return AbsorbanceSeries(ardata, start, stop, blank)	    # use the ones set for threshold finding
	return AbsorbanceSeries(ardata, start, start+30, blank)	    # use the initial 150 points for testing
    end

    concentrations = map(setmeta["Series"]) do series
		series["concentration"] |> float
    end

    menten = setmeta["menten"]

    # just an empty string dictionary to satisfy the constructor
    edict = Dict{String,Any}()

    return SeriesSet(absseries, concentrations, menten, edict)
end

# loop through the seriessets and work on get
# the optimal value of offset for each set
for (arg, seriesset) in zip(arglist, seriessets)
    @info arg

    # blank the seriesset
    (; series, concentrations) = seriesset
    for (series_, conc) in zip(series, concentrations)
		ser = series_ |> blankedframe
		opt_offset = argmax(offsettrange) do test
			blankedtor2(ser, test)
		end

		@show conc
		@show opt_offset
		@show blankedtor2(ser, opt_offset)
    end
end

