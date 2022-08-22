# scripts/setresults.jl -- calculation of set results from raw data

# made so you don't have to work in the repl anymore

using CytCKinetics
import TOML

# const dataloc = joinpath(@__DIR__, "..", "data")
const dataloc = joinpath("/Volumes/ESSENTIALS", "CytCKinetics", "data")

arglist = ["WTC:Set5", "WTC:Set6", "WTC:Set7", "WTC:Set8"]
split_arg(arg) = split(arg, ':') .|> string
function getsetmeta(strain, set)
    metaall = open(joinpath(dataloc, "$strain.toml")) do io
        TOML.parse(io)
    end
    root = metaall["rootfolder"]

    # get the rootfolder of the metadata for the set
    # apparently, I did not include it in the actual set meta for
    # brevity and reusability
    return metaall[set], root
end
function seriesfrominfo(seriesmeta, root, sub)
    ardata = open(joinpath(dataloc, root, sub, seriesmeta["filename"]))do io
        read(io, AbsorbanceRaw)
    end

    framestart, framestop = seriesmeta["viewrange"]
    blanktime = seriesmeta["blankat"]
    return AbsorbanceSeries(ardata, framestart, framestop, blanktime)
end

# load series set information from metadata
metadata = map(arglist) do arg
    strain, set = split_arg(arg)
    return getsetmeta(strain, set)
end
@info "Loaded metadata for each set..."

# abstract series information and create series objects
seriessets = map(metadata) do setmeta
    meta, root = setmeta
    series_ = map(meta["Series"]) do object
        seriesfrominfo(object, root, meta["subfolder"])
    end
    concentrations_ = map(meta["Series"]) do object
        object["concentration"]
    end
    menten_ = get(meta, "menten", false)

    ## for now you probably don't need the other metadata in the
    ## series set definition
    return SeriesSet(series_, concentrations_, menten_, Dict{String,Any}())
end

# fit results using the internal fit in CytCKinetics
fitresults = fit.(SeriesSetResults, seriessets; method=:exp, offset=0)
@info "Finished calculations..."

# serialize into the results toml

## first is you read the results toml and then deserialize all to 
## seriessetresults (probably make handling for setgroups different)
const resultsmetaloc = joinpath(dataloc, "results.toml")

resultsmeta = open(resultsmetaloc) do io
    TOML.parse(io)
end

# add the serialized fitting results to the results meta file
# print in alphabetical order?
cmp_(str1::String, str2::String) = cmp(str1, str2)
cmp_(str::String) = identity(str)

paired = zip(arglist, fitresults)

for (arg, object) in paired
    strain, set = split_arg(arg)
    if !haskey(resultsmeta, strain)
        push!(resultsmeta, strain=>Dict{String,Any}())
    end
    strainmeta = resultsmeta[strain]
    strainmeta[set] = CytCKinetics.serialize(object)
end

open(resultsmetaloc, "w") do io
    TOML.print(io, resultsmeta, sorted=true, by=cmp_)
end
@info "Written results in $resultsmetaloc"
