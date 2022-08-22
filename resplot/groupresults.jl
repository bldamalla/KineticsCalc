# scripts/groupresults.jl --- for obtaining results for set groups

## actually there are two ways to do this:
## 1. a short way using set results
## 2. from scratch from creating the series sets

import TOML
using CytCKinetics

# const dataloc = joinpath(@__DIR__, "..", "data")
const dataloc = joinpath("/Volumes/ESSENTIALS", "CytCKinetics", "data")
arglist = ["WTC:GroupC", "WTC:GroupD", "WTC:GroupCD"]
split_arg(arg) = split(arg, ':') .|> string
function getsetmeta(strain, set)
    metaall = open(joinpath(dataloc, "$strain.toml")) do io
        TOML.parse(io)
    end
    root = metaall["rootfolder"]

    return metaall[set], root
end
function seriesfrominfo(seriesmeta, root, sub; method=:exp)
    ardata = open(joinpath(dataloc, root, sub, seriesmeta["filename"])) do io
        read(io, AbsorbanceRaw)
    end

    framestart, framestop = seriesmeta["viewrange"]
    blanktime = seriesmeta["blankat"]
    method === :line && return AbsorbanceSeries(ardata, framestart, framestop, blanktime)
    method === :exp && return AbsorbanceSeries(ardata, framestart, framestart+30, blanktime)
    throw(ArgumentError("`method`=$(method) not supported."))
end

# open the group meta and look for the groups you are interested in
groupmeta = open(joinpath(dataloc, "Groups.toml")) do io
    TOML.parse(io)
end
groupsets = map(arglist) do arg
    strain, group = split_arg(arg)

    return strain, groupmeta[strain][group]["members"]
end
@info "Retrieved group information"

setgroups = map(groupsets) do (strain, members)
    # get a vector of the setmetadata and the rootfolder of the strain
    request = map(members) do member
        getsetmeta(strain, member)
    end

    # get series information and construct the seriessets
    request2 = map(request) do setmeta
        meta, root = setmeta
        series_ = map(meta["Series"]) do object
            seriesfrominfo(object, root, meta["subfolder"]; method=:exp)
        end
        concentrations_ = map(meta["Series"]) do object
            object["concentration"]
        end
        menten = false  # actually this is just temporary anyway

        return SeriesSet(series_, concentrations_, menten, Dict{String,Any}())
    end

    # now merge the series sets for that group
    return CytCKinetics.setmerge(SetGroup, request2...)
end

# now fit the merged sets
fitresults = fit.(SetGroupResults, setgroups; method=:exp, offset=0)
@info "Fitted calculations..."

# serialize into the results toml

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
    strain, group = split_arg(arg)
    if !haskey(resultsmeta, strain)
        push!(resultsmeta, strain=>Dict{String,Any}())
    end
    strainmeta = resultsmeta[strain]
    strainmeta[group] = CytCKinetics.serialize(object)
end

open(resultsmetaloc, "w") do io
    TOML.print(io, resultsmeta, sorted=true, by=cmp_)
end
@info "Written results in $resultsmetaloc"
