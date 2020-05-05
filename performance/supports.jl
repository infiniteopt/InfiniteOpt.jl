# include("C:/Users/pulsipher/.julia/dev/InfiniteOpt/performance/supports.jl")
using DataStructures, BenchmarkTools

################################################################################
## Test multi-dimensional support performance with various types
################################################################################
## Make the datatypes
# Matrix
struct MatSupps
    supps::Matrix{Float64} # columns are supports
    labels::Vector{Set{Symbol}}
end

# Vector of Vectors
struct VVSupps
    supps::Vector{Vector{Float64}}
    labels::Vector{Set{Symbol}}
end

# Vector of Tuples
struct VTSupps{N}
    supps::Vector{NTuple{N, Float64}}
    labels::Vector{Set{Symbol}}
end

# Dictionary of Vectors
struct DVSupps
    supps::Dict{Vector{Float64}, Set{Symbol}}
end

# Dictionary of Tuples
struct DTSupps{N}
    supps::Dict{NTuple{N, Float64}, Set{Symbol}}
end

# Sorted Dictionary
struct SDSupps
    supps::DataStructures.SortedDict{Vector{Float64}, Set{Symbol}}
end

## Define methods to test unsorted and ununique addition with matrix input
# MatSupps
function add1(d::MatSupps, supps::Matrix, label::Symbol)::MatSupps
    new_supps = hcat(d.supps, supps)
    new_labels = append!(d.labels, [Set([label]) for i in 1:size(supps, 2)])
    return MatSupps(new_supps, new_labels)
end

# VVSupps
function add1(d::VVSupps, supps::Matrix, label::Symbol)::VVSupps
    for i in 1:size(supps, 2)
        push!(d.supps, supps[:, i])
        push!(d.labels, Set([label]))
    end
    return d
end

# VTSupps
function add1(d::VTSupps{N}, supps::Matrix, label::Symbol)::VTSupps{N} where {N}
    for i in 1:size(supps, 2)
        push!(d.supps, Tuple(supps[j, i] for j in 1:size(supps, 1)))
        push!(d.labels, Set([label]))
    end
    return d
end

# DVSupps
function add1(d::DVSupps, supps::Matrix, label::Symbol)::DVSupps
    for i in 1:size(supps, 2)
        s = @view(supps[:, i])
        if haskey(d.supps, s)
            push!(d.supps[s], label)
        else
            d.supps[s] = Set([label])
        end
    end
    return d
end

function add1_2(d::DVSupps, supps::Matrix, label::Symbol)::DVSupps
    d2 = Dict(@views supps[:, i] => Set([label]) for i in 1:size(supps, 2))
    merge!(union!, d.supps, d2)
    return d
end

# DTSupps
function add1(d::DTSupps{N}, supps::Matrix, label::Symbol)::DTSupps where {N}
    for i in 1:size(supps, 2)
        s = Tuple(supps[j, i] for j in 1:size(supps, 1))
        if haskey(d.supps, s)
            push!(d.supps[s], label)
        else
            d.supps[s] = Set([label])
        end
    end
    return d
end

# SDSupps
function add1(d::SDSupps, supps::Matrix, label::Symbol)::SDSupps
    for i in 1:size(supps, 2)
        s = supps[:, i]
        if haskey(d.supps, s)
            push!(d.supps[s], label)
        else
            d.supps[s] = Set([label])
        end
    end
    return d
end

## Define methods to test unsorted and unique addition with matrix input
# unique addition function
function _unique_cols(a::Matrix)
    if size(a, 2) < 2
        return a
    else
        inds = [1]
        prev = a[:, 1]
        for i in 2:size(a, 2)
            addit = true
            curr = a[:, i]
            if curr == prev
                continue
            end
            for j in inds
                if a[:, i] == a[:, j]
                    addit = false
                    break
                end
            end
            if addit
                push!(inds, i)
            end
            prev = curr
        end
        return a[:, inds]
    end
end

# MatSupps
function add2(d::MatSupps, supps::Matrix, label::Symbol)::MatSupps
    new_supps = _unique_cols(hcat(d.supps, supps))
    new_labels = append!(d.labels, [Set([label]) for i in 1:size(new_supps, 2)])
    return MatSupps(new_supps, new_labels)
end

# VVSupps
function add2(d::VVSupps, supps::Matrix, label::Symbol)::VVSupps
    for col in eachcol(supps)
        idx = findfirst(isequal(col), d.supps)
        if idx === nothing
            push!(d.supps, col)
            push!(d.labels, Set([label]))
        else
            push!(d.labels[idx], label)
        end
    end

    return d
end

# VTSupps
function add2(d::VTSupps{N}, supps::Matrix, label::Symbol)::VTSupps{N} where {N}
    for i in 1:size(supps, 2)
        s = Tuple(supps[j, i] for j in 1:size(supps, 1))
        if !(s in d.supps)
            push!(d.supps, s)
            push!(d.labels, Set([label]))
        end
    end
    return d
end

# DVSupps
function add2(d::DVSupps, supps::Matrix, label::Symbol)::DVSupps
    return add1(d, supps, label)
end

# DTSupps
function add2(d::DTSupps{N}, supps::Matrix, label::Symbol)::DTSupps where {N}
    return add1(d, supps, label)
end

# SDSupps
function add2(d::SDSupps, supps::Matrix, label::Symbol)::SDSupps
    return add1(d, supps, label)
end

## Make methods to test iteration speed over entire set
# MatSupps
function loopall(d::MatSupps)::Nothing
    i = 1
    for col in eachcol(d.supps)
        a = (col, d.labels[i])
        i += 1
    end
    return
end

# VVSupps
function loopall(d::VVSupps)::Nothing
    for i in eachindex(d.supps)
        a = (d.supps[i], d.labels[i])
    end
    return
end

# VTSupps
function loopall(d::VTSupps)::Nothing
    for i in eachindex(d.supps)
        a = (d.supps[i], d.labels[i])
    end
    return
end

# DVSupps
function loopall(d::DVSupps)::Nothing
    for (s, l) in d.supps
        a = (s, l)
    end
    return
end

# DTSupps
function loopall(d::DTSupps)::Nothing
    for (s, l) in d.supps
        a = (s, l)
    end
    return
end

# SDSupps
function loopall(d::SDSupps)::Nothing
    for (s, l) in d.supps
        a = (s, l)
    end
    return
end

## Test filtered iterate
# MatSupps
function loopsome(d::MatSupps, label::Symbol)::Nothing
    for i in findall(e -> label in e, d.labels)
        a = (d.supps[:, i], d.labels[i])
    end
    return
end

# VVSupps
function loopsome(d::VVSupps, label::Symbol)::Nothing
    # for i in findall(e -> label in e, d.labels)
    #     a = (d.supps[i], d.labels[i])
    # end
    for i in eachindex(d.supps)
        if label in d.labels[i]
            a = (d.supps[i], d.labels[i])
        end
    end
    return
end

# VTSupps
function loopsome(d::VTSupps, label::Symbol)::Nothing
    for i in findall(e -> label in e, d.labels)
        a = (d.supps[i], d.labels[i])
    end
    return
end

# DVSupps
function loopsome(d::DVSupps, label::Symbol)::Nothing
    # for (s, l) in filter(p -> label in p[2], d.supps)
    #     a = (s, l)
    # end
    for (s, l) in d.supps
        if label in l
            a = (s, l)
        end
    end
    return
end

# DTSupps
function loopsome(d::DTSupps, label::Symbol)::Nothing
    for (s, l) in filter(p -> label in p[2], d.supps)
        a = (s, l)
    end
    return
end

# SDSupps
function loopsome(d::SDSupps, label::Symbol)::Nothing
    for (s, l) in filter(p -> label in p[2], d.supps)
        a = (s, l)
    end
    return
end

## Test extracting all as matrix
# MatSupps
function suppmat(d::MatSupps)::Matrix
    return d.supps
end

# VVSupps
function suppmat(d::VVSupps)::Matrix
    supps = Matrix{Float64}(undef, length(first(d.supps)), length(d.supps))
    i = 1
    for s in d.supps
        supps[:, i] = s
        i += 1
    end
    return supps
end

# VTSupps
function suppmat(d::VTSupps)::Matrix
    supps = Matrix{Float64}(undef, length(first(d.supps)), length(d.supps))
    i = 1
    for s in d.supps
        supps[:, i] = [s...]
        i += 1
    end
    return supps
end

# DVSupps
function suppmat(d::DVSupps)::Matrix
    supps = Matrix{Float64}(undef, length(first(keys(d.supps))), length(d.supps))
    i = 1
    for s in keys(d.supps)
        supps[:, i] = s
        i += 1
    end
    return supps
end

# DTSupps
function suppmat(d::DTSupps)::Matrix
    supps = Matrix{Float64}(undef, length(first(keys(d.supps))), length(d.supps))
    i = 1
    for s in keys(d.supps)
        supps[:, i] = [s...]
        i += 1
    end
    return supps
end

# SDSupps
function suppmat(d::SDSupps)::Matrix
    supps = Matrix{Float64}(undef, length(first(keys(d.supps))), length(d.supps))
    i = 1
    for s in keys(d.supps)
        supps[:, i] = s
        i += 1
    end
    return supps
end

## Test methods for getting just some of the supports
# VVSupps
function somesupps(d::VVSupps, label::Symbol)::Matrix
    inds = findall(e -> label in e, d.labels)
    supps = Matrix{Float64}(undef, length(first(d.supps)), length(inds))
    for i in eachindex(inds)
        supps[:, i] = d.supps[inds[i]]
    end
    return supps
end

# DVSupps
function somesupps(d::DVSupps, label::Symbol)::Matrix
    new_supps = filter(p -> label in p[2], d.supps)
    supps = Matrix{Float64}(undef, length(first(keys(new_supps))), length(new_supps))
    i = 1
    for s in keys(new_supps)
        supps[:, i] = s
        i += 1
    end
    return supps
end

function somesupps2(d::DVSupps, label::Symbol)::Matrix
    new_supps = findall(e -> label in e, d.supps)
    supps = Matrix{Float64}(undef, length(first(new_supps)), length(new_supps))
    i = 1
    for s in new_supps
        supps[:, i] = s
        i += 1
    end
    return supps
end

function somesupps3(d::DVSupps, label::Symbol)::Matrix
    # new_supps = findall(e -> label in e, d.supps)
    return reduce(hcat, findall(e -> label in e, d.supps))
end

function somesupps4(d::DVSupps, label::Symbol)::Matrix
    supp_vect = Float64[]
    counter = 0
    for (supp, labels) in d.supps
        if label in labels
            append!(supp_vect, supp)
            counter += 1
        end
    end
    return reshape(supp_vect, length(first(d.supps)[1]), counter)
end

## Test the performance of each
# Prepare the additions
a1 = round.(rand(10, 10000), sigdigits = 8)
a2 = round.(hcat(a1[:, rand(1:10000, 5000)], rand(10, 5000)), sigdigits = 8)
num_samples = 20

# Test simple addition
add1_table = Matrix{Any}(undef, 7, 4)
add1_table[1, :] = ["Type", "Time", "Memory", "Allocs"]
add1_table[2:end, 1] = ["MatSupps", "VVSupps", "VTSupps", "DVSupps", "DTSupps", "SDSupps"]
t = @benchmark(add1(T, $a1, :a), evals = 1, samples = num_samples, seconds = 1, setup = (T = MatSupps(zeros(10, 1), [Set([:a])])))
med = BenchmarkTools.median(t)
add1_table[2, 2] = med.time
add1_table[2, 3] = med.memory
add1_table[2, 4] = med.allocs
t = @benchmark(add1(T, $a1, :a), evals = 1, samples = num_samples, seconds = 1, setup = (T = VVSupps([zeros(10)], [Set([:a])])))
med = BenchmarkTools.median(t)
add1_table[3, 2] = med.time
add1_table[3, 3] = med.memory
add1_table[3, 4] = med.allocs
t = @benchmark(add1(T, $a1, :a), evals = 1, samples = num_samples, seconds = 1, setup = (T = VTSupps([(zeros(10)...,)], [Set([:a])])))
med = BenchmarkTools.median(t)
add1_table[4, 2] = med.time
add1_table[4, 3] = med.memory
add1_table[4, 4] = med.allocs
t = @benchmark(add1(T, $a1, :a), evals = 1, samples = num_samples, seconds = 1, setup = (T = DVSupps(Dict(zeros(10) => Set([:a])))))
med = BenchmarkTools.median(t)
add1_table[5, 2] = med.time
add1_table[5, 3] = med.memory
add1_table[5, 4] = med.allocs
t = @benchmark(add1(T, $a1, :a), evals = 1, samples = num_samples, seconds = 1, setup = (T = DTSupps(Dict((zeros(10)...,) => Set([:a])))))
med = BenchmarkTools.median(t)
add1_table[6, 2] = med.time
add1_table[6, 3] = med.memory
add1_table[6, 4] = med.allocs
t = @benchmark(add1(T, $a1, :a), evals = 1, samples = num_samples, seconds = 1, setup = (T = SDSupps(SortedDict(zeros(10) => Set([:a])))))
med = BenchmarkTools.median(t)
add1_table[7, 2] = med.time
add1_table[7, 3] = med.memory
add1_table[7, 4] = med.allocs

# Test unique addition
add2_table = Matrix{Any}(undef, 7, 4)
add2_table[1, :] = ["Type", "Time", "Memory", "Allocs"]
add2_table[2:end, 1] = ["MatSupps", "VVSupps", "VTSupps", "DVSupps", "DTSupps", "SDSupps"]
t = @benchmark(add2(T, $a2, :b), evals = 1, samples = num_samples, seconds = 1, setup = (T = add1(MatSupps(zeros(10, 1), [Set([:a])]), $a1, :a)))
med = BenchmarkTools.median(t)
add2_table[2, 2] = med.time
add2_table[2, 3] = med.memory
add2_table[2, 4] = med.allocs
t = @benchmark(add2(T, $a2, :b), evals = 1, samples = num_samples, seconds = 1, setup = (T = add1(VVSupps([zeros(10)], [Set([:a])]), $a1, :a)))
med = BenchmarkTools.median(t)
add2_table[3, 2] = med.time
add2_table[3, 3] = med.memory
add2_table[3, 4] = med.allocs
t = @benchmark(add2(T, $a2, :b), evals = 1, samples = num_samples, seconds = 1, setup = (T = add1(VTSupps([(zeros(10)...,)], [Set([:a])]), $a1, :a)))
med = BenchmarkTools.median(t)
add2_table[4, 2] = med.time
add2_table[4, 3] = med.memory
add2_table[4, 4] = med.allocs
t = @benchmark(add2(T, $a2, :b), evals = 1, samples = num_samples, seconds = 1, setup = (T = add1(DVSupps(Dict(zeros(10) => Set([:a]))), $a1, :a)))
med = BenchmarkTools.median(t)
add2_table[5, 2] = med.time
add2_table[5, 3] = med.memory
add2_table[5, 4] = med.allocs
t = @benchmark(add2(T, $a2, :b), evals = 1, samples = num_samples, seconds = 1, setup = (T = add1(DTSupps(Dict((zeros(10)...,) => Set([:a]))), $a1, :a)))
med = BenchmarkTools.median(t)
add2_table[6, 2] = med.time
add2_table[6, 3] = med.memory
add2_table[6, 4] = med.allocs
t = @benchmark(add2(T, $a2, :b), evals = 1, samples = num_samples, seconds = 1, setup = (T = add1(SDSupps(SortedDict(zeros(10) => Set([:a]))), $a1, :a)))
med = BenchmarkTools.median(t)
add2_table[7, 2] = med.time
add2_table[7, 3] = med.memory
add2_table[7, 4] = med.allocs

# Test looping over all
T1 = add2(add1(MatSupps(zeros(10, 1), [Set([:a])]), a1, :a), a2, :b)
T2 = add2(add1(VVSupps([zeros(10)], [Set([:a])]), a1, :a), a2, :b)
T3 = add2(add1(VTSupps([(zeros(10)...,)], [Set([:a])]), a1, :a), a2, :b)
T4 = add2(add1(DVSupps(Dict(zeros(10) => Set([:a]))), a1, :a), a2, :b)
T5 = add2(add1(DTSupps(Dict((zeros(10)...,) => Set([:a]))), a1, :a), a2, :b)
T6 = add2(add1(SDSupps(SortedDict(zeros(10) => Set([:a]))), a1, :a), a2, :b)
Ts = [T1, T2, T3, T4, T5, T6]

loopall_table = Matrix{Any}(undef, 7, 4)
loopall_table[1, :] = ["Type", "Time", "Memory", "Allocs"]
loopall_table[2:end, 1] = ["MatSupps", "VVSupps", "VTSupps", "DVSupps", "DTSupps", "SDSupps"]
for i in eachindex(Ts)
    T = Ts[i]
    t = @benchmark(loopall($T), evals = 1, samples = num_samples, seconds = 1)
    med = BenchmarkTools.median(t)
    loopall_table[i+1, 2] = med.time
    loopall_table[i+1, 3] = med.memory
    loopall_table[i+1, 4] = med.allocs
end

# Test looping over some
loopsome_table = Matrix{Any}(undef, 7, 4)
loopsome_table[1, :] = ["Type", "Time", "Memory", "Allocs"]
loopsome_table[2:end, 1] = ["MatSupps", "VVSupps", "VTSupps", "DVSupps", "DTSupps", "SDSupps"]
for i in eachindex(Ts)
    T = Ts[i]
    t = @benchmark(loopsome($T, :a), evals = 1, samples = num_samples, seconds = 1)
    med = BenchmarkTools.median(t)
    loopsome_table[i+1, 2] = med.time
    loopsome_table[i+1, 3] = med.memory
    loopsome_table[i+1, 4] = med.allocs
end

# Test extraction
extract_table = Matrix{Any}(undef, 7, 4)
extract_table[1, :] = ["Type", "Time", "Memory", "Allocs"]
extract_table[2:end, 1] = ["MatSupps", "VVSupps", "VTSupps", "DVSupps", "DTSupps", "SDSupps"]
for i in eachindex(Ts)
    T = Ts[i]
    t = @benchmark(suppmat($T), evals = 1, samples = num_samples, seconds = 1)
    med = BenchmarkTools.median(t)
    extract_table[i+1, 2] = med.time
    extract_table[i+1, 3] = med.memory
    extract_table[i+1, 4] = med.allocs
end
