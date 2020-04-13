"""
    DualDict{K1, K2, V}

A convenient dictionary type that efficiently implements a dictionary with
two different key types `K1` and `K2` for one common value type `V`. General
usage is experimental and this is primarily intended to be used to store
`ParmaterBounds`.

**Example**
```julia-repl
julia> DualDict{String, Int, Float64}() # Initialize empty DualDict
DualDict{String,Int64,Float64} with 0 entries:

julia> dd = DualDict{String, Int, Float64}("a" => 2.1, 3 => 5.7, "v" => 3.1)
DualDict{String,Int64,Float64} with 3 entries:
  "v" => 3.1
  "a" => 2.1
  3 => 5.7

julia> dd[3] # Index normally
5.7

julia> dd["v"] = 42 # Set values normally
42

julia> dd[2] = 8.3 # Add new pairs
8.3

julia> dd
DualDict{String,Int64,Float64} with 4 entries:
  "v" => 42.0
  "a" => 2.1
  2 => 8.3
  3 => 5.7

julia> delete!(dd, "a") # Delete pairs as usual
DualDict{String,Int64,Float64} with 3 entries:
  "v" => 42.0
  2 => 8.3
  3 => 5.7

julia> haskey(dd, "a") # Use haskey
false

julia> [k => v for (k, v) in dd] # Iterate as normal Dict
3-element Array{Pair{A,Float64} where A,1}:
 "v" => 42.0
   2 => 8.3
   3 => 5.7

julia> empty!(dd) # Empty it out
DualDict{String,Int64,Float64} with 0 entries:
```
"""
struct DualDict{K1, K2, V}
    dict1::Dict{K1, V}
    dict2::Dict{K2, V}
end

# Default constructor
function DualDict{K1, K2, V}() where {K1, K2, V}
    return DualDict{K1, K2, V}(Dict{K1, V}(), Dict{K2, V}())
end

# Tuple constructor
function DualDict{K1, K2, V}(tup::NTuple{N, Pair{<:Union{K1, K2}, <:V}}
                             ) where {N, K1, K2, V}
    pairs1 = Iterators.filter(x -> x[1] isa K1, tup)
    pairs2 = Iterators.filter(x -> x[1] isa K2, tup)
    return DualDict{K1, K2, V}(Dict{K1, V}(pairs1), Dict{K2, V}(pairs2))
end

# Pair constructor
function DualDict{K1, K2, V}(pairs::Pair{<:Union{K1, K2}, <:V}...
                             ) where {K1, K2, V}
    return DualDict{K1, K2, V}(pairs)
end

# Extend Base.:(==)
function Base.:(==)(dd1::DualDict, dd2::DualDict)::Bool
    return dd1.dict1 == dd2.dict1 && dd1.dict2 == dd2.dict2
end

# Extend Base.copy
function Base.copy(dd::DualDict{K1, K2, V})::DualDict{K1, K2, V} where {K1, K2, V}
    return DualDict(copy(dd.dict1), copy(dd.dict2))
end

# Extend Base.in when checking the first dict
function Base.in(pair::Pair{K1, V}, dd::DualDict{K1, K2, V})::Bool where {K1, K2, V}
    return pair in dd.dict1
end

# Extend Base.in when checking the second dict
function Base.in(pair::Pair{K2, V}, dd::DualDict{K1, K2, V})::Bool where {K1, K2, V}
    return pair in dd.dict2
end

# Extend Base.haskey when the index is that of dict1
function Base.haskey(dd::DualDict{K1, K2, V}, index::K1)::Bool where {K1, K2, V}
    return haskey(dd.dict1, index)
end

# Extend Base.haskey when the index is that of dict1
function Base.haskey(dd::DualDict{K1, K2, V}, index::K2)::Bool where {K1, K2, V}
    return haskey(dd.dict2, index)
end

# Extend Base.length
Base.length(dd::DualDict)::Int = length(dd.dict1) + length(dd.dict2)

# Extend Base.isempty
Base.isempty(dd::DualDict)::Bool = isempty(dd.dict1) && isempty(dd.dict2)

# Extend Base.getindex when the index is that of dict1
function Base.getindex(dd::DualDict{K1, K2, V}, index::K1)::V where {K1, K2, V}
    haskey(dd, index) || throw(BoundsError("Index $index not found in $dd."))
    return dd.dict1[index]
end

# Extend Base.getindex when the index is that of dict2
function Base.getindex(dd::DualDict{K1, K2, V}, index::K2)::V where {K1, K2, V}
    haskey(dd, index) || throw(BoundsError("Index $index not found in $dd."))
    return dd.dict2[index]
end

# Extend Base.setindex! when the index is that of dict1
function Base.setindex!(dd::DualDict{K1, K2, V}, value,
                        index::K1)::V where {K1, K2, V}
    return dd.dict1[index] = value
end

# Extend Base.setindex! when the index is that of dict2
function Base.setindex!(dd::DualDict{K1, K2, V}, value,
                        index::K2)::V where {K1, K2, V}
    return dd.dict2[index] = value
end

# Extend Base.empty!
function Base.empty!(dd::DualDict{K1, K2, V})::DualDict{K1, K2, V} where {K1, K2, V}
    empty!(dd.dict1)
    empty!(dd.dict2)
    return dd
end

# Extend Base.delete! when the index is that of dict1
function Base.delete!(dd::DualDict{K1, K2, V},
                      index::K1)::DualDict{K1, K2, V} where {K1, K2, V}
    delete!(dd.dict1, index)
    return dd
end

# Extend Base.delete! when the index is that of dict2
function Base.delete!(dd::DualDict{K1, K2, V},
                      index::K2)::DualDict{K1, K2, V} where {K1, K2, V}
    delete!(dd.dict2, index)
    return dd
end

## Helper functions for iterate
# iterate returned nothing
_process_iterate(id::Int, result::Nothing)::Nothing = result

# iterate returned a Tuple
function _process_iterate(id::Int, state::Tuple{Pair{K, V}, Int}
                          )::Tuple{Pair{K, V}, Tuple{Int, Int}} where {K, V}
    return state[1], (id, state[2])
end

# Extend Base.iterate (initial)
function Base.iterate(dd::DualDict{K1, K2, V}
                      )::Union{Nothing, Tuple{Pair{K1,V},Tuple{Int, Int}},
                               Tuple{Pair{K2,V},Tuple{Int, Int}}} where {K1, K2, V}
    if !isempty(dd.dict1)
        return _process_iterate(1, iterate(dd.dict1))
    else
        return _process_iterate(2, iterate(dd.dict2))
    end
end

# Extend Base.iterate (later)
function Base.iterate(dd::DualDict{K1, K2, V},
                      state::Tuple{Int, Int}
                      )::Union{Nothing, Tuple{Pair{K1,V},Tuple{Int, Int}},
                               Tuple{Pair{K2,V},Tuple{Int, Int}}} where {K1, K2, V}
    if state[1] == 1
        result = _process_iterate(1, iterate(dd.dict1, state[2]))
        return result !== nothing ? result : _process_iterate(2, iterate(dd.dict2))
    else
        return _process_iterate(2, iterate(dd.dict2, state[2]))
    end
end

# Extend Base.keys
function Base.keys(dd::DualDict)::Vector
    return [collect(keys(dd.dict1))..., collect(keys(dd.dict2))...]
end

# Extend Base.values
function Base.values(dd::DualDict{K1, K2, V})::Vector{V} where {K1, K2, V}
    return [values(dd.dict1)..., values(dd.dict2)...]
end

# Make the show string
function show_string(dd::DualDict)::String
    num_pairs = length(dd)
    if num_pairs == 1
        str = string(typeof(dd), " with 1 entry:")
    else
        str = string(typeof(dd), " with ", num_pairs, " entries:")
    end
    for pair in dd
        str *= string("\n  ", pair)
    end
    return str
end

# Show DualDicts in REPLMode
function Base.show(io::IO, dd::DualDict)
    print(io, show_string(dd))
    return
end

# Show DualDicts in IJuliaMode
function Base.show(io::IO, ::MIME"text/latex", dd::DualDict)
    print(io, show_string(dd))
end
