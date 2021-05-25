import Random # used for example

# Define the new generative support information struct 
struct MyGenInfo <: InfiniteOpt.AbstractGenerativeInfo  # REPLACE WITH ACTUAL NAME
    attr::Int # REPLACE WITH DESIRED ATTRIBUTES
    # ADD MORE ATTRIBUTES AS NEEDED
end

# Make a unique support label (recommended)
struct MyGenLabel <: InternalLabel end

# Extend support_label 
function InfiniteOpt.support_label(info::MyGenInfo)::Type{MyGenLabel} # REPLACE WITH ACTUAL INFO TYPE
    return MyGenLabel # REPLACE WITH ACTUAL MAPPING
end

# Extend make_generative_supports 
function InfiniteOpt.make_generative_supports(
    info::MyGenInfo, 
    pref, 
    supps
    )::Vector{Float64} # REPLACE WITH ACTUAL INFO TYPE
    # REPLACE BELOW WITH ACTUAL CODE TO CREATE THE GENERATIVE SUPPORTS BASED ON THE EXISTING
    num_existing = length(supps)
    num_existing <= 1 && error("`$pref` doesn't have enough supports.")
    num_internal = info.attr
    gen_supps = Float64[]
    for i = 1:num_existing-1 
        lb = supps[i]
        ub = supps[i+1]
        append!(gen_supps, rand(num_internal) * (ub - lb) .+ lb)
    end
    return gen_supps
end
