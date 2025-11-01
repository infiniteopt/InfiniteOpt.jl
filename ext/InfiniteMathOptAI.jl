module InfiniteMathOptAI

import InfiniteOpt, JuMP, MathOptAI

"""
    function MathOptAI.add_variables(
        model::InfiniteModel, x::Vector{GeneralVariableRef},
        n::Int,
        base_name::String,
    )::Vector{GeneralVariableRef}

Extend `MathOptAI.add_variables` to properly support infinite variables (i.e., ensure
the output variables of a predictor have the necessary infinite parameter dependencies).
This method should not be directly used by users, it is used to enable the use of 
[`MathOptAI.add_predictor`](https://lanl-ansi.github.io/MathOptAI.jl/stable/api/#add_predictor).
"""
function MathOptAI.add_variables(
    model::InfiniteOpt.InfiniteModel,
    x::Vector{InfiniteOpt.GeneralVariableRef},
    n::Int,
    base_name::String,
)
    inds = InfiniteOpt.parameter_group_int_indices(x)
    if isempty(inds)
        return JuMP.@variable(model, [1:n], base_name = base_name)
    end
    params = InfiniteOpt.parameter_refs(model)[inds]
    return JuMP.@variable(
        model,
        [1:n],
        base_name = base_name,
        variable_type = InfiniteOpt.Infinite(params...),
    )
end

# TODO implement these to fully support InfiniteOpt variables
MathOptAI.get_variable_bounds(::InfiniteOpt.GeneralVariableRef) = missing, missing
MathOptAI.get_variable_start(::InfiniteOpt.GeneralVariableRef) = missing

end # end of module