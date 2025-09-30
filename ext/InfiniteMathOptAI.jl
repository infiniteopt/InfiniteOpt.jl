module InfiniteMathOptAI

import InfiniteOpt, JuMP, MathOptAI

# Extend MathOptAI to add variables correctly
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

# TODO extend MathOptAI.[get/set]_variable_bounds if/when function bounds are allowed in InfiniteOpt

end # end of module