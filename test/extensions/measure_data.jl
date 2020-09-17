# Define the new measure data type with <: AbstractMeasureData
struct NewMeasureData <: AbstractMeasureData
    attr1::String # REPLACE WITH ACTUAL ATTRIBUTE
    attr2::DiscreteMeasureData # REPLACE WITH ACTUAL ATTRIBUTE
    # ADD MORE ATTRIBUTES AS NEEDED
    # constructor
    function NewMeasureData(attr1::String, attr2::DiscreteMeasureData)
        # INSERT CHECKS AND OTHER CONSTRUCTOR METHODS HERE
        return new(attr1, attr2) # REPLACE WITH ACTUAL ATTRIBUTES
    end
end

# Extend parameter_refs to return the parameter(s) being measured by a measure using NewMeasureData
function InfiniteOpt.parameter_refs(data::NewMeasureData)::Union{GeneralVariableRef, AbstractArray{<:GeneralVariableRef}}
    return data.attr2.parameter_refs # REPLACE WITH ACTUAL PARAMETER LOCATION
end

# Extend expand_measure to return the finite reformulation of a measure using NewMeasureData
function InfiniteOpt.expand_measure(expr::JuMP.AbstractJuMPScalar,
                                    data::NewMeasureData,
                                    write_model::JuMP.AbstractModel
                                    )::JuMP.AbstractJuMPScalar
    # INSERT APPROPRIATE METHODS HERE
    # USING make_point_variable_ref AND make_reduced_variable_ref MAY BE USEFUL
    return expand_measure(expr, data.attr2, write_model) # REPLACE ACTUAL RESULT
end

# Extend add_supports_to_parameters to add supports to parameters according to the NewMeasureData
# This is only optional if the measure is always analytic
function InfiniteOpt.add_supports_to_parameters(data::NewMeasureData)
    return add_supports_to_parameters(data.attr2) # REPLACE WITH ACTUAL SUPPORT ADDITION STEPS
end

# Extend supports to return any infinite parameter supports employed by NewMeasureData
# This is only optional if the new abstraction doesn't use supports at all
function InfiniteOpt.supports(data::NewMeasureData)::Vector
    return data.attr2.supports # REPLACE WITH ACTUAL LOCATION
end

# Extend measure_data_in_hold_bounds to determine if NewMeasureData is in the
# domain of hold variable bounds. (Enables hold variable error checking)
function InfiniteOpt.measure_data_in_hold_bounds(data::NewMeasureData,
                                                 bounds::ParameterBounds)::Bool
    # INSERT ACTUAL CHECK HERE
    in_bounds = measure_data_in_hold_bounds(data.attr2, bounds) # REPLACE WITH ACTUAL RESULT
    return in_bounds
end

# Extend coefficients to return the coefficients stored in NewMeasureData if appropriate
# This is optional (returns empty vector otherwise)
function InfiniteOpt.coefficients(data::NewMeasureData)::Vector{Float64}
    return data.attr2.coefficients # REPLACE WITH ACTUAL LOCATION
end

# Extend weight_function to return the weight function stored in NewMeasureData if appropriate
# This is optional (returns default_weight otherwise)
function InfiniteOpt.weight_function(data::NewMeasureData)::Function
    return data.attr2.weight_function # REPLACE WITH ACTUAL LOCATION
end

# Make a convenient measure constructor function for our new measure type
# This should employ measure(expr, data)
function new_measure(expr::JuMP.AbstractJuMPScalar, param::GeneralVariableRef,
                     lb::Number, ub::Number; name::String = "NewMeas",
                     num_supports::Int = 10)::GeneralVariableRef # REPLACE ARGS WITH ACTUAL DESIRED
    # INSERT RELAVENT CHECKS AND OPERATIONS HERE
    # REPLACE BELOW WITH ACTUAL CONSTRUCTION
    attr2 = DiscreteMeasureData(param, ones(num_supports), [lb + (ub - lb) / (num_supports - 1) * i for i in 1:num_supports]) # just an example
    data = NewMeasureData(name, attr2) # REPLACE WITH ACTUAL
    # built the measure using the built-in constructor
    return measure(expr, data)
end
