################################################################################
#                           ATTRIBUTE CACHE METHODS
################################################################################
## Extend basic Base functions
# Non-model attributes
for (I, A, name) = ((:FiniteParameterIndex, :FiniteParameterAttr, :finite_parameters), 
          (:IndependentParameterIndex, :InfiniteParameterAttr, :independent_parameters),
          (:DependentParametersIndex, :InfiniteParameterAttr, :dependent_parameters),
          (:InfiniteVariableIndex, :VariableAttr, :infinite_variables),
          (:SemiInfiniteVariableIndex, :VariableAttr, :semi_infinite_variables),
          (:PointVariableIndex, :VariableAttr, :point_variables),
          (:DerivativeIndex, :DerivativeAttr, :derivatives),
          (:MeasureIndex, :MeasureAttr, :measures),
          (:InfOptConstraintIndex, :ConstraintAttr, :constraints)
        )
    @eval begin
        function Base.getindex(cache::TransformAttrCache, idx::$I, attr::$A)
            return cache.$(name)[idx, attr]
        end
        function Base.setindex!(cache::TransformAttrCache, value, idx::$I, attr::$A)
            return cache.$(name)[idx, attr] = value
        end
        function Base.haskey(cache::TransformAttrCache, key::Tuple{$I, $A})
            return haskey(cache.$(name), key)
        end
        function Base.get(cache::TransformAttrCache, key::Tuple{$I, $A}, default)
            return Base.get(cache.$(name), key, default)
        end
    end
end

# Model attributes
function Base.getindex(cache::TransformAttrCache, attr::ModelAttr)
    return cache.model[attr]
end
function Base.setindex!(cache::TransformAttrCache, value, attr::ModelAttr)
    return cache.model[attr] = value
end
function Base.haskey(cache::TransformAttrCache, attr::ModelAttr)
    return haskey(cache.model, attr)
end
function Base.get(cache::TransformAttrCache, attr::ModelAttr, default)
    return Base.get(cache.model, attr, default)
end

################################################################################
#                              BASIC ATTRIBUTE API
################################################################################
## Define the basic getters/setters
# Fallbacks
"""
    InfiniteOpt.get(model::InfiniteModel, [index], attribute)::AbstractTransformAttr

Retrieve the value of the transformation backend attribute `attribute` that is 
associated with object `index` in `model`. If `attribute` is a [`ModelAttr`](@ref) 
then the `index` argument is omitted. Errors if the `index` and `model` have 
no such attribute. This is intended for use by those writing transformation 
backends for `InfiniteModel`s. 
"""
function get(model::InfiniteModel, idx, attr)
    error("Objects with indices of type `$(typeof(idx))` are not compatible with ",
          "attributes of type `$(typeof(attr))`.")
end
function get(model::InfiniteModel, attr)
    error("`InfiniteModel`s are not compatible with ",
          "attributes of type `$(typeof(attr))`.")
end

"""
    attribute_value_type(attribute::AbstractTransformAttr)::DataType

Returns the type of a value that an `attribute` can accept. This should be 
extended for new [`AbstractTransformAttr`](@ref)s such that they can be 
checked. Defaults to `Any` such that [`InfiniteOpt.set`](@ref) does not 
check the value type of an attribute. 
"""
attribute_value_type(::AbstractTransformAttr) = Any

"""
    InfiniteOpt.set(model::InfiniteModel, [index], attribute, value)::Nothing

Set the transformation backend attribute `attribute` for object `index` in 
`model` to `value`. If `attribute` is a [`ModelAttr`](@ref) then the `index` 
argument is omitted. Errors if `value` is incompatible with `attribute`. This 
is intended for use by those writing transformation backends for `InfiniteModel`s. 
"""
function set(model::InfiniteModel, idx, attr, value)
    error("Objects with indices of type `$(typeof(idx))` are not compatible with ",
          "attributes of type `$(typeof(attr))`.")
end
function set(model::InfiniteModel, attr, value)
    error("`InfiniteModel`s are not compatible with ",
          "attributes of type `$(typeof(attr))`.")
end

# TODO enable `set` to optionally update the transformation backend incrementally

# Singular dispatch variables (and the dependent parameter & constraint getters)
for (I, A) = ((:FiniteParameterIndex, :FiniteParameterAttr), 
              (:IndependentParameterIndex, :InfiniteParameterAttr),
              (:DependentParametersIndex, :InfiniteParameterAttr),
              (:InfiniteVariableIndex, :VariableAttr),
              (:SemiInfiniteVariableIndex, :VariableAttr),
              (:PointVariableIndex, :VariableAttr),
              (:DerivativeIndex, :DerivativeAttr),
              (:MeasureIndex, :MeasureAttr),
              (:InfOptConstraintIndex, :ConstraintAttr))
    @eval begin
        if !($I in (DependentParametersIndex, InfOptConstraintIndex))
            function get(model::InfiniteModel, idx::$I, attr::$A)
                value = Base.get(model.transform_attrs, (idx, attr), nothing)
                if isnothing(value)
                    error("$(dispatch_variable_ref(model, idx)) does not have transform ",
                          "attribute `$attr`.")
                end
                return value
            end
        end
        function set(model::InfiniteModel, idx::$I, attr::$A, value)
            if !isa(value, attribute_value_type(attr))
                error("Expected a value of type `$(attribute_value_type(attr))` for", 
                      "the attribute `$(attr)`, but got `$(value)` of type ",
                      "`$(typeof(value))`.")
            end
            model.transform_attrs[idx, attr] = value
            return
        end
    end
end

# Additional getters for constraints, models, and dependent parameters
function get(
    model::InfiniteModel, 
    idx::InfOptConstraintIndex, 
    attr::ConstraintAttr
    )
    value = Base.get(model.transform_attrs, (idx, attr), nothing)
    if isnothing(value)
        error("$(InfOptConstraintRef(model, idx)) does not have transform ",
              "attribute `$attr`.")
    end
    return value
end
function get(
    model::InfiniteModel, 
    idx::DependentParametersIndex, 
    attr::InfiniteParameterAttr
    )
    value = Base.get(model.transform_attrs, (idx, attr), nothing)
    if isnothing(value)
        error("The dependent parameter group with index `$idx` does not ",
              "have transform attribute `$attr`.")
    end
    return value
end
function get(model::InfiniteModel, attr::ModelAttr)
    value = Base.get(model.transform_attrs, attr, nothing)
    if isnothing(value)
        error("The model does not have transform attribute `$attr`.")
    end
    return value
end

# Set for models 
function set(model::InfiniteModel, attr::ModelAttr, value)
    if !isa(value, attribute_value_type(attr))
        error("Expected a value of type `$(attribute_value_type(attr))` for", 
              "the attribute `$(attr)`, but got `$(value)` of type ",
              "`$(typeof(value))`.")
    end
    model.transform_attrs[attr] = value
    return
end
