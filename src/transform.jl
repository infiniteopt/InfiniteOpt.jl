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
                error("Expected a value of type `$(attribute_value_type(attr))` for ", 
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
        error("Expected a value of type `$(attribute_value_type(attr))` for ", 
              "the attribute `$(attr)`, but got `$(value)` of type ",
              "`$(typeof(value))`.")
    end
    model.transform_attrs[attr] = value
    return
end

################################################################################
#                              KEYWORD INPUT API
################################################################################
"""
    register_transform_keyword(kw::Symbol, attr::AbstractTransformAttr, 
                               [type = attribute_value_type(attr)])::Nothing

Register a new [`AbstractTransformAttr`](@ref) to be specified via keyword 
argument to the appropriate creation macro (e.g., `@infinite_parameter` for 
[`InfiniteParameterAttr`](@ref)s). Using `type` indicates what type of input 
of input is accepted for `kw`. Note that `type` need not equal 
`attribute_value_type(attr)` if [`process_transform_value`](@ref) is extended 
for that value type. This is intended for advanced users who are implementing 
their own transformation backend. Errors if `kw` is already in use.
"""
function register_transform_keyword end

# Helper function for keyword registration
function _add_transform_keyword(dict, kw, attr, type)
    if haskey(dict, kw)
        error("Cannot register $kw as a keyword argument, since it is already ",
              "in use.")
    end
    if isnothing(type)
        dict[kw] = (attr, attribute_value_type(attr))
    else
        dict[kw] = (attr, type)
    end
    return
end

"""
    process_transform_value(_error::Function, attr, value, object)

Process a transformation backend attribute `attr` with a raw `value` that will 
be associated with `object` once it is added to the model. This is intended 
to handle raw input `value`s collected via keyword arguments when `object` is 
created (i.e., transformation keywords that have been registered via 
[`register_transform_keyword`](@ref)). This provides extra flexibility allow 
more convenient input for users. By default `value` is simply returned. This 
should be extended for registered transformation attribute types that wish to 
process input that doesn't readily match [`attribute_value_type`](@ref). Note 
that this is intended as an advanced function for developers of transformation 
backends.
"""
function process_transform_value(_error::Function, attr, value, object)
    return value
end

# Helper function for checking the transformation keywords in build functions
function _process_transform_kwargs(_error, dict, kwargs, obj)
    for (k, v) in kwargs
        if !haskey(dict, k)
            _error("Unrecognized keyword argument `$k`.")
        elseif !isa(v, dict[k][2])
            _error("The `$k` keyword argument should be of type ",
                   "`$(dict[k][2])` not of type `$(typeof(v))`.")
        end
    end
    if isempty(kwargs)
        return obj
    else
        processed_kwargs = Dict()
        sizehint!(processed_kwargs, length(kwargs))
        for (k, v) in kwargs
            attr = dict[k][1]
            if haskey(processed_kwargs, attr)
                kw_inds = findall(kw -> dict[kw][1] == attr, keys(kwargs))
                kws = keys(kwargs)[kw_inds]
                _error("The following keyword arguments cannot be given ",
                       "simultaneously: '", join(kws[1:end-1], "', "), "' and '",
                       kws[end], "'.")
            end
            processed_kwargs[attr] = process_transform_value(_error, attr, v, obj)
        end
        return ObjectWithAttributes(obj, processed_kwargs)
    end
end

################################################################################
#                         DYNAMIC ATTRIBUTE UPDATE API
################################################################################
# Create storage container for transform attributes to be updated incrementally
const _IncrementalAttributes = Set{AbstractTransformAttr}()

"""
    register_attribute_to_update(attr::AbstractTransformAttr)::Nothing

Register a transformation backend attribute `attr` that can be updated 
incrementally when modeling objects are added to `InfiniteModel`s. This 
must be implemented in combination with [`update_attribute_on_creation`](@ref). 
Only individuals writing transformation backends should use this function.
"""
function register_attribute_to_update(attr::AbstractTransformAttr)
    push!(_IncrementalAttributes, attr)
    return
end

"""
    update_attribute_on_creation(model::InfiniteModel, attr::AbstractTransformAttr, obj)::Nothing

Update a particular transformation backend attribute `attr` when an `InfiniteOpt` 
modeling object `obj` is added to `model`. This enables transformation attributes 
to be built-up incrementally as the `InfiniteModel` is created. This method is 
only invoked for attributes that have been registered via 
[`register_attribute_to_update`](@ref). By default, nothing is updated. Those 
writing transformation backends should extend this based on particular 
attribute-object combinations where an update should be made. For instance, we 
could obtain the discrete supports with the point from a point variable when it 
is created (this is what `TranscriptionOpt` does). 

Note this should NOT be used to modify/create an attribute for `obj` itself, it 
should only be used to modify attributes associated with other modeling objects. 
Instead use [`register_transform_keyword`](@ref) if you what to create an attribute 
for `obj` on creation. 
"""
function update_attribute_on_creation(
    model::InfiniteModel, 
    attr::AbstractTransformAttr, 
    obj
    )
    return 
end

# Helper function for incrementally updating transform attributes when 
# adding modeling objects to the model
function _update_transform_attributes(model, obj)
    for attr in _IncrementalAttributes
        update_attribute_on_creation(model, attr, obj)
    end
    return
end
