################################################################################
#                           ATTRIBUTE CACHE METHODS
################################################################################
## Extend Base.[get/set]index
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
        function Base.setindex(cache::TransformAttrCache, idx::$I, attr::$A, value)
            cache.$(name)[idx, attr] = value
            return
        end
    end
end

# Model attributes
function Base.getindex(cache::TransformAttrCache, attr::ModelAttr)
    return cache.model[attr]
end
function Base.setindex(cache::TransformAttrCache, attr::ModelAttr, value)
    cache.model[attr] = value
    return
end
