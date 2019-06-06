"""
    @infinite_variable(args...)
A wrapper macro for the `JuMP.@variable` macro that behaves the same except that
it defines variables of type `InfiniteVariableRef`.
"""
macro infinite_variable(model, args...)
    code = quote
        @assert isa($model, InfiniteModel)
        JuMP.@variable($model, ($(args...)), variable_type = Infinite)
    end
    return esc(code)
end

"""
    @point_variable(args...)
A wrapper macro for the `JuMP.@variable` macro that behaves the same except that
it defines variables of type `PointVariableRef`.
"""
macro point_variable(model, args...)
    code = quote
        @assert isa($model, InfiniteModel)
        JuMP.@variable($model, ($(args...)), variable_type = Point)
    end
    return esc(code)
end

"""
    @global_variable(args...)
A wrapper macro for the `JuMP.@variable` macro that behaves the same except that
it defines variables of type `GlobalVariableRef`.
"""
macro global_variable(model, args...)
    code = quote
        @assert isa($model, InfiniteModel)
        JuMP.@variable($model, ($(args...)), variable_type = Global)
    end
    return esc(code)
end
