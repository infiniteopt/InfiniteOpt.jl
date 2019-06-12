# TODO Clean up the parameter macro and add errors
"""
    @infinite_parameter(model, set, name_expr)
A macro for the defining parameters of type `ParameterRef`.
"""
macro infinite_parameter(model, name_expr, set)
    if isa(name_expr, Symbol)
        # easy case
        code = quote
            @assert isa($model, InfiniteModel)
            @assert typeof($set) <: AbstractInfiniteSet
            JuMP.@variable($model, ($(name_expr)), variable_type = Parameter, param_set = $set)
        end
        return esc(code)
    else
        if !JuMP.isexpr(name_expr, :ref)
            error("Syntax error: Expected $var to be of form var[...]")
        end

        model = esc(model)
        set = esc(set)

        variable = gensym()
        refcall, idxvars, idxsets, condition = JuMP._build_ref_sets(name_expr, variable)
        varname = JuMP._get_name(name_expr)
        escvarname = esc(varname)

        varstr = :(string($(string(varname)),"["))
        for idxvar in idxvars
            push!(varstr.args,:(string($(esc(idxvar)))))
            push!(varstr.args,",")
        end
        deleteat!(varstr.args,length(varstr.args))
        push!(varstr.args,"]")

        code = :( $(refcall) = JuMP.add_variable($model, InfOptParameter($set), $varstr) )
        looped = JuMP._get_looped_code(variable, code, condition, idxvars, idxsets, :ParameterRef, :Auto)
        return quote
            $looped
            $escvarname = $variable
        end
    end
end

# TODO Figure out a more robust type definition
"""
    @infinite_variable(model, set, args...)
A wrapper macro for the `JuMP.@variable` macro that behaves the same except that
it defines variables of type `InfiniteVariableRef`.
"""
macro infinite_variable(model, params, args...)
    code = quote
        @assert isa($model, InfiniteModel)
        @assert typeof($params) <: Union{ParameterRef, AbstractArray{<:ParameterRef}, Vector{<:AbstractArray{<:ParameterRef}}}
        JuMP.@variable($model, ($(args...)), variable_type = Infinite, parameters = $params)
    end
    return esc(code)
end

"""
    @point_variable(model, args...)
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
    @global_variable(model, args...)
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
