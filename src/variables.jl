"""
    FiniteVariableRef <: JuMP.AbstractVariableRef
An abstract type to define new finite variable types.
"""
abstract type FiniteVariableRef <: JuMP.AbstractVariableRef end

"""
GlobalVariableRef <: JuMP.AbstractVariableRef
A DataType for finite fixed variables (e.g., first stage variables,
steady-state variables).
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of variable in model.
"""
struct GlobalVariableRef <: FiniteVariableRef
    model::InfiniteModel # `model` owning the variable
    index::Int           # Index in `model.variables`
end

"""
PointVariableRef <: JuMP.AbstractVariableRef
A DataType for variables defined at a transcipted point (e.g., second stage
variable at a particular scenario, dynamic variable at a discretized time point).
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of variable in model.
"""
struct PointVariableRef <: FiniteVariableRef
    model::InfiniteModel # `model` owning the variable
    index::Int           # Index in `model.variables`
end

"""
InfiniteVariableRef <: JuMP.AbstractVariableRef
A DataType for untranscripted infinite dimensional variables (e.g., second stage
variables, dynamic variables).
**Fields**
- `model::InfiniteModel` Flexibility model.
- `index::Int` Index of variable in model.
"""
struct InfiniteVariableRef <: JuMP.AbstractVariableRef
    model::InfiniteModel # `model` owning the variable
    index::Int           # Index in `model.variables`
end

# Extend Base.copy for new variable types
Base.copy(v::Union{FiniteVariableRef, InfiniteVariableRef}) = v
Base.copy(v::InfiniteVariableRef, new_model::InfiniteModel) = InfiniteVariableRef(new_model, v.index)
Base.copy(v::GlobalVariableRef, new_model::InfiniteModel) = GlobalVariableRef(new_model, v.index)
Base.copy(v::PointVariableRef, new_model::InfiniteModel) = PointVariableRef(new_model, v.index)

# Extend other Base functions
Base.:(==)(v::Union{FiniteVariableRef, InfiniteVariableRef}, w::Union{FiniteVariableRef, InfiniteVariableRef}) = v.model === w.model && v.index == w.index
Base.broadcastable(v::Union{FiniteVariableRef, InfiniteVariableRef}) = Ref(v)

# Extend JuMP functions
JuMP.isequal_canonical(v::Union{FiniteVariableRef, InfiniteVariableRef}, w::Union{FiniteVariableRef, InfiniteVariableRef}) = v == w
JuMP.variable_type(m::InfiniteModel) = Union{FiniteVariableRef, InfiniteVariableRef}

# TODO setup args for defining the 3 new variables
function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, test_arg::Int; extra_kw_args...)
    for (kwarg, _) in extra_kw_args
        _error("Unrecognized keyword argument $kwarg")
    end
    return JuMP.ScalarVariable(info)
end

"""
    JuMP.add_variable(m::InfiniteModel, v::JuMP.AbstractVariable, name::String="")
Extend the `JuMP.add_variable` function to accomodate our new variable types.
"""
function JuMP.add_variable(m::InfiniteModel, v::JuMP.AbstractVariable, name::String="")
    m.next_var_index += 1
    vref = GlobalVariableRef(m, m.next_var_index)
    m.variables[vref.index] = v
    JuMP.set_name(vref, name)
    return vref
end

function JuMP.delete(model::InfiniteModel, vref::GlobalVariableRef)
    @assert JuMP.is_valid(model, vref)
    delete!(model.variables, vref.index)
    delete!(model.var_to_name, vref.index)
    return
end

function JuMP.is_valid(model::InfiniteModel, vref::GlobalVariableRef)
    return (model === vref.model &&
            vref.index in keys(model.variables))
end

JuMP.num_variables(m::InfiniteModel) = length(m.variables)

# Internal functions
variable_info(vref::GlobalVariableRef) = vref.model.variables[vref.index].info
function update_variable_info(vref::GlobalVariableRef, info::JuMP.VariableInfo)
    vref.model.variables[vref.index] = JuMP.ScalarVariable(info)
end

JuMP.has_lower_bound(vref::GlobalVariableRef) = variable_info(vref).has_lb

function JuMP.lower_bound(vref::GlobalVariableRef)::Float64
    @assert !JuMP.is_fixed(vref)
    return variable_info(vref).lower_bound
end

function JuMP.set_lower_bound(vref::GlobalVariableRef, lower)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(true, lower,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end

function JuMP.delete_lower_bound(vref::GlobalVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(false, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end

JuMP.has_upper_bound(vref::GlobalVariableRef) = variable_info(vref).has_ub

function JuMP.upper_bound(vref::GlobalVariableRef)::Float64
    @assert !JuMP.is_fixed(vref)
    return variable_info(vref).upper_bound
end

function JuMP.set_upper_bound(vref::GlobalVariableRef, upper)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           true, upper,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end

function JuMP.delete_upper_bound(vref::GlobalVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           false, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end

JuMP.is_fixed(vref::GlobalVariableRef) = variable_info(vref).has_fix

function JuMP.fix_value(vref::GlobalVariableRef)::Float64
    return variable_info(vref).fixed_value
end

function JuMP.fix(vref::GlobalVariableRef, value; force::Bool = false)
    info = variable_info(vref)
    if !force && (info.has_lb || info.has_ub)
        error("Unable to fix $(vref) to $(value) because it has existing bounds.")
    end
    update_variable_info(vref, JuMP.VariableInfo(
        false, info.lower_bound, false, info.upper_bound, true, value,
        info.has_start, info.start, info.binary, info.integer)
    )
    return
end

function JuMP.unfix(vref::GlobalVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           false, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, info.integer))
end

function JuMP.start_value(vref::GlobalVariableRef)::Union{Nothing, Float64}
    return variable_info(vref).start
end

function JuMP.set_start_value(vref::GlobalVariableRef, start)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           true, start,
                                           info.binary, info.integer))
end

JuMP.is_binary(vref::GlobalVariableRef) = variable_info(vref).binary

function JuMP.set_binary(vref::GlobalVariableRef)
    @assert !JuMP.is_integer(vref)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           true, info.integer))
end

function JuMP.unset_binary(vref::GlobalVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           false, info.integer))
end

JuMP.is_integer(vref::GlobalVariableRef) = variable_info(vref).integer
function JuMP.set_integer(vref::GlobalVariableRef)
    @assert !JuMP.is_binary(vref)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, true))
end

function JuMP.unset_integer(vref::GlobalVariableRef)
    info = variable_info(vref)
    update_variable_info(vref,
                         JuMP.VariableInfo(info.has_lb, info.lower_bound,
                                           info.has_ub, info.upper_bound,
                                           info.has_fix, info.fixed_value,
                                           info.has_start, info.start,
                                           info.binary, false))
end

JuMP.name(vref::GlobalVariableRef) = vref.model.var_to_name[vref.index]

function JuMP.set_name(vref::GlobalVariableRef, name::String)
    vref.model.var_to_name[vref.index] = name
    vref.model.name_to_var = nothing
end

function JuMP.variable_by_name(model::InfiniteModel, name::String)
    if model.name_to_var === nothing
        # Inspired from MOI/src/Utilities/model.jl
        model.name_to_var = Dict{String, Int}()
        for (var, var_name) in model.var_to_name
            if haskey(model.name_to_var, var_name)
                # -1 is a special value that means this string does not map to
                # a unique variable name.
                model.name_to_var[var_name] = -1
            else
                model.name_to_var[var_name] = var
            end
        end
    end
    index = get(model.name_to_var, name, nothing)
    if index isa Nothing
        return nothing
    elseif index == -1
        error("Multiple variables have the name $name.")
    else
        return GlobalVariableRef(model, index)
    end
end
