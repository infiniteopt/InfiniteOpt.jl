## A template for defining a new transformation backend

# NOTE THAT USING `JuMPBackend` GREATLY SIMPLIFIES THIS (SEE DOCS)

# Define a mutable struct for storing infinite model to backend mappings
# plus other needed information
mutable struct NewReformData
    # Variable mapping (REPLACE AND REFORMAT BELOW AS NEEDED)
    infvar_mappings::Dict{GeneralVariableRef, Vector{JuMP.VariableRef}}
    finvar_mappings::Dict{GeneralVariableRef, JuMP.VariableRef}

    # Map other variable info (REPLACE AND REFORMAT BELOW AS NEEDED)
    infvar_to_supports::Dict{GeneralVariableRef, Vector{<:Tuple}}

    # Measure mappings and other info (ONLY NEEDED TO ENABLE MEASURE QUERIES)
    meas_mappings::Dict{GeneralVariableRef, Vector{JuMP.VariableRef}}
    meas_to_supports::Dict{GeneralVariableRef, Vector{<:Tuple}}

    # Constraint mapping (REPLACE AND REFORMAT BELOW AS NEEDED)
    constr_mappings::Dict{InfOptConstraintRef, Vector{JuMP.ConstraintRef}}

    # Map other constraint info (REPLACE AND REFORMAT BELOW AS NEEDED)
    constr_to_supports::Dict{InfOptConstraintRef, Vector{<:Tuple}}

    # ADD OTHER MAPPING/MODEL INFORMATION HERE

    # default constructor
    function NewReformData() # REFORMAT BELOW IN ACCORDANCE WITH ABOVE ATTRIBUTES
        return new(Dict{GeneralVariableRef, Vector{JuMP.VariableRef}}(),
                   Dict{GeneralVariableRef, JuMP.VariableRef}(),
                   Dict{GeneralVariableRef, Vector{Tuple}}(),
                   Dict{GeneralVariableRef, Vector{JuMP.VariableRef}}(),
                   Dict{GeneralVariableRef, Vector{Tuple}}(),
                   Dict{InfOptConstraintRef, Vector{JuMP.ConstraintRef}}(),
                   Dict{InfOptConstraintRef, Vector{Tuple}}()
                   )
    end
end

# Extend Base.empty!
function Base.empty!(data::NewReformData)
    # REFORMAT BELOW AS APPROPRIATE
    empty!(data.infvar_mappings)
    empty!(data.finvar_mappings)
    empty!(data.infvar_to_supports)
    empty!(data.meas_mappings)
    empty!(data.meas_to_supports)
    empty!(data.constr_mappings)
    empty!(data.constr_to_supports)
    return data
end

# Define the new backend (Note that backends based on JuMP should actually use JuMPBackend)
struct NewReformBackend <: AbstractTransformationBackend
    model::JuMP.Model
    data::NewReformData
end

# Make a constructor for the new backend
function NewReformBackend(args...; kwargs...) # ADD EXPLICT ARGS AS NEEDED
    # CREATE THE BACKEND AS NEEDED
    return NewReformBackend(JuMP.Model(args...; kwargs...), NewReformData())
end

# Extend Base.empty! for the backend (not needed for `JuMPBackend`s)
function Base.empty!(backend::NewReformBackend)
    # REPLACE WITH ACTUAL
    empty!(backend.model)
    empty!(backend.data)
    return backend
end

# Extend transformation_model and transformation_data (not needed for `JuMPBackend`s)
InfiniteOpt.transformation_model(backend::NewReformBackend) = backend.model
InfiniteOpt.transformation_data(backend::NewReformBackend) = backend.data

# Extend JuMP.[get/set]_attribute such to pass and retrieve settings for the backend
# Not necessary for `JuMPBackend`s
function JuMP.get_attribute(backend::NewReformBackend, attr) 
    # REPLACE WITH ACTUAL
    return JuMP.get_attribute(backend.model, attr)
end
function JuMP.set_attribute(backend::NewReformBackend, attr, val) 
    # REPLACE WITH ACTUAL
    JuMP.set_attribute(backend.model, attr, val)
    return
end

# Extend JuMP.show_backend_summary (optional for better printing)
function JuMP.show_backend_summary(io::IO, model::InfiniteModel, backend::NewReformBackend)
    # ADD DESIRED INFORMATION (BE SURE TO INDENT BY 2 SPACES)
    println(io, "  Backend type: NewReformBackend")
    println(io, "  Some useful info: 42")
    return
end

# Extend build_transformation_backend!
function InfiniteOpt.build_transformation_backend!(
    model::InfiniteModel,
    backend::NewReformBackend;
    my_kwarg::Bool = true # ADD KEYWORD ARGUMENTS AS NEEDED
    )
    # clear backend for a build
    empty!(backend)
    backend.model.operator_counter = 0

    # load in user defined nonlinear operators
    # THIS HELPER FUNCTION ONLY IS FOR BACKENDS THAT USE JuMP
    add_operators_to_jump(backend.model, model)

    # IT MAY BE USEFUL TO CALL `expand_all_measures!` TO HANDLE MEASURES FIRST
    # otherwise can extend `add_measure_variable` and `delete_semi_infinite_variable` to
    # expand in place without modifying the infinite model

    # IT MAY BE USEFUL TO CALL `evaluate_all_derivatives!` TO HANDLE DERIVATIVES FIRST
    # otherwise can use `evaluate_derivative` or the combo of `derivative_expr_data` & `make_indexed_derivative_expr`

    # REPLACE BELOW WITH OPERATIONS TO BUILD A `NewReformBackend` BASED ON `backend`
    # these lines just generate artificial data, see `TranscriptionOpt` for a thorough example of implementation
    data = backend.data
    reform_model = backend.model
    for vref in all_variables(model)
        if index(vref) isa Union{InfiniteVariableIndex, SemiInfiniteVariableIndex}
            data.infvar_mappings[vref] = [@variable(reform_model) for _ = 1:2]
            data.infvar_to_supports[vref] = [(0.,), (1.,)]
        else
            data.finvar_mappings[vref] = @variable(reform_model)
        end
    end
    # TODO add derivatives
    for mref in all_measures(model)
        data.meas_mappings[mref] = [@variable(reform_model) for _ = 1:2]
        data.meas_to_supports[mref] = [(-1.,), (-2.,)]
    end
    for cref in all_constraints(model)
        data.constr_mappings[cref] = [@constraint(reform_model, @variable(reform_model) >= 0) for _ = 1:2]
        data.constr_to_supports[cref] = [(2.,), (3.,)]
    end
    set_objective(reform_model, objective_sense(model), @variable(reform_model))
    return
end

# Extend JuMP.optimize! (not needed for `JuMPBackend`s)
function JuMP.optimize!(backend::NewReformBackend)
    # DO WANT NEEDS TO BE DONE TO SOLVE THE TRANSFORMED MODEL
    return JuMP.optimize!(backend.model)
end

# To the extend desired, extend any or all of the JuMP API below
# Not needed for `JuMPBackend`s
for func in (:set_silent, :unset_silent, :bridge_constraints, 
             :unset_time_limit_sec, :time_limit_sec, :solver_name, :backend,
             :mode, :unsafe_backend, :compute_conflict!, :copy_conflict,
             :set_string_names_on_creation)
    @eval begin
        # EXTEND FUNCTION AS APPROPRIATE
        function JuMP.$func(backend::NewReformBackend)
            return JuMP.$func(backend.model)
        end
    end
end
for func in (:set_time_limit_sec, :set_string_names_on_creation, :add_bridge)
    @eval begin
        # EXTEND FUNCTION AS APPROPRIATE
        function JuMP.$func(backend::NewReformBackend, value)
            return JuMP.$func(backend.model, value)
        end
    end
end
function JuMP.print_active_bridges(io::IO, backend::NewReformBackend, args...)
    return JuMP.print_active_bridges(io, backend.model, args...)
end
function JuMP.print_bridge_graph(io::IO, backend::NewReformBackend)
    return JuMP.print_bridge_graph(io, backend.model)
end
function JuMP.set_optimizer(
    backend::NewReformBackend,
    optimizer_constructor;
    add_bridges::Bool = true # ADD KWARGS AS DESIRED
    )
    JuMP.set_optimizer(backend.model, optimizer_constructor, add_bridges = add_bridges)
    return
end

# Extend transformation_variable if appropriate to enable variable related queries
function InfiniteOpt.transformation_variable(
    vref::GeneralVariableRef,
    backend::NewReformBackend;
    my_kwarg::Bool = true # ADD KEY ARGS AS NEEDED
    )
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE OPTIMIZER MODEL VARIABLE(S)
    vindex = index(vref)
    if vindex isa Union{InfiniteVariableIndex, SemiInfiniteVariableIndex, DerivativeIndex}
        map_dict = backend.data.infvar_mappings
    elseif vindex isa MeasureIndex
        map_dict = backend.data.meas_mappings
    else
        map_dict = backend.data.finvar_mappings
    end
    haskey(map_dict, vref) || error("Variable $vref not used in the backend.")
    return map_dict[vref]
end

# Extend transformation_expression if appropriate to enable expression related queries
function InfiniteOpt.transformation_expression(
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, JuMP.GenericNonlinearExpr}, # POSSIBLY BREAK THESE UP INTO 3 SEPARATE FUNCTIONS
    backend::NewReformBackend;
    my_kwarg::Bool = true # ADD KEY ARGS AS NEEDED
    )
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO REFORMULATED EXPRESSIONS
    reform_expr = zero(AffExpr)
    return reform_expr
end

# Extend transformation_constraint if appropriate to enable constraint related queries
function InfiniteOpt.transformation_constraint(
    cref::InfOptConstraintRef,
    backend::NewReformBackend;
    my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE OPTIMIZER MODEL CONSTRAINT(S)
    map_dict = backend.data.constr_mappings
    haskey(map_dict, cref) || error("Constraint $cref not used in the backend.")
    return map_dict[cref]
end

# If appropriate extend variable_supports (enables support queries of infinite variables)
function InfiniteOpt.variable_supports(
    vref::Union{InfiniteVariableRef, SemiInfiniteVariableRef, DerivativeRef},
    backend::NewReformBackend;
    my_kwarg::Bool = true # ADD KEY ARGS AS NEEDED
    )
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE INFINITE VARIABLE SUPPORT VALUES
    map_dict = backend.data.infvar_to_supports
    gvref = InfiniteOpt.GeneralVariableRef(owner_model(vref), index(vref))
    haskey(map_dict, gvref) || error("Variable $gvref not used in the backend.")
    return map_dict[gvref]
end

# If appropriate extend variable_supports for measures (enables support queries of measures)
function InfiniteOpt.variable_supports(
    mref::MeasureRef,
    backend::NewReformBackend;
    my_kwarg::Bool = true # ADD KEY ARGS AS NEEDED
    )
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE INFINITE VARIABLE SUPPORT VALUES
    map_dict = backend.data.meas_to_supports
    gvref = InfiniteOpt.GeneralVariableRef(owner_model(mref), index(mref))
    haskey(map_dict, gvref) || error("Variable $gvref not used in the backend.")
    return map_dict[gvref]
end

# If appropriate extend expression_supports (enables support queries of expressions)
function InfiniteOpt.expression_supports(
    expr::Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr, JuMP.GenericNonlinearExpr},
    backend::NewReformBackend;
    my_kwarg::Bool = true # ADD KEY ARGS AS NEEDED
    )
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO SUPPORT(S) OF THE EXPRESSION(S)
    supps = [(-42.,), (1.,)]
    return supps
end

# If appropriate extend constraint_supports (enables support queries of constraints)
function InfiniteOpt.constraint_supports(
    cref::InfOptConstraintRef,
    backend::NewReformBackend;
    my_kwarg::Bool = true) # ADD KEY ARGS AS NEEDED
    # REPLACE BELOW WITH ACTUAL CORRESPONDENCE TO THE INFINITE VARIABLE SUPPORT VALUES
    map_dict = backend.data.constr_to_supports
    haskey(map_dict, cref) || error("Constraint $cref does not have supports in the backend.")
    return length(map_dict[cref]) == 1 ? first(map_dict[cref]) : map_dict[cref]
end

# As desired and as appropriate extend any or all of the following JuMP query API
# Not required for `JuMPBackend`S
for func in (:termination_status, :raw_status, :solve_time, :simplex_iterations,
             :barrier_iterations, :node_count, :objective_bound, :relative_gap,
             :result_count)
    @eval begin
        # EXTEND AS NEEDED
        function JuMP.$func(backend::NewReformBackend)
            return JuMP.$func(backend.model)
        end
    end
end
# Simple result dependent model queries
for func in (:primal_status, :dual_status, :has_values, :has_duals, 
             :objective_value, :dual_objective_value)
    @eval begin 
        # EXTEND AS NEEDED
        function JuMP.$func(backend::NewReformBackend; kwargs...)
            return JuMP.$func(backend.model; kwargs...)
        end
    end
end

# Extend map_value to query optimal values
# Not needed for `JuMPBackend`s that implement `transformation_[variable/expression/constraint]`
function InfiniteOpt.map_value(
    vref::GeneralVariableRef,
    backend::NewReformBackend;
    result::Int = 1, # ADD ANY DESIRED KEY ARGS
    kwargs... # EXTRA CAN BE PASSED ON TO THE MAPPING
    )
    # REPLACE WITH ACTUAL MAPPING
    opt_vref = transformation_variable(vref, backend; kwargs...)
    if opt_vref isa AbstractArray
        return map(v -> JuMP.value(v, result = result), opt_vref)
    else
        return JuMP.value(opt_vref, result = result)
    end
end
function InfiniteOpt.map_value(
    expr::JuMP.AbstractJuMPScalar,
    backend::NewReformBackend;
    result::Int = 1, # ADD ANY DESIRED KEY ARGS
    kwargs... # EXTRA CAN BE PASSED ON TO THE MAPPING
    )
    # REPLACE WITH ACTUAL MAPPING
    opt_expr = transformation_expression(expr, backend; kwargs...)
    if opt_expr isa AbstractArray
        return map(v -> JuMP.value(v, result = result), opt_expr)
    else
        return JuMP.value(opt_expr, result = result)
    end
end
function InfiniteOpt.map_value(
    cref::InfOptConstraintRef,
    backend::NewReformBackend;
    result::Int = 1, # ADD ANY DESIRED KEY ARGS
    kwargs... # EXTRA CAN BE PASSED ON TO THE MAPPING
    )
    # REPLACE WITH ACTUAL MAPPING
    opt_cref = transformation_constraint(cref, backend; kwargs...)
    if opt_cref isa AbstractArray
        return map(c -> JuMP.value(c, result = result), opt_cref)
    else
        return JuMP.value(opt_cref, result = result)
    end
end

# IF NEEDED, EXTEND `InfiniteOpt.map_infinite_parameter_value` (defaults to using `supports(pref)`)

# Extend `map_dual` to be able to query duals
# Not needed for `JuMPBackend`s that implement `transformation_constraint`
function InfiniteOpt.map_dual(
    cref::InfOptConstraintRef,
    backend::NewReformBackend;
    result::Int = 1, # ADD ANY DESIRED KEY ARGS
    kwargs... # EXTRA CAN BE PASSED ON TO THE MAPPING
    )
    # REPLACE WITH ACTUAL MAPPING
    opt_cref = transformation_constraint(cref, backend; kwargs...)
    if opt_cref isa AbstractArray
        return map(c -> JuMP.dual(c, result = result), opt_cref)
    else
        return JuMP.dual(opt_cref, result = result)
    end
end

# IF APPROPRIATE, EXTEND THE FOLLOWING
# Not required for `JuMPBackend`s that use `transformation_[variable/constraint]`
#   `InfiniteOpt.map_reduced_cost`
#   `InfiniteOpt.map_optimizer_index`
#   `InfiniteOpt.map_shadow_price`
#   `JuMP.lp_sensitivity_report`
