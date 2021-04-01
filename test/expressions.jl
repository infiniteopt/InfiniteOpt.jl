# Test infinite parameter function methods 
@testset "Infinite Parameter Function Methods" begin 
    # setup the needed info 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1])
    func = InfiniteParameterFunction(sin, IC.VectorTuple(t), [1], [1])
    object = ParameterFunctionData(func, "test")
    idx = ParameterFunctionIndex(1)
    fref = ParameterFunctionRef(m, idx)
    gvref = GeneralVariableRef(m, 1, ParameterFunctionIndex)
    bad_vref = ParameterFunctionRef(m, ParameterFunctionIndex(-1))
    # JuMP.owner_model
    @testset "JuMP.owner_model" begin
        @test owner_model(fref) === m
        @test owner_model(gvref) === m
    end
    # JuMP.index
    @testset "JuMP.index" begin
        @test index(fref) == idx
        @test index(gvref) == idx
    end
    # dispatch_variable_ref
    @testset "dispatch_variable_ref" begin
        @test dispatch_variable_ref(m, idx) == fref
        @test dispatch_variable_ref(gvref) == fref
    end
    # _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, object) == idx
    end
    # _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(fref) === m.param_functions
        @test InfiniteOpt._data_dictionary(gvref) === m.param_functions
    end
    # JuMP.is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, fref)
        @test is_valid(m, gvref)
    end
    # _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(fref) === object
        @test InfiniteOpt._data_object(gvref) === object
        @test_throws ErrorException InfiniteOpt._data_object(bad_vref)
    end
    # _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(fref) === func
        @test InfiniteOpt._core_variable_object(gvref) === func
    end
    # _object_numbers
    @testset "_object_numbers" begin
        @test InfiniteOpt._object_numbers(fref) == [1]
    end
    # _parameter_numbers
    @testset "_parameter_numbers" begin
        @test InfiniteOpt._parameter_numbers(fref) == [1]
    end
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(fref) == "test"
        @test name(gvref) == "test"
        @test name(bad_vref) == ""
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(fref, "new"), Nothing)
        @test name(fref) == "new"
        @test isa(set_name(gvref, "new2"), Nothing)
        @test name(fref) == "new2"
        @test_throws ErrorException set_name(bad_vref, "")
    end
    # _measure_dependencies
    @testset "_measure_dependencies" begin
        @test InfiniteOpt._measure_dependencies(fref) == MeasureIndex[]
        @test InfiniteOpt._measure_dependencies(gvref) == MeasureIndex[]
    end
    # _constraint_dependencies
    @testset "_constraint_dependencies" begin
        @test InfiniteOpt._constraint_dependencies(fref) == ConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(gvref) == ConstraintIndex[]
    end
    # _semi_infinite_variable_dependencies
    @testset "_semi_infinite_variable_dependencies" begin
        @test InfiniteOpt._semi_infinite_variable_dependencies(fref) == SemiInfiniteVariableIndex[]
        @test InfiniteOpt._semi_infinite_variable_dependencies(gvref) == SemiInfiniteVariableIndex[]
    end
    # _derivative_dependencies
    @testset "_derivative_dependencies" begin
        @test InfiniteOpt._derivative_dependencies(fref) == DerivativeIndex[]
        @test InfiniteOpt._derivative_dependencies(gvref) == DerivativeIndex[]
    end
    # raw_parameter_refs
    @testset "raw_parameter_refs" begin
        @test raw_parameter_refs(fref) == IC.VectorTuple(t)
        @test raw_parameter_refs(gvref) == IC.VectorTuple(t)
    end
    # parameter_refs
    @testset "parameter_refs" begin
        @test parameter_refs(fref) == (t, )
        @test parameter_refs(gvref) == (t, )
    end
    # parameter_list
    @testset "parameter_list" begin
        @test parameter_list(fref) == [t]
        @test parameter_list(gvref) == [t]
    end
    # raw_function
    @testset "raw_function" begin
        @test raw_function(fref) == sin
        @test raw_function(gvref) == sin
    end
    # test used_by_semi_infinite_variable
    @testset "used_by_semi_infinite_variable" begin
        @test !used_by_semi_infinite_variable(fref)
        push!(InfiniteOpt._semi_infinite_variable_dependencies(fref),
              SemiInfiniteVariableIndex(1))
        @test used_by_semi_infinite_variable(fref)
        empty!(InfiniteOpt._semi_infinite_variable_dependencies(fref))
    end
    # test used_by_derivative
    @testset "used_by_derivative" begin
        @test !used_by_derivative(fref)
        push!(InfiniteOpt._derivative_dependencies(fref), DerivativeIndex(1))
        @test used_by_derivative(fref)
        empty!(InfiniteOpt._derivative_dependencies(fref))
    end
    # test used_by_measure
    @testset "used_by_measure" begin
        @test !used_by_measure(fref)
        push!(InfiniteOpt._measure_dependencies(fref), MeasureIndex(1))
        @test used_by_measure(fref)
        empty!(InfiniteOpt._measure_dependencies(fref))
    end
    # test used_by_constraint
    @testset "used_by_constraint" begin
        @test !used_by_constraint(fref)
        push!(InfiniteOpt._constraint_dependencies(fref), ConstraintIndex(1))
        @test used_by_constraint(fref)
        empty!(InfiniteOpt._constraint_dependencies(fref))
    end
    # test is_used
    @testset "is_used" begin
        @test !is_used(fref)
        push!(InfiniteOpt._constraint_dependencies(fref), ConstraintIndex(1))
        @test is_used(fref)
        empty!(InfiniteOpt._constraint_dependencies(fref))
    end
    # test _update_param_var_mapping
    @testset "_update_param_var_mapping" begin 
        @test InfiniteOpt._update_param_var_mapping(fref, IC.VectorTuple(t)) isa Nothing 
        @test used_by_parameter_function(t)
        empty!(InfiniteOpt._parameter_function_dependencies(t))
    end
    # _delete_data_object
    @testset "_delete_data_object" begin
        @test InfiniteOpt._delete_data_object(fref) isa Nothing
        @test length(InfiniteOpt._data_dictionary(fref)) == 0
        @test !is_valid(m, fref)
    end
    # test build_parameter_function
    @testset "build_parameter_function" begin 
        # test errors  
        @test_throws ErrorException build_parameter_function(error, sin, x[1])
        @test_throws ErrorException build_parameter_function(error, sin, (t, x))
        # test normal  
        @test build_parameter_function(error, (a, b) -> 2, (t, x)) isa InfiniteParameterFunction 
    end
    # test add_parameter_function
    @testset "add_parameter_function" begin 
        # test does not belong 
        @infinite_parameter(InfiniteModel(), t2 in [0, 1])
        func = build_parameter_function(error, (a, b) -> 2, (t2, x))
        @test_throws VariableNotOwned{GeneralVariableRef} add_parameter_function(m, func)
        # test normal 
        func = build_parameter_function(error, (a, b) -> 2, (t, x))
        fref = GeneralVariableRef(m, 2, ParameterFunctionIndex)
        @test add_parameter_function(m, func, "test") == fref
        @test name(fref) == "test"
        @test parameter_refs(fref) == (t, x)
        @test used_by_parameter_function(t)
        # test default name 
        func = build_parameter_function(error, cos, t)
        fref = GeneralVariableRef(m, 3, ParameterFunctionIndex)
        @test add_parameter_function(m, func) == fref
        @test name(fref) == "cos"
    end
    # test parameter_function
    @testset "parameter_function" begin 
        @test parameter_function(sin, t) isa GeneralVariableRef
        @test parameter_function(sin, t, name = "name") isa GeneralVariableRef
    end
    # test making other objects 
    @testset "Other Objects" begin
        f = parameter_function((a,b) -> 2, (t, x))
        # test making semi_infinite variable
        d = Dict{Int, Float64}(1 => 0)
        @test add_variable(m, build_variable(error, f, d)) isa GeneralVariableRef
        # test making derivative 
        @test deriv(f, t) isa GeneralVariableRef
    end
end

# Test _all_function_variables
@testset "_all_function_variables" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    data = TestData(par, 0, 5)
    meas = Measure(hold, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    meas = InfiniteOpt._make_variable_ref(m, mindex)
    dinf = @deriv(inf, par)
    # test for variable reference
    @testset "Variable" begin
        @test InfiniteOpt._all_function_variables(par) == [par]
        @test InfiniteOpt._all_function_variables(inf) == [inf]
        @test InfiniteOpt._all_function_variables(meas) == [meas]
        @test InfiniteOpt._all_function_variables(dinf) == [dinf]
    end
    # test for GenericAffExpr
    @testset "AffExpr" begin
        # make expressions
        aff1 = meas + 2par + hold - dinf
        aff2 = zero(GenericAffExpr{Float64, GeneralVariableRef})
        # test expressions
        @test isempty(setdiff(InfiniteOpt._all_function_variables(aff1),
                              [meas, par, hold, dinf]))
        @test InfiniteOpt._all_function_variables(aff2) == GeneralVariableRef[]
    end
    # test for GenericQuadExpr
    @testset "QuadExpr" begin
        # make expressions
        quad1 = pt^2 + inf * pt - meas + 2par + hold - dinf
        quad2 = pt^2 + inf * pt
        quad3 = zero(GenericQuadExpr{Float64, GeneralVariableRef})
        # test expressions
        @test isempty(setdiff(InfiniteOpt._all_function_variables(quad1),
                      [meas, par, hold, dinf, pt, inf]))
        @test isempty(setdiff(InfiniteOpt._all_function_variables(quad2),
                              [pt, inf]))
        @test InfiniteOpt._all_function_variables(quad3) == GeneralVariableRef[]
    end
    # test backup
    @testset "Fallback" begin
        @variable(Model(), x)
        @test_throws ErrorException InfiniteOpt._all_function_variables(x)
    end
end

# Test comparisons
@testset "Comparisons" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, par2))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    # test AffExpr comparison
    @testset "Base.:(==) AffExpr" begin
        @test par + par2 + inf - 2 == par + (par2 + inf) - 2
        @test 0.25par- inf == 0.25par - inf
        @test inf + 3 - inf != par + inf
        @test inf + 3 - inf != inf + 3 - hold
    end
    # test QuadExpr comparison
    @testset "Base.:(==) QuadExpr" begin
        @test par * inf + par2 + inf - 2 == par * inf + (par2 + inf) - 2
        @test inf * inf - inf == inf * inf - inf
        @test inf * inf + 3 - inf != par * inf2 + inf
        @test par * par2 + inf * inf2 == par * par2 + inf * inf2
        @test par * par2 + inf * inf2 != par * par2
        @test par * par2 + 2 * inf * inf2 !=  par * par2 + inf * inf2
        @test par * inf + par2 + inf - 2 != par * inf + (par2 + inf) - 3
    end
end

# Test _object_numbers
@testset "_object_numbers" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, pars))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    var = build_variable(error, inf2, Dict{Int, Float64}(1 => 0.5), check = false)
    red = add_variable(m, var)
    # test for finite variable reference
    @testset "FiniteVariable" begin
        @test InfiniteOpt._object_numbers(pt) == []
        @test InfiniteOpt._object_numbers(hold) == []
    end
    # test for infinite variable reference
    @testset "InfiniteVariable" begin
        @test InfiniteOpt._object_numbers(inf) == [1]
    end
    # test for parameter reference
    @testset "Parameter" begin
        @test InfiniteOpt._object_numbers(par) == [1]
        @test InfiniteOpt._object_numbers(pars[1]) == [2]
    end
    # test for semi-infinite variable reference
    @testset "SemiInfiniteInfinite" begin
        @test InfiniteOpt._object_numbers(red) == [2]
    end
    # test for GenericAffExpr
    @testset "AffExpr" begin
        # make expressions
        aff1 = inf + inf2 + pt - 3
        aff2 = pt + hold - 2
        # test expressions
        @test sort!(InfiniteOpt._object_numbers(aff1)) == [1, 2]
        @test InfiniteOpt._object_numbers(aff2) == []
    end
    # test for GenericQuadExpr
    @testset "QuadExpr" begin
        # make expressions
        quad1 = inf * inf2 + inf + inf2 + pt - 3 - par
        quad2 = pt * pt + pt + hold - 2
        # test expressions
        @test sort!(InfiniteOpt._object_numbers(quad1)) == [1, 2]
        @test InfiniteOpt._object_numbers(quad2) == []
    end
end

# Test _parameter_numbers
@testset "_parameter_numbers" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= pars[1:2] <= 1)
    @infinite_variable(m, inf(par))
    @infinite_variable(m, inf2(par, pars))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    var = build_variable(error, inf2, Dict{Int, Float64}(1 => 0.5), check = false)
    red = add_variable(m, var)
    # test for finite variable reference
    @testset "FiniteVariable" begin
        @test InfiniteOpt._parameter_numbers(pt) == []
        @test InfiniteOpt._parameter_numbers(hold) == []
    end
    # test for infinite variable reference
    @testset "InfiniteVariable" begin
        @test InfiniteOpt._parameter_numbers(inf) == [1]
        @test InfiniteOpt._parameter_numbers(inf2) == [1, 2, 3]
    end
    # test for parameter reference
    @testset "Parameter" begin
        @test InfiniteOpt._parameter_numbers(par) == [1]
        @test InfiniteOpt._parameter_numbers(pars[2]) == [3]
    end
    # test for semi-infinite variable reference
    @testset "SemiInfiniteInfinite" begin
        @test InfiniteOpt._parameter_numbers(red) == [2, 3]
    end
    # test for GenericAffExpr
    @testset "AffExpr" begin
        # make expressions
        aff1 = inf + inf2 + pt - 3
        aff2 = pt + hold - 2
        # test expressions
        @test sort!(InfiniteOpt._parameter_numbers(aff1)) == [1, 2, 3]
        @test InfiniteOpt._parameter_numbers(aff2) == []
    end
    # test for GenericQuadExpr
    @testset "QuadExpr" begin
        # make expressions
        quad1 = inf * inf2 + inf + inf2 + pt - 3 - par
        quad2 = pt * pt + pt + hold - 2
        # test expressions
        @test sort!(InfiniteOpt._parameter_numbers(quad1)) == [1, 2, 3]
        @test InfiniteOpt._parameter_numbers(quad2) == []
    end
end

# Test _model_from_expr
@testset "_model_from_expr" begin
    # initialize model and references
    m = InfiniteModel()
    @hold_variable(m, hd)
    # test for variable reference
    @testset "Variable" begin
        @test InfiniteOpt._model_from_expr(hd) === m
    end
    # test for GenericAffExpr
    @testset "AffExpr" begin
        # make expressions
        aff1 = hd + 2
        aff2 = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        # test expressions
        @test InfiniteOpt._model_from_expr(aff1) === m
        @test InfiniteOpt._model_from_expr(aff2) isa Nothing
    end
    # test for GenericQuadExpr
    @testset "QuadExpr" begin
        # make expressions
        quad1 = hd * hd + hd + 1
        quad2 = zero(JuMP.GenericQuadExpr{Float64, GeneralVariableRef})
        quad3 = hd * hd + 1
        # test expressions
        @test InfiniteOpt._model_from_expr(quad1) === m
        @test InfiniteOpt._model_from_expr(quad2) isa Nothing
        @test InfiniteOpt._model_from_expr(quad3) === m
    end
    # test for Vector{GeneralVariableRef}
    @testset "Vector{GeneralVariableRef}" begin
        vrefs1 = GeneralVariableRef[]
        vrefs2 = [hd]
        # test expressions
        @test InfiniteOpt._model_from_expr(vrefs1) isa Nothing
        @test InfiniteOpt._model_from_expr(vrefs2) === m
    end
    # test Fallback
    @testset "Fallback" begin
        @test_throws ErrorException InfiniteOpt._model_from_expr(42)
    end
end

# Test _remove_variable
@testset "_remove_variable" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0), pt)
    @hold_variable(m, hold)
    # test for GenericAffExpr
    @testset "AffExpr" begin
        # make expressions
        aff1 = pt + 2par + hold
        aff2 = zero(GenericAffExpr{Float64, GeneralVariableRef})
        # test expressions
        @test isa(InfiniteOpt._remove_variable(aff1, hold), Nothing)
        @test !haskey(aff1.terms, hold)
        @test isa(InfiniteOpt._remove_variable(aff2, inf), Nothing)
    end
    # test for GenericQuadExpr
    @testset "QuadExpr" begin
        # make expressions
        quad1 = pt^2 + inf * pt + 2par + hold
        quad2 = zero(GenericQuadExpr{Float64, GeneralVariableRef})
        # test expressions
        quad = copy(quad1)
        @test isa(InfiniteOpt._remove_variable(quad, inf), Nothing)
        @test !haskey(quad.aff.terms, inf)
        @test !haskey(quad.terms, UnorderedPair{GeneralVariableRef}(inf, pt))
        quad = copy(quad1)
        @test isa(InfiniteOpt._remove_variable(quad, pt), Nothing)
        @test !haskey(quad.terms, UnorderedPair{GeneralVariableRef}(inf, pt))
        @test !haskey(quad.terms, UnorderedPair{GeneralVariableRef}(pt, pt))
        @test isa(InfiniteOpt._remove_variable(quad2, inf), Nothing)
    end
end

# Test _set_variable_coefficient!
@testset "_set_variable_coefficient!" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_parameter(m, 0 <= par2 <= 1)
    @infinite_variable(m, x(par))
    @infinite_variable(m, y(par, par2))
    @hold_variable(m, z)
    # test with GeneralVariableRef
    @testset "GeneralVariableRef" begin
        @test InfiniteOpt._set_variable_coefficient!(x, x, 2) == 2 * x
        @test InfiniteOpt._set_variable_coefficient!(z, x, 2) == z + 2x
    end
    # test with GenericAffExpr
    @testset "AffExpr" begin
        @test InfiniteOpt._set_variable_coefficient!(x + z, x, 2) == 2x + z
        @test InfiniteOpt._set_variable_coefficient!(x + z, y, 2) == x + z + 2y
    end
    # test with GenericQuadExpr
    @testset "QuadExpr" begin
        @test InfiniteOpt._set_variable_coefficient!(y ^2 + x, x, 2) == y^2 + 2x
        @test InfiniteOpt._set_variable_coefficient!(y^2 + x, y, 2) == y^2 + x + 2y
    end
    # test fallabck
    @testset "Fallback" begin
        @test_throws ErrorException InfiniteOpt._set_variable_coefficient!(2, x, 2)
    end
end

# Test parameter reference methods
@testset "Parameter References" begin
    m = InfiniteModel()
    @independent_parameter(m, t in [0, 1])
    @independent_parameter(m, y in [0, 1])
    @dependent_parameters(m, x[1:3] in [0, 1])
    @hold_variable(m, z)
    @expression(m, c1, 2z)
    @expression(m, c2, z + t + x[1])
    # test _make_param_tuple_element (IndependentParameterIndex)
    @testset "_make_param_tuple_element (Independent)" begin
        @test InfiniteOpt._make_param_tuple_element(m, index(t)) == t
        @test InfiniteOpt._make_param_tuple_element(m, index(y)) == y
    end
    # test _make_param_tuple_element (DependentParametersIndex)
    @testset "_make_param_tuple_element (Dependent)" begin
        obj_idx = index(first(x)).object_index
        @test InfiniteOpt._make_param_tuple_element(m, obj_idx) == x
    end
    # test parameter_refs
    @testset "parameter_refs" begin
        @test parameter_refs(c1) == ()
        @test parameter_refs(c2) == (t, x)
        @test parameter_refs(zero(AffExpr)) == ()
    end
end
