# Test parameter function methods 
@testset "Parameter Function Methods" begin 
    # setup the needed info 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [-1, 1])
    func = ParameterFunction(sin, IC.VectorTuple(t), [1], [1])
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
        @test isequal(dispatch_variable_ref(m, idx), fref)
        @test isequal(dispatch_variable_ref(gvref), fref)
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
        @test core_object(fref) === func
        @test core_object(gvref) === func
    end
    # parameter_group_int_indices
    @testset "parameter_group_int_indices" begin
        @test InfiniteOpt.parameter_group_int_indices(fref) == [1]
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
        @test InfiniteOpt._constraint_dependencies(fref) == InfOptConstraintIndex[]
        @test InfiniteOpt._constraint_dependencies(gvref) == InfOptConstraintIndex[]
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
        @test isequal(raw_parameter_refs(fref), IC.VectorTuple(t))
        @test isequal(raw_parameter_refs(gvref), IC.VectorTuple(t))
    end
    # parameter_refs
    @testset "parameter_refs" begin
        @test isequal(parameter_refs(fref), (t, ))
        @test isequal(parameter_refs(gvref), (t, ))
    end
    # parameter_list
    @testset "parameter_list" begin
        @test isequal(parameter_list(fref), [t])
        @test isequal(parameter_list(gvref), [t])
    end
    # raw_function
    @testset "raw_function" begin
        @test raw_function(fref) == sin
        @test raw_function(gvref) == sin
    end
    # call_function
    @testset "call_function" begin
        @test call_function(fref, 0) == 0
        @test call_function(gvref, 0) == 0
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
        push!(InfiniteOpt._constraint_dependencies(fref), InfOptConstraintIndex(1))
        @test used_by_constraint(fref)
        empty!(InfiniteOpt._constraint_dependencies(fref))
    end
    # test is_used
    @testset "is_used" begin
        @test !is_used(fref)
        push!(InfiniteOpt._constraint_dependencies(fref), InfOptConstraintIndex(1))
        @test is_used(fref)
        empty!(InfiniteOpt._constraint_dependencies(fref))
    end
    # test _update_param_var_mapping
    @testset "_update_param_var_mapping" begin 
        @test InfiniteOpt._update_param_var_mapping(fref, IC.VectorTuple(t)) isa Nothing 
        @test used_by_parameter_function(t)
        empty!(InfiniteOpt._parameter_function_dependencies(t))
    end
    # test _delete_data_object
    @testset "_delete_data_object" begin
        @test InfiniteOpt._delete_data_object(fref) isa Nothing
        @test length(InfiniteOpt._data_dictionary(fref)) == 0
        @test !is_valid(m, fref)
    end
    # test build_parameter_function
    @testset "build_parameter_function" begin 
        f3(ts, xs, a...; b...) = 42
        # test errors  
        @test_throws ErrorException build_parameter_function(error, sin, t, bob = 1)
        @test_throws ErrorException build_parameter_function(error, sin, x[1])
        @test_throws ErrorException build_parameter_function(error, sin, (t, x))
        @test_throws ErrorException build_parameter_function(error, sin, 2)
        # test normal  
        @test build_parameter_function(error, (a, b) -> 2, (t, x)) isa ParameterFunction 
        @test build_parameter_function(error, f3, (t, x)) isa ParameterFunction 
        @test build_parameter_function(error, (ts, xs) -> f3(ts, xs, 1), (t, x)) isa ParameterFunction 
        @test build_parameter_function(error, (ts, xs) -> f3(ts, xs, d = 1), (t, x)) isa ParameterFunction 
        @test build_parameter_function(error, sin, t) isa ParameterFunction 
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
        @test isequal(add_parameter_function(m, func, "test"), fref)
        @test name(fref) == "test"
        @test isequal(parameter_refs(fref), (t, x))
        @test used_by_parameter_function(t)
        # test default name 
        func = build_parameter_function(error, cos, t)
        fref = GeneralVariableRef(m, 3, ParameterFunctionIndex)
        @test isequal(add_parameter_function(m, func), fref)
        @test name(fref) == "cos"
    end
    # test parameter_function
    @testset "parameter_function" begin 
        f4(ts, xs, a...; b...) = 42
        # test normal
        @test parameter_function(sin, t) isa GeneralVariableRef
        @test parameter_function(sin, t, name = "name") isa GeneralVariableRef
        @test parameter_function((ts, xs) -> f4(ts, xs, 1, d = 1), (t, x)) isa GeneralVariableRef
        # test errors 
        @test_throws ErrorException parameter_function(sin, (t, x))
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

# Test the parameter function methods for macro definition 
@testset "Parameter Function Macro" begin
    # setup 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x[1:2] in [0, 1], independent = true)
    f5(t, x, a...; b...) = 42
    # test _expr_replace!
    @testset "_expr_replace!" begin
        @test InfiniteOpt._expr_replace!(:(t -> f(t, x)), :t, :y) == :(y -> f(y, x))
        @test InfiniteOpt._expr_replace!(:((t, x[1]) -> f(t, x[1])), :(x[1]), :y) == :((t, y) -> f(t, y))
        @test InfiniteOpt._expr_replace!(:t, :t, :y) == :y
        @test InfiniteOpt._expr_replace!(:t, :f, :y) == :t
    end
    # test _extract_parameters
    @testset "_extract_parameters" begin
        @test InfiniteOpt._extract_parameters(:y) == esc(:y)
        @test InfiniteOpt._extract_parameters(:(x[1])) == esc(:(x[1])) 
        ex = :((t, x[1]))
        @test InfiniteOpt._extract_parameters(ex).args[1] == ex
        @test InfiniteOpt._extract_parameters(ex).args[1] !== ex
    end
    # test _process_func_expr
    @testset "_process_func_expr" begin
        # test normal
        @test InfiniteOpt._process_func_expr(error, :(f(t, x))) == (esc(:f), esc(:(t, x)), false)
        anon = :((t, x) -> f(t, x, 1, d = 1))
        @test InfiniteOpt._process_func_expr(error, anon) == (esc(anon), esc(:(t,x)), true)
        anon = :((t, x[1]) -> sin(t + x[1]))
        # Need invokelatest due to https://github.com/infiniteopt/InfiniteOpt.jl/pull/369
        @test invokelatest(eval(InfiniteOpt._process_func_expr(error, anon)[1].args[1]), 0.5, 0.2) == sin(0.5 + 0.2)
        anon = :(t[i] -> t[i] + 3)
        @test invokelatest(eval(InfiniteOpt._process_func_expr(error, anon)[1].args[1]), 2) == 5
        # test errors
        @test_throws ErrorException InfiniteOpt._process_func_expr(error, :(f(t, d = 2)))
        @test_throws ErrorException InfiniteOpt._process_func_expr(error, :(f(x, t; d = 2)))
        @test_throws ErrorException InfiniteOpt._process_func_expr(error, :(f >= 2))
        @test_throws ErrorException InfiniteOpt._process_func_expr(error, :(f[t, x]))
    end
    # test @parameter_function
    @testset "@parameter_function" begin
        # test errors
        @test_macro_throws ErrorException @parameter_function()
        @test_macro_throws ErrorException @parameter_function(m)
        @test_macro_throws ErrorException @parameter_function(m, func = f5)
        @test_macro_throws ErrorException @parameter_function(m, y == sin(t), Int)
        @test_macro_throws ErrorException @parameter_function(m, 2 == sin(t))
        @test_macro_throws ErrorException @parameter_function(m, [1:2])
        @test_macro_throws ErrorException @parameter_function(m, "a$(1)" == f5(t, x))
        @test_macro_throws ErrorException @parameter_function(m, a[m = 1:2] == f5(t, x))
        @test_macro_throws ErrorException @parameter_function(Model(), sin(t))
        @test_macro_throws ErrorException @parameter_function(m, [a...] == sin(t))
        # test anonymous singular 
        idx = 1
        ref = GeneralVariableRef(m, idx, ParameterFunctionIndex)
        @test isequal(@parameter_function(m, f5(t, x); base_name = "a"), ref)
        @test name(ref) == "a"
        @test raw_function(ref) == f5 
        @test isequal(parameter_refs(ref), (t, x))
        idx += 1
        ref = GeneralVariableRef(m, idx, ParameterFunctionIndex)
        @test isequal(@parameter_function(m, f5(t, x)), ref)
        @test name(ref) == "f5"
        @test raw_function(ref) == f5
        @test isequal(parameter_refs(ref), (t, x))
        idx += 1
        # test anonymous single argument multi-dim 
        refs = [GeneralVariableRef(m, idx + i, ParameterFunctionIndex) for i in 0:1]
        @test isequal(@parameter_function(m, [1:2] == f5(t, x)), refs)
        @test isequal(parameter_refs(refs[1]), (t, x))
        @test raw_function(refs[2]) == f5 
        @test name.(refs) == ["f5", "f5"]
        idx += 2
        refs = [GeneralVariableRef(m, idx + i, ParameterFunctionIndex) for i in 0:1]
        @test @parameter_function(m, [i = 1:2; i >= 1] == (sin, cos)[i](t)) isa JuMPC.SparseAxisArray
        @test isequal(parameter_refs(refs[1]), (t,))
        @test raw_function(refs[2]) == cos
        @test raw_function(refs[1]) == sin
        @test name(refs[1]) == "sin"
        @test name(refs[2]) == "cos"
        idx += 2
        # test explicit single 
        ref = GeneralVariableRef(m, idx, ParameterFunctionIndex)
        @test isequal(@parameter_function(m, a == f5(t, x), base_name = "bob"), ref)
        @test name(ref) == "bob"
        @test raw_function(ref) == f5 
        @test isequal(parameter_refs(ref), (t, x))
        idx += 1
        ref = GeneralVariableRef(m, idx, ParameterFunctionIndex)
        @test isequal(@parameter_function(m, c == (t, x) -> f5(t, x, 2, d = 1)), ref)
        @test name(ref) == "c"
        @test raw_function(ref) != f5 
        @test raw_function(ref)(1, [1, 1]) == 42
        @test isequal(parameter_refs(ref), (t, x))
        idx += 1
        # test explicit multi-dim 
        refs = [GeneralVariableRef(m, idx + i, ParameterFunctionIndex) for i in 0:1]
        @test isequal(@parameter_function(m, d[1:2] == f5(t, x)), refs)
        @test isequal(parameter_refs(refs[1]), (t, x))
        @test raw_function(refs[2]) == f5 
        @test name.(refs) == ["d[1]", "d[2]"]
        idx += 2
        refs = [GeneralVariableRef(m, idx + i, ParameterFunctionIndex) for i in 0:1]
        @test isequal(@parameter_function(m, e[1:2] == (t, x) -> f5(t, x, 2; s = 1)), refs)
        @test isequal(parameter_refs(refs[1]), (t, x))
        @test raw_function(refs[2]) != f5 
        @test name.(refs) == ["e[1]", "e[2]"]
        idx += 2
        # test infinite parameter with reference 
        refs = [GeneralVariableRef(m, idx + i, ParameterFunctionIndex) for i in 0:1]
        @test isequal(@parameter_function(m, [i = 1:2] == (t, x[i]) -> sin(t + x[i])), refs)
        @test isequal(parameter_refs.(refs), [(t, x[1]), (t, x[2])])
        @test call_function.(refs, 0.5, 0.2) == [sin(0.5 + 0.2), sin(0.5 + 0.2)] 
        @test name.(refs) == ["", ""]
        idx += 2
    end
end

# Test _interrogate_variables
@testset "_interrogate_variables" begin 
    # setup model
    m = InfiniteModel()
    @variable(m, z)
    @variable(m, y)
    aff = 2z + 42
    quad = z^2 + 2z
    nlp = sin(z) + aff
    # test constant 
    @testset "Constant" begin 
        @test InfiniteOpt._interrogate_variables(i -> error, 42) isa Nothing
    end
    # test variable
    @testset "Variable" begin 
        a = [0]
        @test InfiniteOpt._interrogate_variables(i -> a[1] += 1, z) isa Nothing
        @test a[1] == 1
    end
    # test AffExpr
    @testset "AffExpr" begin
        a = []
        @test InfiniteOpt._interrogate_variables(i -> push!(a, i), aff) isa Nothing
        @test isequal(a, [z])
    end
    # test QuadExpr
    @testset "QuadExpr" begin
        a = []
        @test InfiniteOpt._interrogate_variables(i -> push!(a, i), quad) isa Nothing
        @test isequal(a, [z, z, z])
    end
    # test GenericNonlinearExpr
    @testset "GenericNonlinearExpr" begin
        a = []
        @test InfiniteOpt._interrogate_variables(i -> push!(a, i), nlp) isa Nothing
        @test isequal(a, [z, z])
    end
    # test array of expressions 
    @testset "AbstractArray" begin
        a = []
        @test InfiniteOpt._interrogate_variables(i -> push!(a, i), [nlp, aff]) isa Nothing
        @test isequal(a, [z, z, z])
    end
end

# Test all_expression_variables
@testset "all_expression_variables" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    data = TestData(par, 0, 5)
    meas = Measure(finite, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    meas = InfiniteOpt.GeneralVariableRef(m, mindex)
    dinf = @deriv(inf, par)
    # test for variable reference
    @testset "Variable" begin
        @test isequal(all_expression_variables(par), [par])
        @test isequal(all_expression_variables(inf), [inf])
        @test isequal(all_expression_variables(meas), [meas])
        @test isequal(all_expression_variables(dinf), [dinf])
    end
    # test for GenericAffExpr
    @testset "AffExpr" begin
        # make expressions
        aff1 = meas + 2par + finite - dinf
        aff2 = zero(GenericAffExpr{Float64, GeneralVariableRef})
        # test expressions
        @test isempty(setdiff(all_expression_variables(aff1),
                              [meas, par, finite, dinf]))
        @test all_expression_variables(aff2) == GeneralVariableRef[]
    end
    # test for GenericQuadExpr
    @testset "QuadExpr" begin
        # make expressions
        quad1 = pt^2 + inf * pt - meas + 2par + finite - dinf
        quad2 = pt^2 + inf * pt
        quad3 = zero(GenericQuadExpr{Float64, GeneralVariableRef})
        # test expressions
        @test isempty(setdiff(all_expression_variables(quad1),
                      [meas, par, finite, dinf, pt, inf]))
        @test isempty(setdiff(all_expression_variables(quad2),
                              [pt, inf]))
        @test all_expression_variables(quad3) == GeneralVariableRef[]
    end
    # test for Array of expressions
    @testset "AbstractArray" begin
        # make expressions
        ex1 = [inf, pt]
        ex2 = [inf + pt, meas + pt]
        # test expressions
        @test isempty(setdiff(all_expression_variables(ex1),
                      [pt, inf]))
        @test isempty(setdiff(all_expression_variables(ex2),
                              [pt, inf, meas]))
    end
    # test for Array of expressions
    @testset "GenericNonlinearExpr" begin
        # make expressions
        nlp = sin(pt) + inf / pt
        # test expressions
        @test isempty(setdiff(all_expression_variables(nlp),
                      [pt, inf]))
    end
    # test backup
    @testset "Fallback" begin
        @variable(Model(), x)
        @test_throws ErrorException all_expression_variables(x)
    end
end

# Test parameter_group_int_indices
@testset "parameter_group_int_indices" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, inf2, Infinite(par, pars))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    var = build_variable(error, inf2, Dict{Int, Float64}(1 => 0.5), check = false)
    red = add_variable(m, var)
    # test for finite variable reference
    @testset "FiniteVariable" begin
        @test InfiniteOpt.parameter_group_int_indices(pt) == []
        @test InfiniteOpt.parameter_group_int_indices(finite) == []
    end
    # test for infinite variable reference
    @testset "InfiniteVariable" begin
        @test InfiniteOpt.parameter_group_int_indices(inf) == [1]
    end
    # test for parameter reference
    @testset "Parameter" begin
        @test InfiniteOpt.parameter_group_int_indices(par) == [1]
        @test InfiniteOpt.parameter_group_int_indices(pars[1]) == [2]
    end
    # test for semi-infinite variable reference
    @testset "SemiInfiniteInfinite" begin
        @test InfiniteOpt.parameter_group_int_indices(red) == [2]
    end
    # test for GenericAffExpr
    @testset "AffExpr" begin
        # make expressions
        aff1 = inf + inf2 + pt - 3
        aff2 = pt + finite - 2
        # test expressions
        @test sort!(InfiniteOpt.parameter_group_int_indices(aff1)) == [1, 2]
        @test InfiniteOpt.parameter_group_int_indices(aff2) == []
    end
    # test for GenericQuadExpr
    @testset "QuadExpr" begin
        # make expressions
        quad1 = inf * inf2 + inf + inf2 + pt - 3 - par
        quad2 = pt * pt + pt + finite - 2
        # test expressions
        @test sort!(InfiniteOpt.parameter_group_int_indices(quad1)) == [1, 2]
        @test InfiniteOpt.parameter_group_int_indices(quad2) == []
    end
    # test for GenericNonlinearExpr
    @testset "GenericNonlinearExpr" begin
        # make expressions
        nlp = sin(inf) / pt
        # test expressions
        @test InfiniteOpt.parameter_group_int_indices(nlp) == [1]
    end
end

# Test _parameter_numbers
@testset "_parameter_numbers" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, inf2, Infinite(par, pars))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    var = build_variable(error, inf2, Dict{Int, Float64}(1 => 0.5), check = false)
    red = add_variable(m, var)
    # test for finite variable reference
    @testset "FiniteVariable" begin
        @test InfiniteOpt._parameter_numbers(pt) == []
        @test InfiniteOpt._parameter_numbers(finite) == []
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
        aff2 = pt + finite - 2
        # test expressions
        @test sort!(InfiniteOpt._parameter_numbers(aff1)) == [1, 2, 3]
        @test InfiniteOpt._parameter_numbers(aff2) == []
    end
    # test for GenericQuadExpr
    @testset "QuadExpr" begin
        # make expressions
        quad1 = inf * inf2 + inf + inf2 + pt - 3 - par
        quad2 = pt * pt + pt + finite - 2
        # test expressions
        @test sort!(InfiniteOpt._parameter_numbers(quad1)) == [1, 2, 3]
        @test InfiniteOpt._parameter_numbers(quad2) == []
    end
    # test for GenericNonlinearExpr
    @testset "GenericNonlinearExpr" begin
        # make expressions
        nlp = sin(inf2)
        # test expressions
        @test sort!(InfiniteOpt._parameter_numbers(nlp)) == [1, 2, 3]
    end
end

# Test _remove_variable
@testset "_remove_variable" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0))
    @variable(m, finite)
    # test for GenericAffExpr
    @testset "AffExpr" begin
        # make expressions
        aff1 = pt + 2par + finite
        aff2 = zero(GenericAffExpr{Float64, GeneralVariableRef})
        # test expressions
        @test isa(InfiniteOpt._remove_variable(aff1, finite), Nothing)
        @test !haskey(aff1.terms, finite)
        @test isa(InfiniteOpt._remove_variable(aff2, inf), Nothing)
    end
    # test for GenericQuadExpr
    @testset "QuadExpr" begin
        # make expressions
        quad1 = pt^2 + inf * pt + 2par + finite
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
    # test for GenericNonlinearExpr 
    @testset "GenericNonlinearExpr" begin 
        # make expressions 
        nlp1 = sin(3pt)
        nlp2 = pt^2.3 + ^(inf, pt)
        nlp3 = cos(pt^2 + pt) / (2pt + 2inf)
        # test expressions 
        @test InfiniteOpt._remove_variable(nlp1, pt) isa Nothing 
        @test isequal(nlp1, sin(zero(GenericAffExpr{Float64, GeneralVariableRef})))
        @test InfiniteOpt._remove_variable(nlp2, pt) isa Nothing 
        @test nlp2.args[1].args == [0.0, 2.3]
        @test nlp2.args[2].args == Any[inf, 0.0]
        @test InfiniteOpt._remove_variable(nlp3, inf) isa Nothing 
        @test nlp3.args[2] == 2pt
    end
    # test for AbstractArray
    @testset "AbstractArray" begin
        # make expressions
        ex = [pt + 2par + finite, 2inf + pt]
        # test expressions
        @test isa(InfiniteOpt._remove_variable(ex, finite), Nothing)
        @test !haskey(ex[1].terms, finite)
    end
end

# Test map_expression
@testset "map_expression" begin 
    # setup model
    m = InfiniteModel()
    @variable(m, z)
    @variable(m, y)
    aff = 2z + 42
    quad = z^2 + 2z
    nlp = (sin(z) + aff) ^ 3.4
    @variable(Model(), x)
    # test variable 
    @testset "Variable" begin
        @test isequal(map_expression(v -> x, z), x)
    end
    # test AffExpr
    @testset "AffExpr" begin
        @test isequal(map_expression(v -> x, aff), 2x + 42)
    end
    # test QuadExpr 
    @testset "QuadExpr" begin
        @test isequal(map_expression(v -> x, quad), x^2 + 2x)
    end
    # test GenericNonlinearExpr
    @testset "GenericNonlinearExpr" begin
        @test isequal(map_expression(v -> y, nlp), (sin(y) + (2y + 42))^3.4)
        @test isequal(map_expression(v -> v^3, sin(y)), sin(y^3))
    end
end

# Test map_expression_to_ast
@testset "map_expression_to_ast" begin 
    # setup model
    m = InfiniteModel()
    @variable(m, z)
    @variable(m, y)
    aff = 2z + y + 42
    quad = z^2 + 3 * z * y + 2z
    nlp = (sin(z) + aff) ^ 3.4
    aff0 = zero(GenericAffExpr{Float64, GeneralVariableRef})
    quad2 = 2 * z^2
    quad0 = zero(GenericQuadExpr{Float64, GeneralVariableRef})
    jm = Model()
    @variable(jm, x)
    @variable(jm, w)
    vmap(v) = v == z ? x : w
    omap(op) = :test
    # test constant
    @testset "Constant" begin
        @test map_expression_to_ast(vmap, 42) == :(42)
    end
    # test variable 
    @testset "Variable" begin
        @test map_expression_to_ast(vmap, z) == :($x)
    end
    # test AffExpr
    @testset "AffExpr" begin
        @test map_expression_to_ast(vmap, aff) == :(2 * $x + $w + 42)
        @test map_expression_to_ast(vmap, aff0) == :(+(0))
    end
    # test QuadExpr 
    @testset "QuadExpr" begin
        @test map_expression_to_ast(vmap, quad) == :($x * $x + 3 * $x * $w + 2 * $x)
        @test map_expression_to_ast(vmap, quad2) == :(+(2 * $x * $x))
        @test map_expression_to_ast(vmap, quad0) == :(+(0))
    end
    # test GenericNonlinearExpr
    @testset "GenericNonlinearExpr" begin
        @test map_expression_to_ast(vmap, omap, nlp) == :(test(test(test($x), (2 * $x + $w + 42)), 3.4))
    end
    # test deprecation 
    @testset "map_nlp_to_ast" begin
        @test (@test_deprecated map_nlp_to_ast(vmap, nlp)) == :((sin($x) + (2.0 * $x + $w + 42.0)) ^ 3.4)
    end
end

# Test restrictions
@testset "restrict" begin 
    # setup model
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, x in [-1, 1])
    @variable(m, y, Infinite(t, x))
    @variable(m, q, Infinite(x, t))
    @variable(m, w, Infinite(t))
    @variable(m, z)
    # test AffExpr
    @testset "AffExpr" begin
        @test (2z + y + 42)(0, -1) == 2z + y(0, -1) + 42
        @test (2z + y + 42)(0, x) == 2z + y(0, x) + 42
        @test_throws ErrorException (y + q)(0, 0)
        @test_throws ErrorException (y + t)(0, x)
    end
    # test QuadExpr 
    @testset "QuadExpr" begin
        @test (z^2 + 3y - 2)(0, -1) == z^2 + 3y(0, -1) - 2
        @test (w^2 - 2)(0) == w(0)^2 - 2
    end
    # test GenericNonlinearExpr
    @testset "GenericNonlinearExpr" begin
        @test isequal_canonical((sin(y) * z)(t, -1), sin(y(t, -1)) * z)
    end
end

# Test _set_variable_coefficient!
@testset "_set_variable_coefficient!" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, par2 in [0, 1])
    @variable(m, x, Infinite(par))
    @variable(m, y, Infinite(par, par2))
    @variable(m, z)
    # test with GeneralVariableRef
    @testset "GeneralVariableRef" begin
        @test isequal_canonical(InfiniteOpt._set_variable_coefficient!(x, x, 2), 2 * x)
        @test isequal_canonical(InfiniteOpt._set_variable_coefficient!(z, x, 2), z + 2x)
    end
    # test with GenericAffExpr
    @testset "AffExpr" begin
        @test isequal_canonical(InfiniteOpt._set_variable_coefficient!(x + z, x, 2), 2x + z)
        @test isequal_canonical(InfiniteOpt._set_variable_coefficient!(x + z, y, 2), x + z + 2y)
    end
    # test with GenericQuadExpr
    @testset "QuadExpr" begin
        @test isequal_canonical(InfiniteOpt._set_variable_coefficient!(y ^2 + x, x, 2), y^2 + 2x)
        @test isequal_canonical(InfiniteOpt._set_variable_coefficient!(y^2 + x, y, 2), y^2 + x + 2y)
    end
    # test fallabck
    @testset "Fallback" begin
        @test_throws ErrorException InfiniteOpt._set_variable_coefficient!(2, x, 2)
    end
end

# Test parameter reference methods
@testset "Parameter References" begin
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, y in [0, 1])
    @infinite_parameter(m, x[1:3] in [0, 1])
    @variable(m, z)
    @expression(m, c1, 2z)
    @expression(m, c2, z + t + x[1])
    # test _make_param_tuple_element (IndependentParameterIndex)
    @testset "_make_param_tuple_element (Independent)" begin
        @test isequal(InfiniteOpt._make_param_tuple_element(m, index(t)), t)
        @test isequal(InfiniteOpt._make_param_tuple_element(m, index(y)), y)
    end
    # test _make_param_tuple_element (DependentParametersIndex)
    @testset "_make_param_tuple_element (Dependent)" begin
        obj_idx = index(first(x)).object_index
        @test isequal(InfiniteOpt._make_param_tuple_element(m, obj_idx), x)
    end
    # test parameter_refs
    @testset "parameter_refs" begin
        @test parameter_refs(c1) == ()
        @test isequal(parameter_refs(c2), (t, x))
        @test parameter_refs(zero(AffExpr)) == ()
    end
end
