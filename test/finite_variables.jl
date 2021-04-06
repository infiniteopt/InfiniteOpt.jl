# Test name methods
@testset "Basics" begin
    # initialize model and variable
    m = InfiniteModel()
    num = Float64(0)
    info = VariableInfo(false, num, false, num, false, num, false, num, false, false)
    new_info = VariableInfo(true, 0., true, 0., true, 0., true, 0., true, false)
    bounds = ParameterBounds()
    var = FiniteVariable(info, bounds)
    object = VariableData(var, "test")
    idx = FiniteVariableIndex(1)
    vref = FiniteVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, FiniteVariableIndex)
    bad_vref = FiniteVariableRef(m, FiniteVariableIndex(-1))
    # JuMP.owner_model
    @testset "JuMP.owner_model" begin
        @test owner_model(vref) === m
        @test owner_model(gvref) === m
    end
    # JuMP.index
    @testset "JuMP.index" begin
        @test index(vref) == idx
        @test index(gvref) == idx
    end
    # dispatch_variable_ref
    @testset "dispatch_variable_ref" begin
        @test dispatch_variable_ref(m, idx) == vref
        @test dispatch_variable_ref(gvref) == vref
    end
    # _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, object) == idx
    end
    # _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(m, FiniteVariable) === m.finite_vars
        @test InfiniteOpt._data_dictionary(vref) === m.finite_vars
        @test InfiniteOpt._data_dictionary(gvref) === m.finite_vars
    end
    # JuMP.is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, vref)
        @test is_valid(m, gvref)
    end
    # _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(vref) === object
        @test InfiniteOpt._data_object(gvref) === object
        @test_throws ErrorException InfiniteOpt._data_object(bad_vref)
    end
    # _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(vref) === var
        @test InfiniteOpt._core_variable_object(gvref) === var
    end
    @testset "_variable_info" begin
        @test InfiniteOpt._variable_info(vref) == info
    end
    # _update_variable_info
    @testset "_update_variable_info" begin
        @test isa(InfiniteOpt._update_variable_info(vref, new_info), Nothing)
        @test InfiniteOpt._variable_info(vref) == new_info
    end
    # _set_core_variable_object
    @testset "_set_core_variable_object" begin
        @test InfiniteOpt._set_core_variable_object(vref, var) isa Nothing
    end
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "test"
        @test name(gvref) == "test"
        @test name(bad_vref) == ""
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
        @test isa(set_name(gvref, "new2"), Nothing)
        @test name(vref) == "new2"
        @test_throws ErrorException set_name(bad_vref, "")
    end
    # test parameter_refs
    @testset "parameter_refs" begin
        @test parameter_refs(vref) == ()
        @test parameter_refs(gvref) == ()
    end
    # _make_variable_ref
    @testset "_make_variable_ref" begin
        @test InfiniteOpt._make_variable_ref(m, idx) == gvref
    end
    # _var_name_dict
    @testset "_var_name_dict" begin
        @test InfiniteOpt._var_name_dict(m) isa Nothing
    end
    # _update_var_name_dict
    @testset "_update_var_name_dict" begin
        m.name_to_var = Dict{String, ObjectIndex}()
        dict = InfiniteOpt._data_dictionary(vref)
        @test InfiniteOpt._update_var_name_dict(m, dict) isa Nothing
        @test InfiniteOpt._var_name_dict(m) == Dict(name(vref) => idx)
        m.name_to_var = nothing
    end
    # parameter_by_name
    @testset "JuMP.variable_by_name" begin
        # test normal
        @test variable_by_name(m, "new2") == gvref
        @test variable_by_name(m, "bob") isa Nothing
        # prepare variable with same name
        idx2 = FiniteVariableIndex(2)
        @test InfiniteOpt._add_data_object(m, object) == idx2
        vref2 = FiniteVariableRef(m, idx2)
        @test set_name(vref2, "new2") isa Nothing
        # test multiple name error
        @test_throws ErrorException variable_by_name(m, "new2")
    end
    # _delete_data_object
    @testset "_delete_data_object" begin
        @test InfiniteOpt._delete_data_object(vref) isa Nothing
        @test length(InfiniteOpt._data_dictionary(vref)) == 1
        @test !is_valid(m, vref)
    end
end

# Test variable definition methods
@testset "Definition" begin
    # initialize model and info
    m = InfiniteModel()
    num = Float64(0)
    info = VariableInfo(false, num, false, num, false, num, false, num, false, false)
    info2 = VariableInfo(true, num, true, num, true, num, true, num, true, false)
    info3 = VariableInfo(true, num, true, num, true, num, true, num, false, true)
    @infinite_parameter(m, 0 <= par <= 10)
    @infinite_parameter(m, 0 <= pars[1:2] <= 10)
    @finite_parameter(m, fin, 42)
    # test _check_bounds
    @testset "_check_bounds" begin
        # test normal
        bounds = ParameterBounds((par => IntervalSet(0, 1),))
        @test InfiniteOpt._check_bounds(bounds) isa Nothing
        bounds = ParameterBounds((pars => IntervalSet(0, 1),))
        @test InfiniteOpt._check_bounds(bounds) isa Nothing
        bounds = ParameterBounds((pars => IntervalSet(0, 0),))
        @test InfiniteOpt._check_bounds(bounds) isa Nothing
        bounds = ParameterBounds((pars[1] => IntervalSet(0, 1),))
        @test InfiniteOpt._check_bounds(bounds) isa Nothing
        # test errors
        bounds = ParameterBounds((fin => IntervalSet(0, 1),))
        @test_throws ErrorException InfiniteOpt._check_bounds(bounds)
        bounds = ParameterBounds((par => IntervalSet(-1, 1),))
        @test_throws ErrorException InfiniteOpt._check_bounds(bounds)
        bounds = ParameterBounds((par => IntervalSet(0, 11),))
        @test_throws ErrorException InfiniteOpt._check_bounds(bounds)
        bounds = ParameterBounds((pars[1] => IntervalSet(0, 0),))
        @test_throws ErrorException InfiniteOpt._check_bounds(bounds)
        bounds = ParameterBounds((pars[1] => IntervalSet(0, 0),
                                 pars[2] => IntervalSet(0, 3)))
        @test_throws ErrorException InfiniteOpt._check_bounds(bounds)
    end
    # _make_variable
    @testset "_make_variable" begin
        # test normal
        @test InfiniteOpt._make_variable(error, info, Finite).info == info
        bounds = ParameterBounds((pars => IntervalSet(0, 0),))
        @test InfiniteOpt._make_variable(error, info, Finite,
                                         parameter_bounds = bounds).info == info
        @test InfiniteOpt._make_variable(error, info, Finite,
                            parameter_bounds = bounds).parameter_bounds == bounds
        # test errors
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                                Finite, parameter_values = 3)
        bounds = ParameterBounds((pars[1] => IntervalSet(0, 0),))
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                           Finite, parameter_bounds = bounds)
        bounds = ParameterBounds((par => IntervalSet(-1, -1),))
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                           Finite, parameter_bounds = bounds)
    end
    # build_variable
    @testset "JuMP.build_variable" begin
        # test normal
        @test build_variable(error, info, Finite).info == info
        @test build_variable(error, info, Finite).parameter_bounds == ParameterBounds()
        bounds = ParameterBounds((pars => IntervalSet(0, 5),))
        @test build_variable(error, info, Finite,
                             parameter_bounds = bounds).info == info
        @test build_variable(error, info, Finite,
                             parameter_bounds = bounds).parameter_bounds == bounds
        # test errors
        @test_throws ErrorException build_variable(error, info, Point,
                                                   parameter_bounds = bounds)
        bounds = ParameterBounds((par => IntervalSet(0, 0),
                                  pars => IntervalSet(-1, 1)))
        @test_throws ErrorException build_variable(error, info, Finite,
                                                   parameter_bounds = bounds)
    end
    # _validate_bounds
    @testset "_validate_bounds" begin
        # test normal
        bounds = ParameterBounds((par => IntervalSet(0, 1),))
        @test InfiniteOpt._validate_bounds(m, bounds) isa Nothing
        @test num_supports(par) == 0
        # test error
        @infinite_parameter(InfiniteModel(), par2 in Normal())
        bounds = ParameterBounds((par2 => IntervalSet(0, 1),))
        @test_throws VariableNotOwned{GeneralVariableRef} InfiniteOpt._validate_bounds(m, bounds)
        # test support addition
        bounds = ParameterBounds((par => IntervalSet(0, 0),))
        @test InfiniteOpt._validate_bounds(m, bounds) isa Nothing
        @test supports(par, label = UserDefined) == [0]
        bounds = ParameterBounds((pars => IntervalSet(0, 0),))
        @test InfiniteOpt._validate_bounds(m, bounds) isa Nothing
        @test supports(pars) == zeros(2, 1)
    end
    # _check_and_make_variable_ref
    @testset "_check_and_make_variable_ref" begin
        # test normal
        v = build_variable(error, info, Finite)
        idx = FiniteVariableIndex(1)
        vref = FiniteVariableRef(m, idx)
        @test InfiniteOpt._check_and_make_variable_ref(m, v, "") == vref
        # test with bounds
        bounds = ParameterBounds((par => IntervalSet(0, 2),))
        v = build_variable(error, info, Finite, parameter_bounds = bounds)
        idx = FiniteVariableIndex(2)
        vref = FiniteVariableRef(m, idx)
        @test InfiniteOpt._check_and_make_variable_ref(m, v, "") == vref
        @test m.has_finite_var_bounds
        m.has_finite_var_bounds = false
        # test bad bounds
        @infinite_parameter(InfiniteModel(), par2 in [0, 2])
        v = build_variable(error, info, Finite,
            parameter_bounds = ParameterBounds((par2 => IntervalSet(0, 1),)))
        @test_throws VariableNotOwned{GeneralVariableRef} InfiniteOpt._check_and_make_variable_ref(m, v, "")
    end
    # add_variable
    @testset "JuMP.add_variable" begin
        idx = FiniteVariableIndex(3)
        vref = FiniteVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        v = build_variable(error, info, Finite)
        @test add_variable(m, v, "name") == gvref
        @test haskey(InfiniteOpt._data_dictionary(vref), idx)
        @test name(vref) == "name"
        # prepare infinite variable with all the possible info additions
        v = build_variable(error, info2, Finite)
        # test info addition functions
        idx = FiniteVariableIndex(4)
        vref = FiniteVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test add_variable(m, v, "name") == gvref
        @test !optimizer_model_ready(m)
        @test !m.has_finite_var_bounds
        # lower bound
        cindex = ConstraintIndex(1)
        cref = InfOptConstraintRef(m, cindex)
        @test has_lower_bound(vref)
        @test InfiniteOpt._lower_bound_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.GreaterThan{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # upper bound
        cindex = ConstraintIndex(2)
        cref = InfOptConstraintRef(m, cindex)
        @test has_upper_bound(vref)
        @test InfiniteOpt._upper_bound_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.LessThan{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # fix
        cindex = ConstraintIndex(3)
        cref = InfOptConstraintRef(m, cindex)
        @test is_fixed(vref)
        @test InfiniteOpt._fix_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.EqualTo{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # binary
        cindex = ConstraintIndex(4)
        cref = InfOptConstraintRef(m, cindex)
        @test is_binary(vref)
        @test InfiniteOpt._binary_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.ZeroOne}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(vref) == [ConstraintIndex(i)
                                                             for i = 1:4]
        # prepare infinite variable with integer info addition
        bounds = ParameterBounds((par => IntervalSet(0, 2),
                                  pars => IntervalSet(1, 1)))
        v = build_variable(error, info3, Finite, parameter_bounds = bounds)
        # test integer addition functions
        idx = FiniteVariableIndex(5)
        vref = FiniteVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test add_variable(m, v, "name") == gvref
        @test !optimizer_model_ready(m)
        @test m.has_finite_var_bounds
        @test supports(par) == [0]
        @test sortcols(supports(pars)) == [0 1; 0 1]
        cindex = ConstraintIndex(8)
        cref = InfOptConstraintRef(m, cindex)
        @test is_integer(vref)
        @test InfiniteOpt._integer_index(vref) == cindex
        @test constraint_object(cref) isa BoundedScalarConstraint{GeneralVariableRef,
                                                            MOI.Integer}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(vref) == [ConstraintIndex(i)
                                                             for i = 5:8]
    end
end

# Test the finite variable macro if no bounds are provided
@testset "Macro" begin
    # initialize model
    m = InfiniteModel()
    @testset "@finite_variable (Not Bounded)" begin
        # test regular
        idx = FiniteVariableIndex(1)
        vref = FiniteVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @finite_variable(m, x >= 1, Bin) == gvref
        @test name(vref) == "x"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        # test anan
        idx = FiniteVariableIndex(2)
        vref = FiniteVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @finite_variable(m, binary = true, lower_bound = 1,
                               base_name = "x") == gvref
        @test name(vref) == "x"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        # test array
        idxs = [FiniteVariableIndex(3), FiniteVariableIndex(4)]
        vrefs = [FiniteVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt._make_variable_ref(m, idx) for idx in idxs]
        @test @finite_variable(m, y[1:2] == 2, Int) == gvrefs
        @test name(vrefs[1]) == "y[1]"
        @test fix_value(vrefs[2]) == 2
        @test is_integer(vrefs[1])
        # test errors
        @test_throws AssertionError @finite_variable(Model(), z >= 1, Bin)
        @test_macro_throws ErrorException @finite_variable(m, x >= 1, Bin)
    end
    # test _make_interval_set
    @testset "_make_interval_set" begin
        expected = :(IntervalSet(0, 1))
        @test InfiniteOpt._make_interval_set(error, :([0, 1])) == expected
        @test_throws ErrorException InfiniteOpt._make_interval_set(error, :([0]))
    end
    # test _make_bound_pair (in)
    @testset "_make_bound_pair (in)" begin
        expr = :(t in [0, 1])
        set = :(IntervalSet(0, 1))
        tp = :(t)
        @test InfiniteOpt._make_bound_pair(error, expr,
                                           Val(:in)) == :($(tp) => $(set))
    end
    # test _make_bound_pair (==)
    @testset "_make_bound_pair (==)" begin
        expr = :(t == 0)
        set = :(IntervalSet(0, 0))
        tp = :(t)
        @test InfiniteOpt._make_bound_pair(error, expr,
                                           Val(:(==))) == :($(tp) => $(set))
    end
    # test _make_bound_pair (comparison with <=)
    @testset "_make_bound_pair (<=)" begin
        expr = :(0 <= t <= 1)
        set = :(IntervalSet(0, 1))
        tp = :(t)
        @test InfiniteOpt._make_bound_pair(error, expr, Val(:(<=)),
                                           Val(:(<=))) == :($(tp) => $(set))
    end
    # test _make_bound_pair (comparison with >=)
    @testset "_make_bound_pair (>=)" begin
        expr = :(1 >= t >= 0)
        set = :(IntervalSet(0, 1))
        tp = :(t)
        @test InfiniteOpt._make_bound_pair(error, expr, Val(:(>=)),
                                           Val(:(>=))) == :($(tp) => $(set))
    end
    # test _make_bound_pair (comparison fallback)
    @testset "_make_bound_pair (Fallback)" begin
        expr = :(t <= 0)
        @test_throws ErrorException InfiniteOpt._make_bound_pair(error, expr,
                                                                 Val(:(<=)))
        @test_throws ErrorException InfiniteOpt._make_bound_pair(error, expr,
                                                         Val(:(<=)), Val(:(>=)))
    end
    # test _make_bound_pair (wrapper)
    @testset "_make_bound_pair (Wrapper)" begin
        # test :call expression
        expr = :(t in [0, 1])
        set = :(IntervalSet(0, 1))
        tp = :(t)
        @test InfiniteOpt._make_bound_pair(error, expr,
                                           Val(:in)) == :($(tp) => $(set))
        # test :compariosn expression
        expr = :(0 <= t <= 1)
        @test InfiniteOpt._make_bound_pair(error, expr, Val(:(<=)),
                                           Val(:(<=))) == :($(tp) => $(set))
        # test error
        expr = :([2, 3])
        @test_throws ErrorException InfiniteOpt._make_bound_pair(error, expr)
    end
    # test _parse_parameter_bounds (Vector)
    @testset "_parse_parameter_bounds (Vector)" begin
        args = [:(t in [0, 1]), :(0 <= t <= 1)]
        set = :(IntervalSet(0, 1))
        tp = :(t)
        dict_arg = :($(tp) => $(set))
        dict = :(($(dict_arg), $(dict_arg)))
        bounds = :(ParameterBounds($dict))
        @test InfiniteOpt._parse_parameter_bounds(error, args) == bounds
    end
    # test _parse_parameter_bounds (Expr)
    @testset "_parse_parameter_bounds (Expr)" begin
        expr = :(t in [0, 1])
        set = :(IntervalSet(0, 1))
        tp = :(t)
        dict_arg = :($(tp) => $(set))
        dict = :(($(dict_arg), ))
        bounds = :(ParameterBounds($dict))
        @test InfiniteOpt._parse_parameter_bounds(error, expr) == bounds
    end
    # test _extract_bounds (:call)
    @testset "_extract_bounds (:call)" begin
        # test single anonymous bound
        args = [:in, :t, :([0, 1])]
        set = :(IntervalSet(0, 1))
        tp = :(t)
        dict_arg = :($(tp) => $(set))
        bounds = :(ParameterBounds(($(dict_arg),)))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:call)) == (nothing, bounds)
    end
    # test _extract_bounds (:tuple)
    @testset "_extract_bounds (:tuple)" begin
        args = [:(t in [0, 1]), :(0 <= t <= 1)]
        set = :(IntervalSet(0, 1))
        tp = :(t)
        dict_arg = :($(tp) => $(set))
        bounds = :(ParameterBounds(($(dict_arg), $(dict_arg))))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:tuple)) == (nothing, bounds)
    end
    # test _extract_bounds (:comparison)
    @testset "_extract_bounds (:comparison)" begin
        args = [0, :(<=), :t, :(<=), 1]
        set = :(IntervalSet(0, 1))
        tp = :(t)
        dict_arg = :($(tp) => $(set))
        bounds = :(ParameterBounds(($(dict_arg),)))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:comparison)) == (nothing, bounds)
    end
    # test @finite_variable with bounds
    @testset "@finite_variable (Bounded)" begin
        # prepare parameters
        @infinite_parameter(m, 0 <= par <= 10)
        @infinite_parameter(m, 0 <= pars[1:2] <= 10)
        # test regular call
        idx = FiniteVariableIndex(5)
        vref = FiniteVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @finite_variable(m, xb >= 1, Bin,
                             parameter_bounds = (par == 0)) == gvref
        @test name(vref) == "xb"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        @test m.has_finite_var_bounds
        @test parameter_bounds(vref) == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        @test supports(par) == [0]
        # test regular tuple
        idx = FiniteVariableIndex(6)
        vref = FiniteVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @finite_variable(m, yb, parameter_bounds = (pars == 0,
                                                       par in [0, 1])) == gvref
        @test name(vref) == "yb"
        dict = Dict(pars[1] => IntervalSet(0, 0), pars[2] => IntervalSet(0, 0),
                    par => IntervalSet(0, 1))
        @test parameter_bounds(vref) == ParameterBounds(dict)
        @test supports(pars[1]) == [0]
        # test comparison
        idx = FiniteVariableIndex(7)
        vref = FiniteVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @finite_variable(m, zb, parameter_bounds = (0 <= par <= 1),
                             base_name = "bob") == gvref
        @test name(vref) == "bob"
        @test parameter_bounds(vref) == ParameterBounds(Dict(par => IntervalSet(0, 1)))
        # test unrecognized format
        @test_macro_throws ErrorException @finite_variable(m, parameter_bounds = par)
        # test container specification
        idxs = [FiniteVariableIndex(8), FiniteVariableIndex(9)]
        vrefs = [FiniteVariableRef(m, idx) for idx in idxs]
        vals = [0, 1]
        @test @finite_variable(m, [i = 1:2], parameter_bounds = (pars == vals[i]),
                             container = SparseAxisArray) isa JuMPC.SparseAxisArray
        @test name(vrefs[1]) == ""
        @test parameter_bounds(vrefs[1]) == ParameterBounds((pars => IntervalSet(0, 0),))
        # test wrong model type
        @test_macro_throws ErrorException @finite_variable(Model(),
                                                  parameter_bounds = (par == 0))
    end
    # test deprecation 
    @testset "@hold_variable" begin
        m = InfiniteModel()
        warn = ("`@hold_variable` in deprecated in favor of `@finite_variable` " * 
                "and will be dropped from future verions and may now behave unexpectedly.")
        @test_logs (:warn, warn) @hold_variable(m, z, binary = true)
    end
end

# test usage methods
@testset "Usage" begin
    # initialize model and stuff
    m = InfiniteModel()
    @finite_variable(m, y0)
    vref = dispatch_variable_ref(y0)
    # test used_by_measure
    @testset "used_by_measure" begin
        @test !used_by_measure(vref)
        push!(InfiniteOpt._measure_dependencies(vref), MeasureIndex(1))
        @test used_by_measure(y0)
        @test used_by_measure(vref)
        empty!(InfiniteOpt._measure_dependencies(vref))
    end
    # test used_by_constraint
    @testset "used_by_constraint" begin
        @test !used_by_constraint(vref)
        push!(InfiniteOpt._constraint_dependencies(vref), ConstraintIndex(1))
        @test used_by_constraint(y0)
        @test used_by_constraint(vref)
        empty!(InfiniteOpt._constraint_dependencies(vref))
    end
    # test used_by_objective
    @testset "used_by_objective" begin
        @test !used_by_objective(y0)
        @test !used_by_objective(vref)
        InfiniteOpt._data_object(vref).in_objective = true
        @test used_by_objective(vref)
        InfiniteOpt._data_object(vref).in_objective = false
    end
    # test is_used
    @testset "is_used" begin
        # test not used
        @test !is_used(vref)
        # test used by constraint and/or measure
        push!(InfiniteOpt._constraint_dependencies(vref), ConstraintIndex(1))
        @test is_used(y0)
        empty!(InfiniteOpt._constraint_dependencies(vref))
        # test used by objective
        InfiniteOpt._data_object(vref).in_objective = true
        @test is_used(vref)
    end
end

# Test parameter bound queries for finite variables
@testset "Parameter Bounds" begin
    # initialize the model and other needed information
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 10])
    @infinite_parameter(m, pars[1:2] in [0, 10])
    @infinite_variable(m, inf(par))
    dinf = dispatch_variable_ref(inf)
    @finite_variable(m, x == 0)
    dx = dispatch_variable_ref(x)
    @finite_variable(m, y >= 0, parameter_bounds = par == 0)
    dy = dispatch_variable_ref(y)
    data = TestData(par, 0, 5)
    meas = Measure(x + inf, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    dmref = MeasureRef(m, mindex)
    mref = InfiniteOpt._make_variable_ref(m, mindex)
    push!(InfiniteOpt._measure_dependencies(dx), mindex)
    # test parameter_bounds
    @testset "parameter_bounds" begin
        @test parameter_bounds(dx) == ParameterBounds()
        @test parameter_bounds(y) == ParameterBounds((par => IntervalSet(0, 0),))
    end
    # test has_parameter_bounds
    @testset "has_parameter_bounds" begin
        @test !has_parameter_bounds(x)
        @test has_parameter_bounds(dy)
        @test !has_parameter_bounds(dinf)
    end
    # test _update_variable_param_bounds
    @testset "_update_variable_param_bounds" begin
        bounds = ParameterBounds((par => IntervalSet(0, 1),))
        @test isa(InfiniteOpt._update_variable_param_bounds(dy, bounds), Nothing)
        @test parameter_bounds(dy) == bounds
    end
    # test _check_meas_bounds
    @testset "_check_meas_bounds" begin
        # test ok
        bounds = ParameterBounds((par => IntervalSet(0, 5),))
        @test isa(InfiniteOpt._check_meas_bounds(bounds, data), Nothing)
        # test lower bad
        bounds = ParameterBounds((par => IntervalSet(2, 5),))
        @test_throws ErrorException InfiniteOpt._check_meas_bounds(bounds, data)
        # test upper bad
        bounds = ParameterBounds((par => IntervalSet(0, 2),))
        @test_throws ErrorException InfiniteOpt._check_meas_bounds(bounds, data)
    end
    # test _check_meas_bounds (Fallback)
    @testset "_check_meas_bounds (Fallback)" begin
        warn = "Unable to check if finite variables bounds are valid in measure " *
               "with measure data type `BadData`. This can be resolved by " *
               "extending `measure_data_in_finite_var_bounds`."
        @test_logs (:warn, warn) InfiniteOpt._check_meas_bounds(ParameterBounds(),
                                                                BadData())
    end
    # test _update_bounds
    @testset "_update_bounds" begin
        # test normal
        bounds1 = ParameterBounds((pars[1] => IntervalSet(0, 1),
                                   pars[2] => IntervalSet(-1, 2)))
        bounds2 = ParameterBounds((par => IntervalSet(0, 0),
                                   pars[1] => IntervalSet(0.5, 1.5),
                                   pars[2] => IntervalSet(0, 1)))
        @test isa(InfiniteOpt._update_bounds(bounds1, bounds2), Nothing)
        @test length(bounds1) == 3
        @test bounds1[par] == IntervalSet(0, 0)
        @test bounds1[pars[1]] == IntervalSet(0.5, 1)
        @test bounds1[pars[2]] == IntervalSet(0, 1)
        # test error
        new_bounds = ParameterBounds((par => IntervalSet(1, 1),))
        @test_throws ErrorException InfiniteOpt._update_bounds(bounds1, new_bounds)
        new_bounds = ParameterBounds((par => IntervalSet(-1, -1),))
        @test_throws ErrorException InfiniteOpt._update_bounds(bounds1, new_bounds)
    end
    # test _update_var_bounds (DispatchVariableRef)
    @testset "_update_var_bounds (DispatchVariableRef)" begin
        @test isa(InfiniteOpt._update_var_bounds(dinf, ParameterBounds()),
                  Nothing)
    end
    # test _update_var_bounds (FiniteVariableRef)
    @testset "_update_var_bounds (FiniteVariableRef)" begin
        bounds = ParameterBounds()
        @test InfiniteOpt._update_variable_param_bounds(dx,
                      ParameterBounds((par => IntervalSet(0, 6),))) isa Nothing
        @test isa(InfiniteOpt._update_var_bounds(x, bounds), Nothing)
        @test bounds[par] == IntervalSet(0, 6)
    end
    # test _update_var_bounds (MeasureRef)
    @testset "_update_var_bounds (MeasureRef)" begin
        bounds = ParameterBounds()
        @test isa(InfiniteOpt._update_var_bounds(dmref, bounds), Nothing)
        @test bounds[par] == IntervalSet(0, 6)
        @test InfiniteOpt._update_variable_param_bounds(dx, ParameterBounds()) isa Nothing
    end
    # test _update_var_bounds (GeneralVariableRef)
    @testset "_update_var_bounds (GeneralVariableRef)" begin
        @test isa(InfiniteOpt._update_var_bounds(inf, ParameterBounds()),
                  Nothing)
    end
    # test _rebuild_constr_bounds (BoundedScalarConstraint)
    @testset "_rebuild_constr_bounds (Bounded)" begin
        # test addition
        bounds = ParameterBounds((par => IntervalSet(0, 2),))
        c = BoundedScalarConstraint(y + inf, MOI.EqualTo(0.), copy(bounds),
                                    copy(bounds))
        expected = ParameterBounds((par => IntervalSet(0, 1),))
        @test InfiniteOpt._rebuild_constr_bounds(c, bounds).bounds == expected
        # test deletion
        bounds = ParameterBounds((par => IntervalSet(0, 0),))
        c = BoundedScalarConstraint(y + inf, MOI.EqualTo(0.), bounds,
                                    ParameterBounds())
        InfiniteOpt._update_variable_param_bounds(dy, ParameterBounds())
        @test isa(InfiniteOpt._rebuild_constr_bounds(c, bounds),
                  ScalarConstraint)
        InfiniteOpt._update_variable_param_bounds(dy, bounds)
    end
    # test _rebuild_constr_bounds (ScalarConstraint)
    @testset "_rebuild_constr_bounds (Scalar)" begin
        bounds = ParameterBounds((par => IntervalSet(0, 0),))
        c = ScalarConstraint(y + inf, MOI.EqualTo(0.))
        @test InfiniteOpt._rebuild_constr_bounds(c, bounds).bounds == bounds
    end
    # test set_parameter_bounds
    @testset "set_parameter_bounds" begin
        bounds1 = ParameterBounds((pars[1] => IntervalSet(0, 1),))
        # test force error
        @test_throws ErrorException set_parameter_bounds(y, bounds1)
        # prepare for main test
        cref1 = FixRef(dx)
        bconstr = BoundedScalarConstraint(mref, MOI.EqualTo(0.0), bounds1, bounds1)
        cref2 = InfOptConstraintRef(m, ConstraintIndex(3))
        @test add_constraint(m, bconstr) == cref2
        # test measure error
        bounds = ParameterBounds((par => IntervalSet(1, 2),))
        @test_throws ErrorException set_parameter_bounds(dx, bounds)
        # test constr error
        bounds = ParameterBounds((pars[1] => IntervalSet(2, 3), ))
        @test_throws ErrorException set_parameter_bounds(dx, bounds)
        @test InfiniteOpt._update_variable_param_bounds(dx, ParameterBounds()) isa Nothing
        # test normal
        bounds = ParameterBounds((par => IntervalSet(0, 5),
                                  pars[1] => IntervalSet(0.5, 2)))
        @test isa(set_parameter_bounds(dx, bounds), Nothing)
        @test parameter_bounds(x) == bounds
        @test m.has_finite_var_bounds
        @test !optimizer_model_ready(m)
        expected = ParameterBounds((par => IntervalSet(0, 5),
                                    pars[1] => IntervalSet(0.5, 1)))
        @test parameter_bounds(cref2) == expected
        @test parameter_bounds(cref1) == bounds
        # test forced
        bounds = ParameterBounds((par => IntervalSet(1, 1),))
        @test isa(set_parameter_bounds(y, bounds, force = true), Nothing)
        @test parameter_bounds(dy) == bounds
        @test parameter_bounds(LowerBoundRef(dy)) == bounds
    end
    # test _update_constr_bounds (BoundedScalarConstraint)
    @testset "_update_constr_bounds (Bounded)" begin
        dict1 = Dict(pars[1] => IntervalSet(0, 1))
        dict2 = Dict(par => IntervalSet(0, 0), pars[1] => IntervalSet(1.1, 1.5))
        constr = BoundedScalarConstraint(par, MOI.EqualTo(0.0), ParameterBounds(copy(dict1)),
                                         ParameterBounds(copy(dict1)))
        # test error
        @test_throws ErrorException InfiniteOpt._update_constr_bounds(ParameterBounds(dict2), constr)
        @test constr.bounds == ParameterBounds(dict1)
        # test normal
        dict2 = Dict(par => IntervalSet(0, 0), pars[1] => IntervalSet(0.5, 1.5))
        dict = Dict(par => IntervalSet(0, 0), pars[1] => IntervalSet(0.5, 1))
        @test InfiniteOpt._update_constr_bounds(ParameterBounds(dict2), constr).bounds == ParameterBounds(dict)
    end
    # test _update_constr_bounds (ScalarConstraint)
    @testset "_update_constr_bounds (Scalar)" begin
        bounds = ParameterBounds((pars[1] => IntervalSet(0, 1),))
        constr = JuMP.ScalarConstraint(par, MOI.EqualTo(0.0))
        @test InfiniteOpt._update_constr_bounds(bounds, constr).bounds == bounds
    end
    # test add_parameter_bound
    @testset "add_parameter_bounds" begin
        cref1 = FixRef(dx)
        cref2 = InfOptConstraintRef(m, ConstraintIndex(3))
        # test adding normally
        new_bounds = ParameterBounds((pars[2] => IntervalSet(0, 5),))
        @test isa(add_parameter_bounds(x, new_bounds), Nothing)
        bounds = ParameterBounds((par => IntervalSet(0, 5),
                                  pars[1] => IntervalSet(0.5, 2),
                                  pars[2] => IntervalSet(0, 5)))
        @test parameter_bounds(dx) == bounds
        expected = ParameterBounds((par => IntervalSet(0, 5),
                                    pars[1] => IntervalSet(0.5, 1),
                                    pars[2] => IntervalSet(0, 5)))
        @test parameter_bounds(cref2) == expected
        @test parameter_bounds(cref1) == bounds
        # test measure error --> these errors are properly tested in the macros
        new_bounds = ParameterBounds((par => IntervalSet(1, 2),))
        @test_throws ErrorException add_parameter_bounds(dx, new_bounds)
        # test constr error
        new_bounds = ParameterBounds((pars[1] => IntervalSet(2, 3),))
        @test_throws ErrorException add_parameter_bounds(x, new_bounds)
        # test overwritting existing
        new_bounds = ParameterBounds((par => IntervalSet(0, 0),))
        @test isa(add_parameter_bounds(y, new_bounds), Nothing)
        @test parameter_bounds(y) == new_bounds
    end
    # test delete_parameter_bounds
    @testset "delete_parameter_bounds" begin
        cref1 = FixRef(dx)
        cref2 = InfOptConstraintRef(m, ConstraintIndex(3))
        # test normal
        @test isa(delete_parameter_bounds(x), Nothing)
        @test parameter_bounds(dx) == ParameterBounds()
        expected = ParameterBounds((pars[1] => IntervalSet(0, 1),))
        @test parameter_bounds(cref2) == expected
        @test !has_parameter_bounds(cref1)
        # test already gone
        @test isa(delete_parameter_bounds(dx), Nothing)
        @test parameter_bounds(dx) == ParameterBounds()
        @test parameter_bounds(cref2) == expected
        @test !has_parameter_bounds(cref1)
    end
end

# Test Parameter Bound Macros
@testset "Parameter Bound Macros" begin
    # initialize the model and other needed information
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 10])
    @infinite_parameter(m, pars[1:2] in [0, 10])
    @finite_variable(m, x == 0)
    @finite_variable(m, y, parameter_bounds = par == 0)
    data = TestData(par, 0, 5)
    meas = Measure(x + par, data, Int[], Int[], true)
    object = MeasureData(meas, "test")
    mindex = MeasureIndex(1)
    @test InfiniteOpt._add_data_object(m, object) == mindex
    dmref = MeasureRef(m, mindex)
    mref = InfiniteOpt._make_variable_ref(m, mindex)
    push!(InfiniteOpt._measure_dependencies(x), mindex)
    # test @add_parameter_bounds
    @testset "@add_parameter_bounds" begin
        bounds1 = ParameterBounds(Dict(pars[1] => IntervalSet(0, 1)))
        # prepare for main test
        cref1 = FixRef(x)
        bconstr = BoundedScalarConstraint(mref, MOI.EqualTo(0.0), bounds1, bounds1)
        cref2 = InfOptConstraintRef(m, ConstraintIndex(2))
        @test add_constraint(m, bconstr) == cref2
        # test adding normally
        @test isa(@add_parameter_bounds(x, pars[2] in [0, 5]), Nothing)
        dict = Dict(pars[2] => IntervalSet(0, 5))
        @test parameter_bounds(x) == ParameterBounds(dict)
        expected = Dict(pars[1] => IntervalSet(0, 1),
                        pars[2] => IntervalSet(0, 5))
        @test parameter_bounds(cref2) == ParameterBounds(expected)
        @test parameter_bounds(cref1) == ParameterBounds(dict)
        # test measure error
        @test_macro_throws ErrorException @add_parameter_bounds(x, par in [1, 2])
        # test constr error
        @test_macro_throws ErrorException @add_parameter_bounds(x, pars[1] in [2, 3])
        # test bad input format
        @test_macro_throws ErrorException @add_parameter_bounds(x, pars[1] = 0)
        @test_macro_throws ErrorException @add_parameter_bounds(par, pars[1] == 0)
        # test overwritting existing
        @test isa(@add_parameter_bounds(y, par == 1), Nothing)
        @test parameter_bounds(y) == ParameterBounds((par => IntervalSet(1, 1),))
    end
    # test @set_parameter_bounds
    @testset "@set_parameter_bounds" begin
        @test delete_parameter_bounds(x) isa Nothing
        bounds1 = ParameterBounds((pars[1] => IntervalSet(0, 1),))
        # test force error
        @test_macro_throws ErrorException @set_parameter_bounds(y, pars[1] in [0, 1])
        # test bad arguments
        @test_macro_throws ErrorException @set_parameter_bounds(x, x, 23)
        @test_macro_throws ErrorException @set_parameter_bounds(x, pars[1] == 0,
                                                                bad = 42)
        # test bad input format
        @test_macro_throws ErrorException @set_parameter_bounds(x, pars[1] = 0)
        @test_macro_throws ErrorException @set_parameter_bounds(par, pars[1] == 0)
        # prepare for main test
        cref1 = FixRef(x)
        cref2 = InfOptConstraintRef(m, ConstraintIndex(2))
        # test measure error
        @test_macro_throws ErrorException @set_parameter_bounds(x, par in [1, 2])
        # test constr error
        @test_macro_throws ErrorException @set_parameter_bounds(x, pars[1] in [2, 3])
        # test normal
        bounds = ParameterBounds((par => IntervalSet(0, 5),
                                      pars[1] => IntervalSet(0.5, 2)))
        @test isa(@set_parameter_bounds(x, (par in [0, 5], pars[1] in [0.5, 2])),
                  Nothing)
        @test parameter_bounds(x) == bounds
        @test m.has_finite_var_bounds
        @test !optimizer_model_ready(m)
        expected = ParameterBounds((par => IntervalSet(0, 5),
                                    pars[1] => IntervalSet(0.5, 1)))
        @test parameter_bounds(cref2) == expected
        @test parameter_bounds(cref1) == bounds
        # test forced
        @test isa(@set_parameter_bounds(y, par == 0, force = true), Nothing)
        @test parameter_bounds(y) == ParameterBounds((par => IntervalSet(0, 0),))
    end
end
