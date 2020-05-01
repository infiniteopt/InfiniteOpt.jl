# Test name methods
@testset "Hold Variable Name" begin
    # initialize model and variable
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    bounds = ParameterBounds()
    var = HoldVariable(info, bounds)
    m.vars[1] = var
    m.var_to_name[1] = "test"
    vref = HoldVariableRef(m, 1)
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "test"
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
    end
    # _make_variable_ref
    @testset "_make_variable_ref" begin
        @test InfiniteOpt._make_variable_ref(m, 1) == vref
    end
    # parameter_by_name
    @testset "JuMP.variable_by_name" begin
        # test normal
        @test variable_by_name(m, "new") == vref
        @test isa(variable_by_name(m, "test2"), Nothing)
        # prepare variable with duplicate name
        m.vars[2] = var
        m.var_to_name[2] = "new"
        m.name_to_var = nothing
        # test for duplciate name error
        @test_throws ErrorException variable_by_name(m, "new")
    end
end

# Test variable definition methods
@testset "Hold Variable Definition" begin
    # initialize model and info
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    info2 = VariableInfo(true, 0, true, 0, true, 0, true, 0, true, false)
    info3 = VariableInfo(true, 0, true, 0, true, 0, true, 0, false, true)
    bounds = ParameterBounds()
    @infinite_parameter(m, 0 <= par <= 10)
    @infinite_parameter(m, 0 <= pars[1:2] <= 10)
    # test _check_bounds
    @testset "_check_bounds" begin
        # test normal
        @test isa(InfiniteOpt._check_bounds(ParameterBounds(Dict(par => IntervalSet(0, 1)))),
                                                                        Nothing)
        # test errors
        @test_throws ErrorException InfiniteOpt._check_bounds(
                                ParameterBounds(Dict(par => IntervalSet(-1, 1))))
        @test_throws ErrorException InfiniteOpt._check_bounds(
                                ParameterBounds(Dict(par => IntervalSet(0, 11))))
    end
    # _make_variable
    @testset "_make_variable" begin
        # test normal
        expected = HoldVariable(info, bounds)
        @test InfiniteOpt._make_variable(error, info, Val(Hold)).info == expected.info
        # test errors
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                                Val(Hold), parameter_values = 3)
    end
    # build_variable
    @testset "JuMP.build_variable" begin
        # test normal
        expected = HoldVariable(info, bounds)
        @test build_variable(error, info, Hold).info == expected.info
        # test errors
        @test_throws ErrorException build_variable(error, info, Point,
                                                   parameter_bounds = bounds)
    end
    # _validate_bounds
    @testset "_validate_bounds" begin
        # test normal
        @test isa(InfiniteOpt._validate_bounds(m,
                      ParameterBounds(Dict(par => IntervalSet(0, 1)))), Nothing)
        # test error
        par2 = ParameterRef(InfiniteModel(), 1)
        @test_throws ErrorException InfiniteOpt._validate_bounds(m,
                               ParameterBounds(Dict(par2 => IntervalSet(0, 1))))
        # test support addition
        @test isa(InfiniteOpt._validate_bounds(m,
                      ParameterBounds(Dict(par => IntervalSet(0, 0)))), Nothing)
        @test supports(par) == [0]
    end
    # _check_and_make_variable_ref
    @testset "_check_and_make_variable_ref" begin
        # test normal
        v = build_variable(error, info, Hold)
        @test InfiniteOpt._check_and_make_variable_ref(m, v) == HoldVariableRef(m, 0)
        # test with bounds
        bounds = ParameterBounds(Dict(par => IntervalSet(0, 2)))
        v = build_variable(error, info, Hold, parameter_bounds = bounds)
        @test InfiniteOpt._check_and_make_variable_ref(m, v) == HoldVariableRef(m, 0)
        @test m.has_hold_bounds
        m.has_hold_bounds = false
        # test bad bounds
        @infinite_parameter(InfiniteModel(), par2 in [0, 2])
        v = build_variable(error, info, Hold,
            parameter_bounds = ParameterBounds(Dict(par2 => IntervalSet(0, 1))))
        @test_throws ErrorException InfiniteOpt._check_and_make_variable_ref(m, v)
    end
    # add_variable
    @testset "JuMP.add_variable" begin
        v = build_variable(error, info, Hold)
        @test add_variable(m, v, "name") == HoldVariableRef(m, 1)
        @test haskey(m.vars, 1)
        @test m.var_to_name[1] == "name"
        # prepare infinite variable with all the possible info additions
        v = build_variable(error, info2, Hold)
        # test info addition functions
        vref = HoldVariableRef(m, 2)
        @test add_variable(m, v, "name") == vref
        @test !optimizer_model_ready(m)
        # lower bound
        @test has_lower_bound(vref)
        @test JuMP._lower_bound_index(vref) == 1
        @test isa(m.constrs[1], ScalarConstraint{HoldVariableRef,
                                                 MOI.GreaterThan{Float64}})
        @test m.constr_in_var_info[1]
        # upper bound
        @test has_upper_bound(vref)
        @test JuMP._upper_bound_index(vref) == 2
        @test isa(m.constrs[2], ScalarConstraint{HoldVariableRef,
                                                 MOI.LessThan{Float64}})
        @test m.constr_in_var_info[2]
        # fix
        @test is_fixed(vref)
        @test JuMP._fix_index(vref) == 3
        @test isa(m.constrs[3], ScalarConstraint{HoldVariableRef,
                                                 MOI.EqualTo{Float64}})
        @test m.constr_in_var_info[3]
        # binary
        @test is_binary(vref)
        @test JuMP._binary_index(vref) == 4
        @test isa(m.constrs[4], ScalarConstraint{HoldVariableRef,
                                                 MOI.ZeroOne})
        @test m.constr_in_var_info[4]
        @test m.var_to_constrs[2] == [1, 2, 3, 4]
        # prepare infinite variable with integer info addition
        v = build_variable(error, info3, Hold)
        # test integer addition functions
        vref = HoldVariableRef(m, 3)
        @test add_variable(m, v, "name") == vref
        @test !optimizer_model_ready(m)
        @test is_integer(vref)
        @test JuMP._integer_index(vref) == 8
        @test isa(m.constrs[8], ScalarConstraint{HoldVariableRef,
                                                 MOI.Integer})
        @test m.constr_in_var_info[8]
        @test m.var_to_constrs[3] == [5, 6, 7, 8]
    end
end

# Test the hold variable macro if no bounds are provided
@testset "Hold" begin
    # initialize model
    m = InfiniteModel()
    @testset "@hold_variable (Not Bounded)" begin
        # test regular
        vref = HoldVariableRef(m, 1)
        @test @hold_variable(m, x >= 1, Bin) == vref
        @test name(vref) == "x"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        # test anan
        vref = HoldVariableRef(m, 2)
        @test @hold_variable(m, binary = true, lower_bound = 1,
                               base_name = "x") == vref
        @test name(vref) == "x"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        # test array
        vrefs = [HoldVariableRef(m, 3), HoldVariableRef(m, 4)]
        @test @hold_variable(m, y[1:2] == 2, Int) == vrefs
        @test name(vrefs[1]) == "y[1]"
        @test fix_value(vrefs[2]) == 2
        @test is_integer(vrefs[1])
        # test errors
        @test_throws AssertionError @hold_variable(Model(), z >= 1, Bin)
        @test_macro_throws ErrorException @hold_variable(m, x >= 1, Bin)
    end
    # test _make_interval_set
    @testset "_make_interval_set" begin
        expected = :(IntervalSet(0, 1))
        @test InfiniteOpt._make_interval_set(error, :([0, 1])) == expected
        @test_throws ErrorException InfiniteOpt._make_interval_set(error, :([0]))
    end
    # test _process_parameter
    @testset "_process_parameter" begin
        @test InfiniteOpt._process_parameter(:x) == :(InfiniteOpt._make_vector(x))
    end
    # test _make_bound_pair (in)
    @testset "_make_bound_pair (in)" begin
        expr = :(t in [0, 1])
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr,
                                           Val(:in)) == :($(tp) => $(set))
    end
    # test _make_bound_pair (==)
    @testset "_make_bound_pair (==)" begin
        expr = :(t == 0)
        set = :(IntervalSet(0, 0))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr,
                                           Val(:(==))) == :($(tp) => $(set))
    end
    # test _make_bound_pair (comparison with <=)
    @testset "_make_bound_pair (<=)" begin
        expr = :(0 <= t <= 1)
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        @test InfiniteOpt._make_bound_pair(error, expr, Val(:(<=)),
                                           Val(:(<=))) == :($(tp) => $(set))
    end
    # test _make_bound_pair (comparison with >=)
    @testset "_make_bound_pair (>=)" begin
        expr = :(1 >= t >= 0)
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
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
        tp = :(InfiniteOpt._make_vector(t))
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
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        dict = :(Dict($(dict_arg), $(dict_arg)))
        bounds = :(ParameterBounds($dict))
        @test InfiniteOpt._parse_parameter_bounds(error, args) == bounds
    end
    # test _parse_parameter_bounds (Expr)
    @testset "_parse_parameter_bounds (Expr)" begin
        expr = :(t in [0, 1])
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        dict = :(Dict($(dict_arg)))
        bounds = :(ParameterBounds($dict))
        @test InfiniteOpt._parse_parameter_bounds(error, expr) == bounds
    end
    # test _extract_bounds (:call)
    @testset "_extract_bounds (:call)" begin
        # test single anonymous bound
        args = [:in, :t, :([0, 1])]
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        bounds = :(ParameterBounds(Dict($(dict_arg))))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:call)) == (nothing, bounds)
    end
    # test _extract_bounds (:tuple)
    @testset "_extract_bounds (:tuple)" begin
        args = [:(t in [0, 1]), :(0 <= t <= 1)]
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        bounds = :(ParameterBounds(Dict($(dict_arg), $(dict_arg))))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:tuple)) == (nothing, bounds)
    end
    # test _extract_bounds (:comparison)
    @testset "_extract_bounds (:comparison)" begin
        args = [0, :(<=), :t, :(<=), 1]
        set = :(IntervalSet(0, 1))
        tp = :(InfiniteOpt._make_vector(t))
        dict_arg = :($(tp) => $(set))
        bounds = :(ParameterBounds(Dict($(dict_arg))))
        @test InfiniteOpt._extract_bounds(error, args,
                                          Val(:comparison)) == (nothing, bounds)
    end
    # test @hold_variable with bounds
    @testset "@hold_variable (Bounded)" begin
        # prepare parameters
        @infinite_parameter(m, 0 <= par <= 10)
        @infinite_parameter(m, 0 <= pars[1:2] <= 10)
        # test regular call
        vref = HoldVariableRef(m, 5)
        @test @hold_variable(m, xb >= 1, Bin,
                             parameter_bounds = (par == 0)) == vref
        @test name(vref) == "xb"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        @test m.has_hold_bounds
        @test m.vars[5].parameter_bounds == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        @test supports(par) == [0]
        # test regular tuple
        vref = HoldVariableRef(m, 6)
        @test @hold_variable(m, yb, parameter_bounds = (pars == 0,
                                                       par in [0, 1])) == vref
        @test name(vref) == "yb"
        dict = Dict(pars[1] => IntervalSet(0, 0), pars[2] => IntervalSet(0, 0),
                    par => IntervalSet(0, 1))
        @test m.vars[6].parameter_bounds == ParameterBounds(dict)
        @test supports(pars[1]) == [0]
        # test comparison
        vref = HoldVariableRef(m, 7)
        @test @hold_variable(m, zb, parameter_bounds = (0 <= par <= 1),
                             base_name = "bob") == vref
        @test name(vref) == "bob"
        @test m.vars[7].parameter_bounds == ParameterBounds(Dict(par => IntervalSet(0, 1)))
        # test unrecognized format
        @test_macro_throws ErrorException @hold_variable(m, parameter_bounds = par)
        # test container specification
        vrefs = [HoldVariableRef(m, 8), HoldVariableRef(m, 9)]
        vrefs = convert(JuMPC.SparseAxisArray, vrefs)
        @test @hold_variable(m, [i = 1:2], parameter_bounds = (pars[i] == 0),
                             container = SparseAxisArray) == vrefs
        @test name(vrefs[1]) == ""
        @test m.vars[9].parameter_bounds == ParameterBounds(Dict(pars[2] => IntervalSet(0, 0)))
        # test wrong model type
        @test_macro_throws ErrorException @hold_variable(Model(),
                                                  parameter_bounds = (par == 0))
    end
end
