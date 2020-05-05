# Test name methods
@testset "Basics" begin
    # initialize model and variable
    m = InfiniteModel()
    num = Float64(0)
    info = VariableInfo(false, num, false, num, false, num, false, num, false, false)
    new_info = VariableInfo(true, 0., true, 0., true, 0., true, 0., true, false)
    bounds = ParameterBounds()
    var = HoldVariable(info, bounds)
    object = VariableData(var, "test")
    idx = HoldVariableIndex(1)
    vref = HoldVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, HoldVariableIndex)
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
        @test InfiniteOpt._data_dictionary(m, HoldVariable) === m.hold_vars
        @test InfiniteOpt._data_dictionary(vref) === m.hold_vars
        @test InfiniteOpt._data_dictionary(gvref) === m.hold_vars
    end
    # JuMP.is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, vref)
        @test is_valid(m, gvref)
    end
    # _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(vref) == object
        @test InfiniteOpt._data_object(gvref) == object
    end
    # _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(vref) == var
        @test InfiniteOpt._core_variable_object(gvref) == var
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
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
        @test isa(set_name(gvref, "new2"), Nothing)
        @test name(vref) == "new2"
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
        idx2 = HoldVariableIndex(2)
        @test InfiniteOpt._add_data_object(m, object) == idx2
        vref2 = HoldVariableRef(m, idx2)
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
        @test InfiniteOpt._make_variable(error, info, Val(Hold)).info == info
        bounds = ParameterBounds((pars => IntervalSet(0, 0),))
        @test InfiniteOpt._make_variable(error, info, Val(Hold),
                                         parameter_bounds = bounds).info == info
        @test InfiniteOpt._make_variable(error, info, Val(Hold),
                            parameter_bounds = bounds).parameter_bounds == bounds
        # test errors
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                                Val(Hold), parameter_values = 3)
        bounds = ParameterBounds((pars[1] => IntervalSet(0, 0),))
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                           Val(Hold), parameter_bounds = bounds)
        bounds = ParameterBounds((par => IntervalSet(-1, -1),))
        @test_throws ErrorException InfiniteOpt._make_variable(error, info,
                                           Val(Hold), parameter_bounds = bounds)
    end
    # build_variable
    @testset "JuMP.build_variable" begin
        # test normal
        @test build_variable(error, info, Hold).info == info
        @test build_variable(error, info, Hold).parameter_bounds == ParameterBounds()
        bounds = ParameterBounds((pars => IntervalSet(0, 5),))
        @test build_variable(error, info, Hold,
                             parameter_bounds = bounds).info == info
        @test build_variable(error, info, Hold,
                             parameter_bounds = bounds).parameter_bounds == bounds
        # test errors
        @test_throws ErrorException build_variable(error, info, Point,
                                                   parameter_bounds = bounds)
        bounds = ParameterBounds((par => IntervalSet(0, 0),
                                  pars => IntervalSet(-1, 1)))
        @test_throws ErrorException build_variable(error, info, Hold,
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
        v = build_variable(error, info, Hold)
        idx = HoldVariableIndex(1)
        vref = HoldVariableRef(m, idx)
        @test InfiniteOpt._check_and_make_variable_ref(m, v) == vref
        # test with bounds
        bounds = ParameterBounds((par => IntervalSet(0, 2),))
        v = build_variable(error, info, Hold, parameter_bounds = bounds)
        idx = HoldVariableIndex(2)
        vref = HoldVariableRef(m, idx)
        @test InfiniteOpt._check_and_make_variable_ref(m, v) == vref
        @test m.has_hold_bounds
        m.has_hold_bounds = false
        # test bad bounds
        @infinite_parameter(InfiniteModel(), par2 in [0, 2])
        v = build_variable(error, info, Hold,
            parameter_bounds = ParameterBounds((par2 => IntervalSet(0, 1),)))
        @test_throws VariableNotOwned{GeneralVariableRef} InfiniteOpt._check_and_make_variable_ref(m, v)
    end
    # add_variable
    @testset "JuMP.add_variable" begin
        idx = HoldVariableIndex(3)
        vref = HoldVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        v = build_variable(error, info, Hold)
        @test add_variable(m, v, "name") == gvref
        @test haskey(InfiniteOpt._data_dictionary(vref), idx)
        @test name(vref) == "name"
        # prepare infinite variable with all the possible info additions
        v = build_variable(error, info2, Hold)
        # test info addition functions
        idx = HoldVariableIndex(4)
        vref = HoldVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test add_variable(m, v, "name") == gvref
        @test !optimizer_model_ready(m)
        @test !m.has_hold_bounds
        # lower bound
        cindex = ConstraintIndex(1)
        cref = FiniteConstraintRef(m, cindex, ScalarShape())
        @test has_lower_bound(vref)
        @test JuMP._lower_bound_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.GreaterThan{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # upper bound
        cindex = ConstraintIndex(2)
        cref = FiniteConstraintRef(m, cindex, ScalarShape())
        @test has_upper_bound(vref)
        @test JuMP._upper_bound_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.LessThan{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # fix
        cindex = ConstraintIndex(3)
        cref = FiniteConstraintRef(m, cindex, ScalarShape())
        @test is_fixed(vref)
        @test JuMP._fix_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.EqualTo{Float64}}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        # binary
        cindex = ConstraintIndex(4)
        cref = FiniteConstraintRef(m, cindex, ScalarShape())
        @test is_binary(vref)
        @test JuMP._binary_index(vref) == cindex
        @test constraint_object(cref) isa ScalarConstraint{GeneralVariableRef,
                                                 MOI.ZeroOne}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(vref) == [ConstraintIndex(i)
                                                             for i = 1:4]
        # prepare infinite variable with integer info addition
        bounds = ParameterBounds((par => IntervalSet(0, 2),
                                  pars => IntervalSet(1, 1)))
        v = build_variable(error, info3, Hold, parameter_bounds = bounds)
        # test integer addition functions
        idx = HoldVariableIndex(5)
        vref = HoldVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test add_variable(m, v, "name") == gvref
        @test !optimizer_model_ready(m)
        @test m.has_hold_bounds
        @test supports(par) == [0]
        @test sortcols(supports(pars)) == [0 1; 0 1]
        cindex = ConstraintIndex(8)
        cref = FiniteConstraintRef(m, cindex, ScalarShape())
        @test is_integer(vref)
        @test JuMP._integer_index(vref) == cindex
        @test constraint_object(cref) isa BoundedScalarConstraint{GeneralVariableRef,
                                                            MOI.Integer}
        @test InfiniteOpt._data_object(cref).is_info_constraint
        @test InfiniteOpt._constraint_dependencies(vref) == [ConstraintIndex(i)
                                                             for i = 5:8]
    end
end

# Test the hold variable macro if no bounds are provided
@testset "Macro" begin
    # initialize model
    m = InfiniteModel()
    @testset "@hold_variable (Not Bounded)" begin
        # test regular
        idx = HoldVariableIndex(1)
        vref = HoldVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @hold_variable(m, x >= 1, Bin) == gvref
        @test name(vref) == "x"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        # test anan
        idx = HoldVariableIndex(2)
        vref = HoldVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @hold_variable(m, binary = true, lower_bound = 1,
                               base_name = "x") == gvref
        @test name(vref) == "x"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        # test array
        idxs = [HoldVariableIndex(3), HoldVariableIndex(4)]
        vrefs = [HoldVariableRef(m, idx) for idx in idxs]
        gvrefs = [InfiniteOpt._make_variable_ref(m, idx) for idx in idxs]
        @test @hold_variable(m, y[1:2] == 2, Int) == gvrefs
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
    # test @hold_variable with bounds
    @testset "@hold_variable (Bounded)" begin
        # prepare parameters
        @infinite_parameter(m, 0 <= par <= 10)
        @infinite_parameter(m, 0 <= pars[1:2] <= 10)
        # test regular call
        idx = HoldVariableIndex(5)
        vref = HoldVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @hold_variable(m, xb >= 1, Bin,
                             parameter_bounds = (par == 0)) == gvref
        @test name(vref) == "xb"
        @test lower_bound(vref) == 1
        @test is_binary(vref)
        @test m.has_hold_bounds
        @test parameter_bounds(vref) == ParameterBounds(Dict(par => IntervalSet(0, 0)))
        @test supports(par) == [0]
        # test regular tuple
        idx = HoldVariableIndex(6)
        vref = HoldVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @hold_variable(m, yb, parameter_bounds = (pars == 0,
                                                       par in [0, 1])) == gvref
        @test name(vref) == "yb"
        dict = Dict(pars[1] => IntervalSet(0, 0), pars[2] => IntervalSet(0, 0),
                    par => IntervalSet(0, 1))
        @test parameter_bounds(vref) == ParameterBounds(dict)
        @test supports(pars[1]) == [0]
        # test comparison
        idx = HoldVariableIndex(7)
        vref = HoldVariableRef(m, idx)
        gvref = InfiniteOpt._make_variable_ref(m, idx)
        @test @hold_variable(m, zb, parameter_bounds = (0 <= par <= 1),
                             base_name = "bob") == gvref
        @test name(vref) == "bob"
        @test parameter_bounds(vref) == ParameterBounds(Dict(par => IntervalSet(0, 1)))
        # test unrecognized format
        @test_macro_throws ErrorException @hold_variable(m, parameter_bounds = par)
        # test container specification
        idxs = [HoldVariableIndex(8), HoldVariableIndex(9)]
        vrefs = [HoldVariableRef(m, idx) for idx in idxs]
        vals = [0, 1]
        @test @hold_variable(m, [i = 1:2], parameter_bounds = (pars == vals[i]),
                             container = SparseAxisArray) isa JuMPC.SparseAxisArray
        @test name(vrefs[1]) == ""
        @test parameter_bounds(vrefs[1]) == ParameterBounds((pars => IntervalSet(0, 0),))
        # test wrong model type
        @test_macro_throws ErrorException @hold_variable(Model(),
                                                  parameter_bounds = (par == 0))
    end
end
