# Test macro methods
@testset "Macro Helpers" begin
    @testset "Symbol Methods" begin
        @test Parameter == :Parameter
        @test InfiniteOpt._is_set_keyword(:(lower_bound = 0))
    end
    # ParameterInfoExpr datatype
    @testset "_ParameterInfoExpr" begin
        @test InfiniteOpt._ParameterInfoExpr isa DataType
        @test InfiniteOpt._ParameterInfoExpr(ones(Bool, 8)...).has_lb
        @test !InfiniteOpt._ParameterInfoExpr().has_lb
    end
    # JuMP._set_lower_bound_or_error
    @testset "JuMP._set_lower_bound_or_error" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(JuMP._set_lower_bound_or_error(error, info, 0), Nothing)
        @test info.has_lb && info.lower_bound == 0
        @test_throws ErrorException JuMP._set_lower_bound_or_error(error,
                                                                   info, 0)
        info.has_lb = false; info.has_dist = true
        @test_throws ErrorException JuMP._set_lower_bound_or_error(error,
                                                                   info, 0)
        info.has_dist = false; info.has_set = true
        @test_throws ErrorException JuMP._set_lower_bound_or_error(error,
                                                                   info, 0)
    end
    # JuMP._set_upper_bound_or_error
    @testset "JuMP._set_upper_bound_or_error" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(JuMP._set_upper_bound_or_error(error, info, 0), Nothing)
        @test info.has_ub && info.upper_bound == 0
        @test_throws ErrorException JuMP._set_upper_bound_or_error(error,
                                                                   info, 0)
        info.has_ub = false; info.has_dist = true
        @test_throws ErrorException JuMP._set_upper_bound_or_error(error,
                                                                   info, 0)
        info.has_dist = false; info.has_set = true
        @test_throws ErrorException JuMP._set_upper_bound_or_error(error,
                                                                   info, 0)
    end
    # _dist_or_error
    @testset "_dist_or_error" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(InfiniteOpt._dist_or_error(error, info, 0), Nothing)
        @test info.has_dist && info.distribution == 0
        @test_throws ErrorException InfiniteOpt._dist_or_error(error, info, 0)
        info.has_dist = false; info.has_lb = true
        @test_throws ErrorException InfiniteOpt._dist_or_error(error, info, 0)
        info.has_lb = false; info.has_set = true
        @test_throws ErrorException InfiniteOpt._dist_or_error(error, info, 0)
    end
    # _set_or_error
    @testset "_set_or_error" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(InfiniteOpt._set_or_error(error, info, 0), Nothing)
        @test info.has_set && info.set == 0
        @test_throws ErrorException InfiniteOpt._set_or_error(error, info, 0)
        info.has_set = false; info.has_lb = true
        @test_throws ErrorException InfiniteOpt._set_or_error(error, info, 0)
        info.has_lb = false; info.has_dist = true
        @test_throws ErrorException InfiniteOpt._set_or_error(error, info, 0)
    end
    # _constructor_set
    @testset "InfiniteOpt._constructor_set" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test_throws ErrorException InfiniteOpt._constructor_set(error, info)
        info.has_lb = true; info.lower_bound = 0
        @test_throws ErrorException InfiniteOpt._constructor_set(error, info)
        info.has_ub = true; info.upper_bound = 1
        check = :(isa($(info.lower_bound), Number))
        expected = :($(check) ? IntervalSet($(info.lower_bound), $(info.upper_bound)) : error("Bounds must be a number."))
        @test InfiniteOpt._constructor_set(error, info) == expected
        info = InfiniteOpt._ParameterInfoExpr(distribution = Normal())
        check = :(isa($(info.distribution), Distributions.NonMatrixDistribution))
        expected = :($(check) ? DistributionSet($(info.distribution)) : error("Distribution must be a subtype of Distributions.NonMatrixDistribution."))
        @test InfiniteOpt._constructor_set(error, info) == expected
        info = InfiniteOpt._ParameterInfoExpr(set = IntervalSet(0, 1))
        check = :(isa($(info.set), AbstractInfiniteSet))
        expected = :($(check) ? $(info.set) : error("Set must be a subtype of AbstractInfiniteSet."))
        @test InfiniteOpt._constructor_set(error, info) == expected
    end
    # _parse_one_operator_parameter
    @testset "_parse_one_operator_parameter" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(InfiniteOpt._parse_one_operator_parameter(error, info,
                                                          Val(:<=), 0), Nothing)
        @test info.has_ub && info.upper_bound == 0
        @test isa(InfiniteOpt._parse_one_operator_parameter(error, info,
                                                          Val(:>=), 0), Nothing)
        @test info.has_lb && info.lower_bound == 0
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(InfiniteOpt._parse_one_operator_parameter(error, info,
                                                          Val(:in), 0), Nothing)
        @test info.has_dist && info.distribution == 0
        @test_throws ErrorException InfiniteOpt._parse_one_operator_parameter(error, info,
                                                                              Val(:d), 0)
    end
    # _parse_ternary_parameter
    @testset "_parse_ternary_parameter" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(InfiniteOpt._parse_ternary_parameter(error, info, Val(:<=), 0,
                                                       Val(:<=), 0), Nothing)
        @test info.has_ub && info.upper_bound == 0
        @test info.has_lb && info.lower_bound == 0
        info = InfiniteOpt._ParameterInfoExpr()
        @test isa(InfiniteOpt._parse_ternary_parameter(error, info, Val(:>=), 0,
                                                       Val(:>=), 0), Nothing)
        @test info.has_ub && info.upper_bound == 0
        @test info.has_lb && info.lower_bound == 0
        @test_throws ErrorException InfiniteOpt._parse_ternary_parameter(error,
                                                 info, Val(:<=), 0, Val(:>=), 0)
    end
    # _parse_parameter
    @testset "_parse_parameter" begin
        info = InfiniteOpt._ParameterInfoExpr()
        @test InfiniteOpt._parse_parameter(error, info, :<=, :x, 0) == :x
        @test info.has_ub && info.upper_bound == 0
        @test InfiniteOpt._parse_parameter(error, info, :<=, 0, :x) == :x
        @test info.has_lb && info.lower_bound == 0
        info = InfiniteOpt._ParameterInfoExpr()
        @test InfiniteOpt._parse_parameter(error, info, 0, :<=, :x, :<=, 0) == :x
        @test info.has_ub && info.upper_bound == 0
        @test info.has_lb && info.lower_bound == 0
    end
end

# Test precursor functions needed for add_parameter
@testset "Basic Reference Queries" begin
    m = InfiniteModel()
    pref = ParameterRef(m, 1)
    # JuMP.index
    @testset "JuMP.index" begin
        @test index(pref) == 1
    end
    # JuMP.owner_model
    @testset "JuMP.owner_model" begin
        @test owner_model(pref) == m
    end
end

# Test parameter definition methods
@testset "Definition" begin
    # _check_supports_in_bounds
    @testset "_check_supports_in_bounds" begin
        set = IntervalSet(0, 1)
        @test isa(InfiniteOpt._check_supports_in_bounds(error, 0, set), Nothing)
        @test_throws ErrorException InfiniteOpt._check_supports_in_bounds(error,
                                                                        -1, set)
        @test_throws ErrorException InfiniteOpt._check_supports_in_bounds(error,
                                                                         2, set)
        set = DistributionSet(Uniform())
        @test isa(InfiniteOpt._check_supports_in_bounds(error, 0, set), Nothing)
        @test_throws ErrorException InfiniteOpt._check_supports_in_bounds(error,
                                                                        -1, set)
        @test_throws ErrorException InfiniteOpt._check_supports_in_bounds(error,
                                                                         2, set)
    end
    # build_parameter
    @testset "build_parameter" begin
        set = DistributionSet(Multinomial(3, 2))
        supps = Vector(0:1)
        expected = InfOptParameter(set, supps, false)
        @test build_parameter(error, set, 2, supports = supps) == expected
        expected = InfOptParameter(set, Number[], true)
        @test build_parameter(error, set, 2, independent = true) == expected
        @test_throws ErrorException build_parameter(error, set, 2, bob = 42)
        @test_throws ErrorException build_parameter(error, set, 1)
        warn = "Support points are not unique, eliminating redundant points."
        @test_logs (:warn, warn) build_parameter(error, set, 2,
                                                 supports = ones(3))
    end
    # add_parameter
    @testset "add_parameter" begin
        m = InfiniteModel()
        param = InfOptParameter(IntervalSet(0, 1), Number[], false)
        expected = ParameterRef(m, 1)
        @test add_parameter(m, param) == expected
        @test m.params[1] isa InfOptParameter
        @test m.param_to_group_id[1] == 0
    end
end

# Test basic get/set methods
@testset "Name" begin
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(pref) == "test"
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(pref, "new"), Nothing)
        @test name(pref) == "new"
    end
    # parameter_by_name
    @testset "parameter_by_name" begin
        @test parameter_by_name(m, "new") == pref
        @test isa(parameter_by_name(m, "test2"), Nothing)
        pref2 = add_parameter(m, param, "new")
        @test_throws ErrorException parameter_by_name(m, "new")
    end
end

# Test the parameter macro
@testset "Macro" begin
    m = InfiniteModel()
    # single parameter
    @testset "Single" begin
        pref = ParameterRef(m, 1)
        @test @infinite_parameter(m, 0 <= a <= 1) == pref
        @test m.params[1].set == IntervalSet(0, 1)
        @test name(pref) == "a"
        pref = ParameterRef(m, 2)
        @test @infinite_parameter(m, b in Normal(), supports = [1; 2]) == pref
        @test m.params[2].set == DistributionSet(Normal())
        @test m.params[2].supports == [1, 2]
        pref = ParameterRef(m, 3)
        @test @infinite_parameter(m, c, set = IntervalSet(0, 1)) == pref
        pref = ParameterRef(m, 4)
        @test @infinite_parameter(m, set = IntervalSet(0, 1),
                                  base_name = "d") == pref
        @test name(pref) == "d"
        pref = ParameterRef(m, 5)
        @test @infinite_parameter(m, set = IntervalSet(0, 1)) == pref
        @test name(pref) == ""
    end
    # multiple parameters
    @testset "Array" begin
        prefs = [ParameterRef(m, 6), ParameterRef(m, 7)]
        @test @infinite_parameter(m, 0 <= e[1:2] <= 1) == prefs
        @test m.params[6].set == IntervalSet(0, 1)
        @test m.params[7].set == IntervalSet(0, 1)
        prefs = [ParameterRef(m, 8), ParameterRef(m, 9)]
        @test @infinite_parameter(m, [1:2], set = IntervalSet(0, 1)) == prefs
        @test m.params[8].set == IntervalSet(0, 1)
        @test m.params[9].set == IntervalSet(0, 1)
        prefs = [ParameterRef(m, 10), ParameterRef(m, 11)]
        sets = [IntervalSet(0, 1), IntervalSet(-1, 2)]
        @test @infinite_parameter(m, f[i = 1:2], set = sets[i]) == prefs
        @test m.params[10].set == IntervalSet(0, 1)
        @test m.params[11].set == IntervalSet(-1, 2)
        prefs = [ParameterRef(m, 12), ParameterRef(m, 13)]
        @test @infinite_parameter(m, [i = 1:2], set = sets[i]) == prefs
        @test m.params[12].set == IntervalSet(0, 1)
        @test m.params[13].set == IntervalSet(-1, 2)
        prefs = [ParameterRef(m, 14), ParameterRef(m, 15)]
        @test @infinite_parameter(m, [0, -1][i] <= g[i = 1:2] <= [1, 2][i]) == prefs
        @test m.params[14].set == IntervalSet(0, 1)
        @test m.params[15].set == IntervalSet(-1, 2)
        prefs = [ParameterRef(m, 16), ParameterRef(m, 17)]
        @test @infinite_parameter(m, 0 <= h[1:2] <= 1,
                                  independent = true) == prefs
        @test m.params[16].independent
        @test m.params[17].independent
        prefs = [ParameterRef(m, 18), ParameterRef(m, 19)]
        prefs = convert(JuMP.Containers.SparseAxisArray, prefs)
        @test @infinite_parameter(m, 0 <= i[1:2] <= 1,
                                  container = SparseAxisArray) == prefs
        @test m.params[18].set == IntervalSet(0, 1)
        @test m.params[19].set == IntervalSet(0, 1)
    end
    # test for errors
    @testset "Errors" begin
        @test_macro_throws ErrorException @infinite_parameter(m, 0 <= [1:2] <= 1)
        @test_macro_throws ErrorException @infinite_parameter(m, 0 <= "bob" <= 1)
        @test_macro_throws ErrorException @infinite_parameter(m, 0 <= a <= 1)
        @test_macro_throws ErrorException @infinite_parameter(m, 0 <= j)
        @test_macro_throws ErrorException @infinite_parameter(m, j)
        @test_macro_throws ErrorException @infinite_parameter(m, j, foo = 42)
        @test_macro_throws ErrorException @infinite_parameter(m, j in Multinomial(3, 2))
    end
end

# Test if used
@testset "Used" begin
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    # used_by_constraint
    @testset "used_by_constraint" begin
        @test !used_by_constraint(pref)
        m.param_to_constrs[index(pref)] = [1]
        @test used_by_constraint(pref)
        delete!(m.param_to_constrs, index(pref))
    end
    # used_by_measure
    @testset "used_by_measure" begin
        @test !used_by_measure(pref)
        m.param_to_meas[index(pref)] = [1]
        @test used_by_measure(pref)
        delete!(m.param_to_meas, index(pref))
    end
    # used_by_measure
    @testset "used_by_variable" begin
        @test !used_by_variable(pref)
        m.param_to_vars[index(pref)] = [1]
        @test used_by_variable(pref)
        delete!(m.param_to_vars, index(pref))
    end
    # is_used
    @testset "is_used" begin
        @test !is_used(pref)
        m.param_to_constrs[index(pref)] = [1]
        @test is_used(pref)
        delete!(m.param_to_constrs, index(pref))
        m.param_to_meas[index(pref)] = [1]
        @test is_used(pref)
        delete!(m.param_to_meas, index(pref))
        m.param_to_vars[index(pref)] = [1]
        @test is_used(pref)
        delete!(m.param_to_vars, index(pref))
    end
end


# Test parameter set methods
@testset "Infinite Set" begin
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    # _parameter_set
    @testset "_parameter_set" begin
        @test InfiniteOpt._parameter_set(pref) == IntervalSet(0, 1)
    end
    # _update_parameter_set
    @testset "_update_parameter_set " begin
        @test isa(InfiniteOpt._update_parameter_set(pref,
                                                    IntervalSet(1, 2)), Nothing)
        @test InfiniteOpt._parameter_set(pref) == IntervalSet(1, 2)
    end
    # infinite_set
    @testset "infinite_set" begin
        @test infinite_set(pref) == IntervalSet(1, 2)
    end
    # set_infinite_set
    @testset "set_infinite_set" begin
        @test isa(set_infinite_set(pref, IntervalSet(1, 3)), Nothing)
        @test infinite_set(pref) == IntervalSet(1, 3)
    end
end

# Test parameter support methods
@testset "Supports" begin
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    # _parameter_supports
    @testset "_parameter_supports" begin
        @test InfiniteOpt._parameter_supports(pref) == Number[]
    end
    # _update_parameter_supports
    @testset "_update_parameter_supports " begin
        @test isa(InfiniteOpt._update_parameter_supports(pref, [1]), Nothing)
        @test InfiniteOpt._parameter_supports(pref) == [1]
    end
    # num_supports
    @testset "num_supports" begin
        @test num_supports(pref) == 1
    end
    # has_supports
    @testset "has_supports" begin
        @test has_supports(pref)
        InfiniteOpt._update_parameter_supports(pref, Number[])
        @test !has_supports(pref)
    end
    # supports
    @testset "supports" begin
        @test_throws ErrorException supports(pref)
        InfiniteOpt._update_parameter_supports(pref, [1])
        @test supports(pref) == [1]
    end
    # set_supports
    @testset "set_supports" begin
        @test isa(set_supports(pref, [0, 1]), Nothing)
        @test supports(pref) == [0, 1]
        @test_throws ErrorException set_supports(pref, [2, 3])
        warn = "Support points are not unique, eliminating redundant points."
        @test_logs (:warn, warn) set_supports(pref, [1, 1])
    end
    # add_supports
    @testset "add_supports" begin
        @test isa(add_supports(pref, 0.5), Nothing)
        @test supports(pref) == [0.5, 1]
        @test isa(add_supports(pref, [0, 0.25, 1]), Nothing)
        @test supports(pref) == [0, 0.25, 0.5, 1]
    end
    # delete_supports
    @testset "delete_supports" begin
        @test isa(delete_supports(pref), Nothing)
        @test_throws ErrorException supports(pref)
    end
end

# Test lower bound functions
@testset "Lower Bound" begin
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    # JuMP.has_lower_bound
    @testset "JuMP.has_lower_bound" begin
        @test has_lower_bound(pref)
        set_infinite_set(pref, DistributionSet(Normal()))
        @test has_lower_bound(pref)
        set_infinite_set(pref, DistributionSet(Multinomial(3, 2)))
        @test !has_lower_bound(pref)
        struct test <: AbstractInfiniteSet end
        set_infinite_set(pref, test())
        @test_throws ErrorException has_lower_bound(pref)
    end
    # JuMP.lower_bound
    @testset "JuMP.lower_bound" begin
        set_infinite_set(pref, DistributionSet(Multinomial(3, 2)))
        @test_throws ErrorException lower_bound(pref)
        set_infinite_set(pref, IntervalSet(0, 1))
        @test lower_bound(pref) == 0
        set_infinite_set(pref, DistributionSet(Normal()))
        @test lower_bound(pref) == -Inf
    end
    # JuMP.set_lower_bound
    @testset "JuMP.set_lower_bound" begin
        @test_throws ErrorException set_lower_bound(pref, 2)
        struct test <: AbstractInfiniteSet end
        set_infinite_set(pref, test())
        @test_throws ErrorException set_lower_bound(pref, 2)
        set_infinite_set(pref, IntervalSet(0, 1))
        @test isa(set_lower_bound(pref, -1), Nothing)
        @test lower_bound(pref) == -1
    end
end

# Test upper bound functions
@testset "Upper Bound" begin
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    # JuMP.has_upper_bound
    @testset "JuMP.has_upper_bound" begin
        @test has_upper_bound(pref)
        set_infinite_set(pref, DistributionSet(Normal()))
        @test has_upper_bound(pref)
        set_infinite_set(pref, DistributionSet(Multinomial(3, 2)))
        @test !has_upper_bound(pref)
        struct test <: AbstractInfiniteSet end
        set_infinite_set(pref, test())
        @test_throws ErrorException has_upper_bound(pref)
    end
    # JuMP.upper_bound
    @testset "JuMP.upper_bound" begin
        set_infinite_set(pref, DistributionSet(Multinomial(3, 2)))
        @test_throws ErrorException upper_bound(pref)
        set_infinite_set(pref, IntervalSet(0, 1))
        @test upper_bound(pref) == 1
        set_infinite_set(pref, DistributionSet(Normal()))
        @test upper_bound(pref) == Inf
    end
    # JuMP.set_upper_bound
    @testset "JuMP.set_upper_bound" begin
        @test_throws ErrorException set_upper_bound(pref, 2)
        struct test <: AbstractInfiniteSet end
        set_infinite_set(pref, test())
        @test_throws ErrorException set_upper_bound(pref, 2)
        set_infinite_set(pref, IntervalSet(0, 1))
        @test isa(set_upper_bound(pref, 2), Nothing)
        @test upper_bound(pref) == 2
    end
end

# Test everything else
@testset "Other Queries" begin
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    # JuMP.is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, pref)
        pref2 = ParameterRef(InfiniteModel(), 1)
        @test !is_valid(m, pref2)
        pref3 = ParameterRef(m, 2)
        @test !is_valid(m, pref3)
    end
    # num_parameters
    @testset "num_parameters" begin
        @test num_parameters(m) == 1
        delete!(m.params, JuMP.index(pref))
        @test num_parameters(m) == 0
        m.params[JuMP.index(pref)] = param
    end
    # all_parameters
    @testset "all_parameters" begin
        @test length(all_parameters(m)) == 1
        @test all_parameters(m)[1] == pref
    end
    # is_independent
    @testset "is_indepentent" begin
        @test !is_independent(pref)
    end
    # group_id
    @testset "group_id" begin
        @test group_id(pref) == 0
    end

    # _groups
    @testset "_groups" begin

    end
end
