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
end

# TODO test owner_model, index, set_name, and similar first
@testset "Basic Reference Queries" begin

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
        param = InfOptParameter(set, Number[], false)
        expected = ParameterRef(m, 1)
        @test add_parameter(m, param) == expected
    end
end
