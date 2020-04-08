# Test methods for reduced infinite variables
@testset "Reduced Variable" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 1 <= par1 <= 2)
    @infinite_parameter(m, 1 <= par2 <= 2)
    @infinite_variable(m, 0 <= inf1(par1, par2) <= 1, Int)
    @infinite_variable(m, inf2(par1, par2) == 1, Bin, start = 0)
    index = m.next_var_index + 1
    rvref1 = ReducedInfiniteVariableRef(m, index)
    m.reduced_info[index] = ReducedInfiniteInfo(inf1, Dict(1 => 1))
    rvref2 = ReducedInfiniteVariableRef(m, index + 1)
    m.reduced_info[index + 1] = ReducedInfiniteInfo(inf2, Dict(2 => 1))
    tm = optimizer_model(m)
    rvref3 = ReducedInfiniteVariableRef(m, -1)
    tm.ext[:TransData].reduced_info[-1] = ReducedInfiniteInfo(inf2, Dict(2 => 1))
    # test reduction_info
    @testset "reduction_info" begin
        @test_throws ErrorException reduction_info(rvref1, Val(:bad))
    end
    # test _reduced_info
    @testset "_reduced_info" begin
        # test reduced in the infinite model
        @test InfiniteOpt._reduced_info(rvref1).infinite_variable_ref == inf1
        @test InfiniteOpt._reduced_info(rvref2).infinite_variable_ref == inf2
        @test InfiniteOpt._reduced_info(rvref1).eval_supports == Dict(1 => 1)
        @test InfiniteOpt._reduced_info(rvref2).eval_supports == Dict(2 => 1)
        # test using reduction info
        @test InfiniteOpt._reduced_info(rvref3).infinite_variable_ref == inf2
        @test InfiniteOpt._reduced_info(rvref3).eval_supports == Dict(2 => 1)
    end
    # test infinite_variable_ref
    @testset "infinite_variable_ref" begin
        @test infinite_variable_ref(rvref1) == inf1
        @test infinite_variable_ref(rvref2) == inf2
    end
    # test eval_supports
    @testset "eval_supports" begin
        @test eval_supports(rvref1) == Dict(1 => 1)
        @test eval_supports(rvref2) == Dict(2 => 1)
    end
    # test name
    @testset "JuMP.name" begin
        @test name(rvref1) == "inf1(1, par2)"
        @test name(rvref2) == "inf2(par1, 1)"
    end
    # raw_parameter_refs
    @testset "raw_parameter_refs" begin
        # test the references
        @test raw_parameter_refs(rvref1) == VectorTuple(par2)
        @test raw_parameter_refs(rvref2) == VectorTuple(par1)
    end
    # parameter_refs (reduced infinite variable)
    @testset "parameter_refs" begin
        # test the references
        @test parameter_refs(rvref1) == (par2, )
        @test parameter_refs(rvref2) == (par1, )
    end
    # parameter_list (reduced infinite variable)
    @testset "parameter_list" begin
        # test the references
        @test parameter_list(rvref1) == [par2]
        @test parameter_list(rvref2) == [par1]
    end
    # has_lower_bound
    @testset "JuMP.has_lower_bound" begin
        @test has_lower_bound(rvref1)
        @test !has_lower_bound(rvref2)
    end
    # lower_bound
    @testset "JuMP.lower_bound" begin
        @test lower_bound(rvref1) == 0
        @test_throws ErrorException lower_bound(rvref2)
    end
    # _lower_bound_index
    @testset "JuMP._lower_bound_index" begin
        @test JuMP._lower_bound_index(rvref1) == JuMP._lower_bound_index(inf1)
        @test_throws ErrorException JuMP._lower_bound_index(rvref2)
    end
    # LowerBoundRef
    @testset "JuMP.LowerBoundRef" begin
        @test LowerBoundRef(rvref1) == LowerBoundRef(inf1)
        @test_throws ErrorException LowerBoundRef(rvref2)
    end
    # has_upper_bound
    @testset "JuMP.has_upper_bound" begin
        @test has_upper_bound(rvref1)
        @test !has_upper_bound(rvref2)
    end
    # upper_bound
    @testset "JuMP.upper_bound" begin
        @test upper_bound(rvref1) == 1
        @test_throws ErrorException upper_bound(rvref2)
    end
    # _upper_bound_index
    @testset "JuMP._upper_bound_index" begin
        @test JuMP._upper_bound_index(rvref1) == JuMP._upper_bound_index(inf1)
        @test_throws ErrorException JuMP._upper_bound_index(rvref2)
    end
    # UpperBoundRef
    @testset "JuMP.UpperBoundRef" begin
        @test UpperBoundRef(rvref1) == UpperBoundRef(inf1)
        @test_throws ErrorException UpperBoundRef(rvref2)
    end
    # is_fixed
    @testset "JuMP.is_fixed" begin
        @test is_fixed(rvref2)
        @test !is_fixed(rvref1)
    end
    # fix_value
    @testset "JuMP.fix_value" begin
        @test fix_value(rvref2) == 1
        @test_throws ErrorException fix_value(rvref1)
    end
    # _fix_index
    @testset "JuMP._fix_index" begin
        @test JuMP._fix_index(rvref2) == JuMP._fix_index(inf2)
        @test_throws ErrorException JuMP._fix_index(rvref1)
    end
    # FixRef
    @testset "JuMP.FixRef" begin
        @test FixRef(rvref2) == FixRef(inf2)
        @test_throws ErrorException FixRef(rvref1)
    end
    # start_value
    @testset "JuMP.start_value" begin
        @test start_value(rvref2) == 0
    end
    # is_binary
    @testset "JuMP.is_binary" begin
        @test is_binary(rvref2)
        @test !is_binary(rvref1)
    end
    # _binary_index
    @testset "JuMP._binary_index" begin
        @test JuMP._binary_index(rvref2) == JuMP._binary_index(inf2)
        @test_throws ErrorException JuMP._binary_index(rvref1)
    end
    # BinaryRef
    @testset "JuMP.BinaryRef" begin
        @test BinaryRef(rvref2) == BinaryRef(inf2)
        @test_throws ErrorException BinaryRef(rvref1)
    end
    # is_integer
    @testset "JuMP.is_integer" begin
        @test is_integer(rvref1)
        @test !is_integer(rvref2)
    end
    # _integer_index
    @testset "JuMP._integer_index" begin
        @test JuMP._integer_index(rvref1) == JuMP._integer_index(inf1)
        @test_throws ErrorException JuMP._integer_index(rvref2)
    end
    # IntegerRef
    @testset "JuMP.IntegerRef" begin
        @test IntegerRef(rvref1) == IntegerRef(inf1)
        @test_throws ErrorException IntegerRef(rvref2)
    end
    # used_by_constraint
    @testset "used_by_constraint" begin
        # test not used
        @test !used_by_constraint(rvref1)
        # prepare use case
        m.reduced_to_constrs[JuMP.index(rvref1)] = [1]
        # test used
        @test used_by_constraint(rvref1)
    end
    # used_by_measure
    @testset "used_by_measure" begin
        # test not used
        @test !used_by_measure(rvref1)
        # prepare use case
        m.reduced_to_meas[JuMP.index(rvref1)] = [1]
        # test used
        @test used_by_measure(rvref1)
    end
    # test is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, rvref1)
        @test is_valid(m, rvref2)
        @test !is_valid(m, ReducedInfiniteVariableRef(m, 100))
    end
end
