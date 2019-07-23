# Test constraint deletion
@testset "JuMP.delete (Constraints)" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0.5), pt)
    @global_variable(m, x)
    data = DiscreteMeasureData(par, [1], [1])
    meas = measure(inf + par - x, data)
    rv = ReducedInfiniteVariableRef(m, 42)
    @constraint(m, cref, par - inf + pt + 2x - rv + meas <= 1)
    # test normal deletion
    @test isa(delete(m, cref), Nothing)
    @test !is_valid(m, cref)
    @test !haskey(m.param_to_constrs, JuMP.index(par))
    @test !haskey(m.var_to_constrs, JuMP.index(inf))
    @test !haskey(m.var_to_constrs, JuMP.index(pt))
    @test !haskey(m.var_to_constrs, JuMP.index(x))
    @test !haskey(m.reduced_to_constrs, JuMP.index(rv))
    @test !haskey(m.meas_to_constrs, JuMP.index(meas))
    @test !haskey(m.constrs, JuMP.index(cref))
    @test !haskey(m.constr_to_name, JuMP.index(cref))
    @test !haskey(m.constr_in_var_info, JuMP.index(cref))
    # test assertion error
    @test_throws AssertionError delete(m, cref)
end

# Test parameter deletion
@testset "Parameters" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par <= 1)
    @infinite_variable(m, inf(par))
    @point_variable(m, inf(0.5), pt)
    @global_variable(m, x)
    data = DiscreteMeasureData(par, [1], [1])
    meas = measure(inf + par - x, data)
    rv = ReducedInfiniteVariableRef(m, 42)
end
