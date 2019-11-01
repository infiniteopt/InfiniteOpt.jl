# Test Infinite Sets
@testset "Infinite Sets" begin
    # Abstract types
    @test AbstractInfiniteSet isa DataType
    # IntervalSet
    @test IntervalSet <: AbstractInfiniteSet
    @test IntervalSet(0, 1) == IntervalSet(0.0, 1.0)
    @test_throws ErrorException IntervalSet(1, 0)
    # DistributionSet
    @test DistributionSet <: AbstractInfiniteSet
    @test DistributionSet(Normal()) isa DistributionSet{<:Normal}
end

# Test core parameter/variable data structures
sample_info = VariableInfo(zeros(Bool, 10)...)
@testset "Core Parameter/Variables" begin
    # Parameters
    @test InfOptParameter <: AbstractVariable
    @test InfOptParameter(IntervalSet(0, 1), Int[], true).set isa IntervalSet
    p1 = InfOptParameter(IntervalSet(0, 1), Int[], true)
    p2 = InfOptParameter(IntervalSet(0, 1), Int[], true)
    @test p1 == p2
    # Abstract variable
    @test InfOptVariable <: AbstractVariable
    @test AbstractReducedInfo isa DataType
    # Infinite variable
    @test InfiniteVariable <: InfOptVariable
    @test InfiniteVariable(sample_info, (1, 2)).parameter_refs isa Tuple
    # Point variable
    @test PointVariable <: InfOptVariable
    # Hold variable
    @test ParameterBounds isa DataType
    @test ParameterBounds().intervals == Dict{ParameterRef, IntervalSet}()
    @test HoldVariable <: InfOptVariable
    @test HoldVariable(sample_info,
                       ParameterBounds()).info isa VariableInfo
    # Reduced variable info
    @test ReducedInfiniteInfo <: AbstractReducedInfo
end

# Test the InfiniteModel datatype
@testset "InfiniteModel" begin
    @test InfiniteModel <: JuMP.AbstractModel
    @test InfiniteModel().next_var_index == 0
    mockoptimizer = with_optimizer(MOIU.MockOptimizer,
                                   MOIU.Model{Float64}(),
                                   eval_objective_value=false)
    @test InfiniteModel(mockoptimizer).optimizer_factory.constructor == MOIU.MockOptimizer
    m = InfiniteModel();
    @test isa(Base.broadcastable(m), Base.RefValue{InfiniteModel})
    @test length(JuMP.object_dictionary(m)) == 0
    @test InfiniteModel(caching_mode = MOIU.MANUAL) isa JuMP.AbstractModel
    @test InfiniteModel(mockoptimizer,
                        caching_mode = MOIU.MANUAL) isa JuMP.AbstractModel
    m = InfiniteModel(seed = true)
    #@test reinterpret(Int32, Random.GLOBAL_RNG.seed)[1] == 0 # this test does not work on travis nightly julia
    rand_num1 = rand()
    m = InfiniteModel(seed = true)
    rand_num2 = rand()
    @test rand_num1 == rand_num2
end

# Test reference variable datatypes
@testset "References" begin
    m = InfiniteModel();
    # Abstract types
    @test GeneralVariableRef <: JuMP.AbstractVariableRef
    @test MeasureFiniteVariableRef <: GeneralVariableRef
    @test FiniteVariableRef <: MeasureFiniteVariableRef
    # Hold variable refs
    @test HoldVariableRef <: FiniteVariableRef
    @test HoldVariableRef(m, 1).index == 1
    # Point variable refs
    @test PointVariableRef <: FiniteVariableRef
    @test PointVariableRef(m, 1).index == 1
    # Infinite variable refs
    @test InfiniteVariableRef <: GeneralVariableRef
    @test InfiniteVariableRef(m, 1).index == 1
    ivref = InfiniteVariableRef(m, 1)
    # Reduced infinite variable refs
    @test ReducedInfiniteVariableRef <: GeneralVariableRef
    @test ReducedInfiniteVariableRef(m, 1).index == 1
    # Parameter refs
    @test ParameterRef <: GeneralVariableRef
    @test ParameterRef(m, 1).index == 1
    pref = ParameterRef(m, 1)
    @test copy(pref, m).index == 1
    # Measure refs
    @test MeasureRef <: MeasureFiniteVariableRef
    @test MeasureRef(m, 1).index == 1
end

# Test the constraint datatypes
@testset "Unions" begin
    @test InfOptVariableRef <: GeneralVariableRef
    @test InfiniteExpr <: AbstractJuMPScalar
    @test ParameterExpr <: AbstractJuMPScalar
    @test MeasureExpr <: AbstractJuMPScalar
end

# Test the constraint datatypes
@testset "Constraints" begin
    m = InfiniteModel()
    # Bounded constraints
    @test BoundedScalarConstraint <: JuMP.AbstractConstraint
    pref = ParameterRef(m, 1)
    bounds = ParameterBounds(Dict(pref => IntervalSet(0, 1)))
    @test BoundedScalarConstraint(zero(AffExpr), MOI.Integer(),
                bounds, bounds).bounds.intervals[pref].lower_bound == 0.0
    # Abstract cosntraint refs
    @test GeneralConstraintRef isa DataType
    # Infinite constraint refs
    @test InfiniteConstraintRef <: GeneralConstraintRef
    @test InfiniteConstraintRef(m, 1, JuMP.ScalarShape()).index == 1
    # Finite constraint refs
    @test FiniteConstraintRef <: GeneralConstraintRef
    @test FiniteConstraintRef(m, 1, JuMP.ScalarShape()).index == 1
    # Measure constraint refs
    @test MeasureConstraintRef <: GeneralConstraintRef
    @test MeasureConstraintRef(m, 1, JuMP.ScalarShape()).index == 1
end

# Test ParameterBounds
@testset "Parameter Bounds" begin
    # setup
    m = InfiniteModel()
    m.params[1] = InfOptParameter(IntervalSet(0, 10), Number[], false)
    m.params[2] = InfOptParameter(IntervalSet(0, 10), Number[], false)
    m.params[3] = InfOptParameter(IntervalSet(0, 10), Number[], false)
    par = ParameterRef(m, 1)
    pars = [ParameterRef(m, 2), ParameterRef(m, 3)]
    # test _expand_parameter_dict(Dict{ParameterRef,IntervalSet}))
    @testset "_expand_parameter_dict (acceptable Form)" begin
        d = Dict(par => IntervalSet(0, 1))
        @test InfiniteOpt._expand_parameter_dict(d) == d
    end
    # test _expand_parameter_dict(Dict{Any,IntervalSet}))
    @testset "_expand_parameter_dict (Array Form)" begin
        d = Dict(pars => IntervalSet(0, 1), par => IntervalSet(0, 1))
        @test isa(InfiniteOpt._expand_parameter_dict(d),
                  Dict{ParameterRef, IntervalSet})
    end
    # test _expand_parameter_dict(Dict))
    @testset "_expand_parameter_dict (Fallback)" begin
        d = Dict(pars => 1, par => 2)
        @test_throws ErrorException InfiniteOpt._expand_parameter_dict(d)
    end
    # test expansion definition
    @testset "ParameterBounds Expansion" begin
        d = Dict(par => IntervalSet(0, 1))
        @test ParameterBounds(d).intervals == d
        d = Dict(pars => IntervalSet(0, 1), par => IntervalSet(0, 1))
        @test isa(ParameterBounds(d).intervals, Dict{ParameterRef, IntervalSet})
    end
end
