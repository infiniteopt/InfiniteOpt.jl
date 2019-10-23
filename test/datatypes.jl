# Test Infinite Sets
@testset "Infinite Sets" begin
    # Abstract types
    @test AbstractInfiniteSet isa DataType
    # IntervalSet
    @test IntervalSet <: AbstractInfiniteSet
    @test IntervalSet(0, 1) == IntervalSet(0.0, 1.0)
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
    # Global variable
    @test GlobalVariable <: InfOptVariable
    @test GlobalVariable(sample_info).info isa VariableInfo
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
    @test reinterpret(Int32, Random.GLOBAL_RNG.seed)[1] == 0
end

# Test reference variable datatypes
@testset "References" begin
    m = InfiniteModel();
    # Abstract types
    @test GeneralVariableRef <: JuMP.AbstractVariableRef
    @test MeasureFiniteVariableRef <: GeneralVariableRef
    @test FiniteVariableRef <: MeasureFiniteVariableRef
    # Global variable refs
    @test GlobalVariableRef <: FiniteVariableRef
    @test GlobalVariableRef(m, 1).index == 1
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
    m = InfiniteModel();
    # Bounded constraints
    @test BoundedScalarConstraint <: JuMP.AbstractConstraint
    pref = ParameterRef(m, 1);
    dict = Dict(pref => IntervalSet(0, 1));
    @test BoundedScalarConstraint(zero(AffExpr), MOI.Integer(),
                                  dict).bounds[pref].lower_bound == 0.0
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
