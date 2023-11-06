# Test indexes
@testset "Indexes" begin
    # Abstract types
    @test AbstractInfOptIndex isa DataType
    @test ObjectIndex isa DataType
    @test ObjectIndex <: AbstractInfOptIndex
    # IndependentParameterIndex
    @test IndependentParameterIndex <: ObjectIndex
    @test IndependentParameterIndex(1).value == 1
    # DependentParametersIndex
    @test DependentParametersIndex <: ObjectIndex
    @test DependentParametersIndex(1).value == 1
    # DependentParameterIndex
    @test DependentParameterIndex <: AbstractInfOptIndex
    idx = DependentParametersIndex(1)
    @test DependentParameterIndex(idx, 1).object_index == idx
    @test DependentParameterIndex(idx, 1).param_index == 1
    # InfiniteParameterIndex
    @test InfiniteParameterIndex <: AbstractInfOptIndex
    # FiniteParameterIndex
    @test FiniteParameterIndex <: ObjectIndex
    @test FiniteParameterIndex(1).value == 1
    # ParameterFunctionIndex
    @test ParameterFunctionIndex <: ObjectIndex
    @test ParameterFunctionIndex(1).value == 1
    # InfiniteVariableIndex
    @test InfiniteVariableIndex <: ObjectIndex
    @test InfiniteVariableIndex(1).value == 1
    # SemiInfiniteVariableIndex
    @test SemiInfiniteVariableIndex <: ObjectIndex
    @test SemiInfiniteVariableIndex(1).value == 1
    # PointVariableIndex
    @test PointVariableIndex <: ObjectIndex
    @test PointVariableIndex(1).value == 1
    # FiniteVariableIndex
    @test FiniteVariableIndex <: ObjectIndex
    @test FiniteVariableIndex(1).value == 1
    # DerivativeIndex
    @test DerivativeIndex <: ObjectIndex
    @test DerivativeIndex(1).value == 1
    # MeasureIndex
    @test MeasureIndex <: ObjectIndex
    @test MeasureIndex(1).value == 1
    # InfOptConstraintIndex
    @test InfOptConstraintIndex <: ObjectIndex
    @test InfOptConstraintIndex(1).value == 1
    # test CleverDict extensions
    @testset "CleverDict Extensions" begin
        # index_to_key
        @test MOIUC.index_to_key(IndependentParameterIndex, Int64(42)) == IndependentParameterIndex(42)
        @test MOIUC.index_to_key(DependentParametersIndex, Int64(42)) == DependentParametersIndex(42)
        @test MOIUC.index_to_key(FiniteParameterIndex, Int64(42)) == FiniteParameterIndex(42)
        @test MOIUC.index_to_key(ParameterFunctionIndex, Int64(42)) == ParameterFunctionIndex(42)
        @test MOIUC.index_to_key(InfiniteVariableIndex, Int64(42)) == InfiniteVariableIndex(42)
        @test MOIUC.index_to_key(SemiInfiniteVariableIndex, Int64(42)) == SemiInfiniteVariableIndex(42)
        @test MOIUC.index_to_key(PointVariableIndex, Int64(42)) == PointVariableIndex(42)
        @test MOIUC.index_to_key(FiniteVariableIndex, Int64(42)) == FiniteVariableIndex(42)
        @test MOIUC.index_to_key(DerivativeIndex, Int64(42)) == DerivativeIndex(42)
        @test MOIUC.index_to_key(MeasureIndex, Int64(42)) == MeasureIndex(42)
        @test MOIUC.index_to_key(InfOptConstraintIndex, Int64(42)) == InfOptConstraintIndex(42)
        # index_to_key
        @test MOIUC.key_to_index(MeasureIndex(42)) == 42
    end
    # test Base.length
    @testset "Base.length" begin
        @test length(MeasureIndex(1)) == 1
    end
    # test Base.broadcastable
    @testset "Base.broadcastable" begin
        @test Base.broadcastable(MeasureIndex(1)) isa Ref
    end
    # test Base.iterate
    @testset "Base.iterate" begin
        @test iterate(MeasureIndex(1)) == (MeasureIndex(1), true)
        @test iterate(MeasureIndex(1), true) isa Nothing
    end
    # test Base.hash
    @testset "Base.hash" begin
        @test Base.hash(MeasureIndex(1), UInt(1)) isa UInt
    end
end

# Test Infinite Domains
@testset "Infinite Domains" begin
    # Abstract types
    @test AbstractInfiniteDomain isa DataType
    @test InfiniteScalarDomain <: AbstractInfiniteDomain
    @test InfiniteArrayDomain <: AbstractInfiniteDomain
    # IntervalDomain
    @test IntervalDomain <: InfiniteScalarDomain
    @test IntervalDomain(0, 1) == IntervalDomain(0.0, 1.0)
    @test_throws ErrorException IntervalDomain(1, 0)
    # UniDistributionDomain
    @test UniDistributionDomain <: InfiniteScalarDomain
    @test UniDistributionDomain(Normal()) isa UniDistributionDomain{<:Normal}
    # NonUnivariateDistribution
    @test InfiniteOpt.NonUnivariateDistribution <: Distribution
    # MultiDistributionDomain
    dist = Dirichlet([0.5, 0.5])
    @test MultiDistributionDomain <: InfiniteArrayDomain
    @test MultiDistributionDomain(dist) isa MultiDistributionDomain{<:Dirichlet}
    @test MultiDistributionDomain(dist) == MultiDistributionDomain(dist)
    @test MultiDistributionDomain(dist) != MultiDistributionDomain(MvNormal([1, 1], [1 0; 0 1]))
    # CollectionDomain
    @test CollectionDomain <: InfiniteArrayDomain
    @test CollectionDomain([IntervalDomain(0, 1), IntervalDomain(2, 3)]) isa CollectionDomain{IntervalDomain}
end

# Test Generative Support Info 
@testset "Generative Support Info" begin 
    # Abstract types
    @test AbstractGenerativeInfo isa DataType
    # NoGenerativeSupports
    @test NoGenerativeSupports <: AbstractGenerativeInfo
    @test NoGenerativeSupports() isa NoGenerativeSupports
    @test NoGenerativeSupports() == NoGenerativeSupports()
    # UniformGenerativeInfo
    @test UniformGenerativeInfo <: AbstractGenerativeInfo
    @test UniformGenerativeInfo([0.0, 1.0], InternalLabel) isa UniformGenerativeInfo
    @test UniformGenerativeInfo([-2, 0, 1], InternalLabel, -2, 1).support_basis == [0, 2.0/3.0, 1]
    @test UniformGenerativeInfo([0.0, 1.0], InternalLabel) == UniformGenerativeInfo([2, 6], InternalLabel, 2, 6)
    @test UniformGenerativeInfo([0.0, 1.0], InternalLabel) != UniformGenerativeInfo([2, 5], InternalLabel, 2, 6)
    @test_throws ErrorException UniformGenerativeInfo([2, 6], InternalLabel)
    @test_throws ErrorException UniformGenerativeInfo(Float64[], InternalLabel)
end

# Test Derivative Evaluation Methods 
@testset "Derivative Methods" begin
    # abstract types
    @test AbstractDerivativeMethod isa DataType
    @test GenerativeDerivativeMethod <: AbstractDerivativeMethod
    @test NonGenerativeDerivativeMethod <: AbstractDerivativeMethod
    # test orthogonal collocation (defined in derivative_evaluations.jl)
    @test OrthogonalCollocation <: GenerativeDerivativeMethod
    @test OrthogonalCollocation{GaussLobatto}(3, GaussLobatto()) isa OrthogonalCollocation
    @test OrthogonalCollocation(3, GaussLobatto()) == OrthogonalCollocation{GaussLobatto}(3, GaussLobatto())
    @test OrthogonalCollocation(3) == OrthogonalCollocation(3, GaussLobatto())
    @test_throws ErrorException OrthogonalCollocation(1)
    # test finite difference
    @test FDTechnique isa DataType
    @test Forward <: FDTechnique
    @test Central <: FDTechnique
    @test Backward <: FDTechnique
    @test FiniteDifference <: NonGenerativeDerivativeMethod
    @test FiniteDifference(Forward()) isa FiniteDifference
    @test FiniteDifference() == FiniteDifference(Backward())
    @test FiniteDifference(Forward(), false) isa FiniteDifference
end

# Test parameter datatypes
@testset "Parameters" begin
    # abstract types
    @test InfOptParameter <: AbstractVariable
    @test ScalarParameter <: InfOptParameter
    @test AbstractDataObject isa DataType
    # test IndependentParameter
    @test IndependentParameter <: ScalarParameter
    dict = SortedDict{Float64, Set{DataType}}(2 => Set([All]))
    method = FiniteDifference()
    info = NoGenerativeSupports()
    @test IndependentParameter(IntervalDomain(0, 1), dict, 6, method, info).domain isa IntervalDomain
    # test FiniteParameter
    @test FiniteParameter <: ScalarParameter
    @test FiniteParameter(1) == FiniteParameter(1)
    # test DependentParameters
    @test DependentParameters <: InfOptParameter
    @test DependentParameters(CollectionDomain([IntervalDomain(0, 1)]),
                              Dict(zeros(1) => Set([All])), 6, [method]).domain isa CollectionDomain
    # test ScalarParameterData
    @test ScalarParameterData <: AbstractDataObject
    @test ScalarParameterData(FiniteParameter(42), 1, 1, "bob").name == "bob"
    # test MultiParameterData
    @test MultiParameterData <: AbstractDataObject
    params = DependentParameters(CollectionDomain([IntervalDomain(0, 1)]),
                                 Dict(zeros(1) => Set([All])), 6, [method])
    @test MultiParameterData(params, 1, 1:1, ["par[1]"]) isa MultiParameterData
end

# Test the InfiniteModel datatype
@testset "InfiniteModel" begin
    # test basic
    @test InfiniteModel <: JuMP.AbstractModel
    @test InfiniteModel().ext isa Dict
    # prepare optimizer constructor
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    mockattributes = MOI.OptimizerWithAttributes(mockoptimizer, MOI.Silent() => true)
    # test optimizer constructors
    @test InfiniteModel(mockoptimizer).optimizer_constructor == mockoptimizer
    @test InfiniteModel(mockattributes).optimizer_constructor == mockoptimizer
    m = InfiniteModel();
    @test isa(Base.broadcastable(m), Base.RefValue{InfiniteModel})
    @test length(JuMP.object_dictionary(m)) == 0
    @test InfiniteModel() isa JuMP.AbstractModel
    @test InfiniteModel(mockoptimizer, add_bridges = false) isa JuMP.AbstractModel
    # test accessors
    @test InfiniteOpt._last_param_num(m) == 0
    @test InfiniteOpt._param_object_indices(m) isa Vector{Union{IndependentParameterIndex, DependentParametersIndex}}
    # test other methods 
    @test empty!(InfiniteModel(mockoptimizer)).optimizer_constructor == mockoptimizer
    @test variable_ref_type(InfiniteModel) == GeneralVariableRef
    @test variable_ref_type(InfiniteModel()) == GeneralVariableRef
end

# Test reference variable datatypes
@testset "References" begin
    m = InfiniteModel();
    # Abstract types
    @test DispatchVariableRef <: AbstractVariableRef
    @test FiniteRef <: DispatchVariableRef
    # test GeneralVariableRef
    @test GeneralVariableRef <: AbstractVariableRef
    @test GeneralVariableRef(m, 1, MeasureIndex) isa GeneralVariableRef
    # test IndependentParameterRef
    @test IndependentParameterRef <: DispatchVariableRef
    idx = IndependentParameterIndex(1)
    @test IndependentParameterRef(m, idx) isa IndependentParameterRef
    # test DependentParameterRef
    @test DependentParameterRef <: DispatchVariableRef
    idx = DependentParameterIndex(DependentParametersIndex(1), 1)
    @test DependentParameterRef(m, idx) isa DependentParameterRef
    # test FiniteParameterRef
    @test FiniteParameterRef <: FiniteRef
    idx = FiniteParameterIndex(1)
    @test FiniteParameterRef(m, idx) isa FiniteParameterRef
    # InfiniteVariableRef
    @test InfiniteVariableRef <: DispatchVariableRef
    idx = InfiniteVariableIndex(1)
    @test InfiniteVariableRef(m, idx) isa InfiniteVariableRef
    # SemiInfiniteVariableRef
    @test SemiInfiniteVariableRef <: DispatchVariableRef
    idx = SemiInfiniteVariableIndex(1)
    @test SemiInfiniteVariableRef(m, idx) isa SemiInfiniteVariableRef
    # ParameterFunctionRef
    @test ParameterFunctionRef <: DispatchVariableRef
    idx = ParameterFunctionIndex(1)
    @test ParameterFunctionRef(m, idx) isa ParameterFunctionRef
    # PointVariableRef
    @test PointVariableRef <: FiniteRef
    idx = PointVariableIndex(1)
    @test PointVariableRef(m, idx) isa PointVariableRef
    # FiniteVariableRef
    @test FiniteVariableRef <: FiniteRef
    idx = FiniteVariableIndex(1)
    @test FiniteVariableRef(m, idx) isa FiniteVariableRef
    # DerivativeRef
    @test DerivativeRef <: DispatchVariableRef
    idx = DerivativeIndex(1)
    @test DerivativeRef(m, idx) isa DerivativeRef
    # MeasureRef
    @test MeasureRef <: DispatchVariableRef
    idx = MeasureIndex(1)
    @test MeasureRef(m, idx) isa MeasureRef
    # InfOptConstraintRef
    @test InfOptConstraintRef isa DataType
    idx = InfOptConstraintIndex(1)
    @test InfOptConstraintRef(m, idx) isa InfOptConstraintRef
end

# Test DomainRestrictions
@testset "Domain Restrictions" begin
    # prepare the data
    m = InfiniteModel()
    idx = DependentParametersIndex(1)
    idx1 = DependentParameterIndex(idx, 1)
    idx2 = DependentParameterIndex(idx, 2)
    idx3 = IndependentParameterIndex(1)
    par1 = GeneralVariableRef(m, 1, DependentParameterIndex, 1)
    par2 = GeneralVariableRef(m, 1, DependentParameterIndex, 2)
    par3 = GeneralVariableRef(m, 1, IndependentParameterIndex, -1)
    pars = [par1, par2]
    # test datatype
    @testset "DataType" begin
        @test DomainRestrictions isa UnionAll
        d = Dict(par3 => IntervalDomain(0, 1))
        @test DomainRestrictions(d) isa DomainRestrictions{GeneralVariableRef}
        @test DomainRestrictions() isa DomainRestrictions{GeneralVariableRef}
    end
    # test _expand_parameter_tuple
    @testset "_expand_parameter_tuple" begin
        # test good to go
        d = (par3 => IntervalDomain(0, 1),)
        @test InfiniteOpt._expand_parameter_tuple(d) == Dict(d...)
        # test with array of parameters
        d = (pars => IntervalDomain(0, 1), par3 => IntervalDomain(0, 1))
        @test InfiniteOpt._expand_parameter_tuple(d) isa Dict
        # test symbolic domains 
        d = (pars => [0, 1], par3 => 0)
        @test InfiniteOpt._expand_parameter_tuple(d) isa Dict
        # test bad vector domain
        d = (pars => [0, 1, 3],)
        @test_throws ErrorException InfiniteOpt._expand_parameter_tuple(d)
        # test fallback 
        d = (pars => :a,)
        @test_throws ErrorException InfiniteOpt._expand_parameter_tuple(d)
    end
    # test expansion definition
    @testset "DomainRestrictions Expansion" begin
        d = (par3 => IntervalDomain(0, 1),)
        @test DomainRestrictions(d).intervals == Dict(d...)
        d = (pars => IntervalDomain(0, 1), par3 => IntervalDomain(0, 1))
        @test DomainRestrictions(d).intervals isa Dict{GeneralVariableRef, IntervalDomain}
        @test DomainRestrictions(par3 => 0).intervals isa Dict
    end
    dr = DomainRestrictions(par3 => [0, 1])
    # test intervals
    @testset "intervals" begin
        @test intervals(dr) == dr.intervals
    end
    # test Base.length
    @testset "Base.length" begin
        @test length(dr) == 1
    end
    # test Base.isempty
    @testset "Base.isempty" begin
        @test !isempty(dr)
        @test isempty(DomainRestrictions())
    end
    # test Base.:(==)
    @testset "Base.:(==)" begin
        @test dr == DomainRestrictions((par3 => IntervalDomain(0, 1),))
        @test dr != DomainRestrictions()
    end
    # test Base.copy
    @testset "Base.copy" begin
        @test copy(dr) == dr
    end
    # test Base.getindex
    @testset "Base.getindex" begin
        @test dr[par3] == IntervalDomain(0, 1)
    end
    # test Base.setindex!
    @testset "Base.setindex!" begin
        @test (dr[par3] = IntervalDomain(0, 2)) == IntervalDomain(0, 2)
    end
    # test Base.haskey
    @testset "Base.haskey" begin
        @test haskey(dr, par3)
        @test !haskey(dr, par2)
    end
    # test Base.keys
    @testset "Base.keys" begin
        @test keys(dr) == keys(intervals(dr))
    end
    # test Base.iterate
    @testset "Base.iterate" begin
        @test [r for r in dr] == [par3 => IntervalDomain(0, 2)]
    end
    # test Base.merge
    @testset "Base.merge" begin
        new_dr = DomainRestrictions()
        @test merge(new_dr, dr) == dr
    end
    # test Base.merge!
    @testset "Base.merge!" begin
        new_dr = DomainRestrictions()
        @test merge!(new_dr, dr) isa DomainRestrictions
        @test new_dr == dr
    end
    # test Base.filter
    @testset "Base.filter" begin
        @test filter(e -> e[1] != par3, dr) == DomainRestrictions()
    end
    # test Base.delete!
    @testset "Base.delete!" begin
        @test delete!(dr, par3) isa DomainRestrictions
        @test length(dr) == 0
    end
end

# Test parameter functions 
@testset "Parameter Functions" begin 
    # useful data
    m = InfiniteModel()
    pref = GeneralVariableRef(m, 1, IndependentParameterIndex, -1)
    vt = IC.VectorTuple(pref)
    # test ParameterFunction
    @test ParameterFunction(sin, vt, [1], [1]) isa ParameterFunction
    # test ParameterFunctionData
    @test ParameterFunctionData(ParameterFunction(sin, vt, [1], [1])) isa ParameterFunctionData
end

# Test variable datatypes
@testset "Variables" begin
    # useful data
    m = InfiniteModel()
    num = Float64(0)
    sample_info = VariableInfo(true, num, true, num, true, num, true, num, true, true)
    func = (x) -> NaN
    inf_info = VariableInfo(true, num, true, num, true, num, false, func, true, true)
    pref = GeneralVariableRef(m, 1, IndependentParameterIndex, -1)
    vref = GeneralVariableRef(m, 1, InfiniteVariableIndex)
    # Infinite variable
    @test InfiniteVariable <: JuMP.AbstractVariable
    @test InfiniteVariable(inf_info, IC.VectorTuple(pref), [1], [1], true) isa InfiniteVariable
    # Semi-Infinite variable
    @test SemiInfiniteVariable <: JuMP.AbstractVariable
    @test SemiInfiniteVariable(vref, Dict(1 => Float64(2)), [1], [1]) isa SemiInfiniteVariable
    # Point variable
    @test PointVariable <: JuMP.AbstractVariable
    @test PointVariable(sample_info, vref, Float64[1]) isa PointVariable
    # VariableData
    @test VariableData <: AbstractDataObject
    @test VariableData(ScalarVariable(sample_info)) isa VariableData
end

# Test derivative datatypes
@testset "Derivatives" begin
    # useful data
    m = InfiniteModel()
    num = Float64(0)
    sample_info = VariableInfo(true, num, true, num, true, num, true, num, true, true)
    func = (x) -> NaN
    inf_info = VariableInfo(true, num, true, num, true, num, false, func, true, true)
    pref = GeneralVariableRef(m, 1, IndependentParameterIndex, -1)
    vref = GeneralVariableRef(m, 1, InfiniteVariableIndex)
    dref = GeneralVariableRef(m, 1, DerivativeIndex)
    # derivative
    @test Derivative <: JuMP.AbstractVariable
    @test Derivative(inf_info, true, vref, pref) isa Derivative
    # Semi-infinite derivative
    @test SemiInfiniteVariable(dref, Dict(1 => Float64(2)), [1], [1]) isa SemiInfiniteVariable
    # Point derivative
    @test PointVariable(sample_info, dref, Float64[1]) isa PointVariable
    # VariableData
    @test VariableData(Derivative(inf_info, true, vref, pref)) isa VariableData
end

# Test the measure datatypes
@testset "Measures" begin
    # Helper data
    m = InfiniteModel()
    pref = GeneralVariableRef(m, 2, IndependentParameterIndex)
    w(t) = 1
    info = NoGenerativeSupports()
    # abstract types
    @test AbstractMeasureData isa DataType
    # DiscreteMeasureData
    @test DiscreteMeasureData <: AbstractMeasureData
    @test DiscreteMeasureData(pref, ones(2), ones(2), All, w, NaN, NaN, false) isa DiscreteMeasureData
    @test DiscreteMeasureData([pref], ones(2), ones(1, 2), All, w, [NaN], [NaN], false) isa DiscreteMeasureData
    # FunctionalDistributionDomain
    @test FunctionalDiscreteMeasureData <: AbstractMeasureData
    @test FunctionalDiscreteMeasureData(pref, w, 2, All, info, w, NaN, NaN, false) isa FunctionalDiscreteMeasureData
    @test FunctionalDiscreteMeasureData([pref], w, 2, All, w, [NaN], [NaN], false) isa FunctionalDiscreteMeasureData
    @test FunctionalDiscreteMeasureData([pref], w, 2, All, info, w, [NaN], [NaN], false) isa FunctionalDiscreteMeasureData
    # Measure
    @test Measure isa UnionAll
    @test Measure(zero(AffExpr),
                  DiscreteMeasureData(pref, ones(2), ones(2), All, w, NaN, NaN, false),
                  [1], [1], false) isa Measure
    # MeasureData
    @test MeasureData <: AbstractDataObject
    @test MeasureData(Measure(zero(AffExpr), DiscreteMeasureData(pref,
                      ones(2), ones(2), All, w, NaN, NaN, false), [1], [1], true)) isa MeasureData
end

# Test the constraint datatypes
@testset "Constraints" begin
    # Setup
    m = InfiniteModel()
    pref = GeneralVariableRef(m, 2, IndependentParameterIndex, -1)
    rs = DomainRestrictions(pref => [0, 1])
    con = ScalarConstraint(zero(AffExpr), MOI.Integer())
    # DomainRestrictedConstraint
    @test DomainRestrictedConstraint <: JuMP.AbstractConstraint
    @test DomainRestrictedConstraint(con, rs).restrictions == rs
    # ConstraintData
    @test ConstraintData <: AbstractDataObject
    @test ConstraintData(con, [1], "", MeasureIndex[], false) isa ConstraintData
end

# Test the operator constructor
@testset "Nonlinear Operators" begin
    # Setup
    f(a) = a^3
    g(a::Int) = 42
    h(a, b) = 42
    f1(a) = 32
    f2(a) = 10
    h1(a, b) = 13
    function hg(v::AbstractVector, a, b)
        v[1] = 1
        v[2] = 2
        return 
    end
    function ∇²h(H, x...)
        H[1, 1] = 1200 * x[1]^2 - 400 * x[2] + 2
        H[2, 1] = -400 * x[1]
        H[2, 2] = 200.0
        return
    end
    # Test errors
    @test_throws ErrorException NLPOperator(:a, 1, f, g)
    @test_throws ErrorException NLPOperator(:a, 2, f, g)
    @test_throws ErrorException NLPOperator(:a, 1, f, g, f)
    @test_throws ErrorException NLPOperator(:a, 1, f, f, g)
    @test_throws ErrorException NLPOperator(:a, 2, h, hg, hg)
    # Test regular builds
    @test NLPOperator(:a, 1, f).name == :a
    @test NLPOperator(:a, 1, f).dim == 1
    @test NLPOperator(:a, 1, f).f == f
    @test NLPOperator(:a, 1, f, f) isa NLPOperator{typeof(f), typeof(f), Nothing}
    @test NLPOperator(:a, 1, f, f).∇f == f
    @test NLPOperator(:a, 1, f, f, f) isa NLPOperator{typeof(f), typeof(f), typeof(f)}
    @test NLPOperator(:a, 1, f, f, f).∇²f == f
    @test NLPOperator(:a, 2, h, hg) isa NLPOperator{typeof(h), typeof(hg), Nothing}
    @test NLPOperator(:a, 2, h, hg).∇f == hg
    @test NLPOperator(:a, 2, h, hg, ∇²h) isa NLPOperator{typeof(h), typeof(hg), typeof(∇²h)}
    @test NLPOperator(:a, 2, h, hg, ∇²h).∇²f == ∇²h
end
