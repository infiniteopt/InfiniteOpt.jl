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
    # InfiniteVariableIndex
    @test InfiniteVariableIndex <: ObjectIndex
    @test InfiniteVariableIndex(1).value == 1
    # ReducedVariableIndex
    @test ReducedVariableIndex <: ObjectIndex
    @test ReducedVariableIndex(1).value == 1
    # PointVariableIndex
    @test PointVariableIndex <: ObjectIndex
    @test PointVariableIndex(1).value == 1
    # HoldVariableIndex
    @test HoldVariableIndex <: ObjectIndex
    @test HoldVariableIndex(1).value == 1
    # DerivativeIndex
    @test DerivativeIndex <: ObjectIndex
    @test DerivativeIndex(1).value == 1
    # MeasureIndex
    @test MeasureIndex <: ObjectIndex
    @test MeasureIndex(1).value == 1
    # ConstraintIndex
    @test ConstraintIndex <: ObjectIndex
    @test ConstraintIndex(1).value == 1
    # test CleverDict extensions
    @testset "CleverDict Extensions" begin
        # index_to_key
        @test MOIUC.index_to_key(IndependentParameterIndex, Int64(42)) == IndependentParameterIndex(42)
        @test MOIUC.index_to_key(DependentParametersIndex, Int64(42)) == DependentParametersIndex(42)
        @test MOIUC.index_to_key(FiniteParameterIndex, Int64(42)) == FiniteParameterIndex(42)
        @test MOIUC.index_to_key(InfiniteVariableIndex, Int64(42)) == InfiniteVariableIndex(42)
        @test MOIUC.index_to_key(ReducedVariableIndex, Int64(42)) == ReducedVariableIndex(42)
        @test MOIUC.index_to_key(PointVariableIndex, Int64(42)) == PointVariableIndex(42)
        @test MOIUC.index_to_key(HoldVariableIndex, Int64(42)) == HoldVariableIndex(42)
        @test MOIUC.index_to_key(DerivativeIndex, Int64(42)) == DerivativeIndex(42)
        @test MOIUC.index_to_key(MeasureIndex, Int64(42)) == MeasureIndex(42)
        @test MOIUC.index_to_key(ConstraintIndex, Int64(42)) == ConstraintIndex(42)
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

# Test Infinite Sets
@testset "Infinite Sets" begin
    # Abstract types
    @test AbstractInfiniteSet isa DataType
    @test InfiniteScalarSet <: AbstractInfiniteSet
    @test InfiniteArraySet <: AbstractInfiniteSet
    # IntervalSet
    @test IntervalSet <: InfiniteScalarSet
    @test IntervalSet(0, 1) == IntervalSet(0.0, 1.0)
    @test_throws ErrorException IntervalSet(1, 0)
    # UniDistributionSet
    @test UniDistributionSet <: InfiniteScalarSet
    @test UniDistributionSet(Normal()) isa UniDistributionSet{<:Normal}
    # NonUnivariateDistribution
    @test InfiniteOpt.NonUnivariateDistribution <: Distribution
    # MultiDistributionSet
    @test MultiDistributionSet <: InfiniteArraySet
    @test MultiDistributionSet(Dirichlet([0.5, 0.5])) isa MultiDistributionSet{<:Dirichlet}
    # CollectionSet
    @test CollectionSet <: InfiniteArraySet
    @test CollectionSet([IntervalSet(0, 1), IntervalSet(2, 3)]) isa CollectionSet{IntervalSet}
end

# Test Derivative Evaluation Methods 
@testset "Derivative Methods" begin
    # abstract types
    @test AbstractDerivativeMethod isa DataType
    @test GenerativeDerivativeMethod <: AbstractDerivativeMethod
    @test NonGenerativeDerivativeMethod <: AbstractDerivativeMethod
    # test orthogonal collocation 
    @test OrthogonalCollocation <: GenerativeDerivativeMethod
    @test OrthogonalCollocation(2, Int, Int) isa OrthogonalCollocation
    # test finite difference
    @test FiniteDifference <: NonGenerativeDerivativeMethod
    @test FiniteDifference(Int) isa FiniteDifference
    # test support_label 
    @test_throws ErrorException support_label(FiniteDifference(Int))
    @test support_label(OrthogonalCollocation(2, Int, Int)) == Int
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
    method = FiniteDifference(Int)
    @test IndependentParameter(IntervalSet(0, 1), dict, 6, method).set isa IntervalSet
    # test FiniteParameter
    @test FiniteParameter <: ScalarParameter
    @test FiniteParameter(1) == FiniteParameter(1)
    # test DependentParameters
    @test DependentParameters <: InfOptParameter
    @test DependentParameters(CollectionSet([IntervalSet(0, 1)]),
                              Dict(zeros(1) => Set([All])), 6, [method]).set isa CollectionSet
    # test ScalarParameterData
    @test ScalarParameterData <: AbstractDataObject
    @test ScalarParameterData(FiniteParameter(42), 1, 1, "bob").name == "bob"
    # test MultiParameterData
    @test MultiParameterData <: AbstractDataObject
    params = DependentParameters(CollectionSet([IntervalSet(0, 1)]),
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
    @test InfiniteModel(caching_mode = MOIU.MANUAL) isa JuMP.AbstractModel
    @test InfiniteModel(mockoptimizer,
                        caching_mode = MOIU.MANUAL) isa JuMP.AbstractModel
    # test accessors
    @test InfiniteOpt._last_param_num(m) == 0
    @test InfiniteOpt._param_object_indices(m) isa Vector{Union{IndependentParameterIndex, DependentParametersIndex}}
end

# Test reference variable datatypes
@testset "References" begin
    m = InfiniteModel();
    # Abstract types
    @test DispatchVariableRef <: AbstractVariableRef
    @test MeasureFiniteVariableRef <: DispatchVariableRef
    @test FiniteVariableRef <: MeasureFiniteVariableRef
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
    @test FiniteParameterRef <: FiniteVariableRef
    idx = FiniteParameterIndex(1)
    @test FiniteParameterRef(m, idx) isa FiniteParameterRef
    # InfiniteVariableRef
    @test InfiniteVariableRef <: DispatchVariableRef
    idx = InfiniteVariableIndex(1)
    @test InfiniteVariableRef(m, idx) isa InfiniteVariableRef
    # ReducedVariableRef
    @test ReducedVariableRef <: DispatchVariableRef
    idx = ReducedVariableIndex(1)
    @test ReducedVariableRef(m, idx) isa ReducedVariableRef
    # PointVariableRef
    @test PointVariableRef <: FiniteVariableRef
    idx = PointVariableIndex(1)
    @test PointVariableRef(m, idx) isa PointVariableRef
    # HoldVariableRef
    @test HoldVariableRef <: FiniteVariableRef
    idx = HoldVariableIndex(1)
    @test HoldVariableRef(m, idx) isa HoldVariableRef
    # DerivativeRef
    @test DerivativeRef <: DispatchVariableRef
    idx = DerivativeIndex(1)
    @test DerivativeRef(m, idx) isa DerivativeRef
    # MeasureRef
    @test MeasureRef <: MeasureFiniteVariableRef
    idx = MeasureIndex(1)
    @test MeasureRef(m, idx) isa MeasureRef
    # InfOptConstraintRef
    @test InfOptConstraintRef isa UnionAll
    idx = ConstraintIndex(1)
    @test InfOptConstraintRef(m, idx, JuMP.ScalarShape()) isa InfOptConstraintRef
end

# Test ParameterBounds
@testset "Parameter Bounds" begin
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
        @test ParameterBounds isa UnionAll
        d = Dict(par3 => IntervalSet(0, 1))
        @test ParameterBounds(d) isa ParameterBounds{GeneralVariableRef}
        @test ParameterBounds() isa ParameterBounds{GeneralVariableRef}
    end
    # test _expand_parameter_tuple
    @testset "_expand_parameter_tuple" begin
        d = (par3 => IntervalSet(0, 1),)
        @test InfiniteOpt._expand_parameter_tuple(d) == Dict(d...)
    end
    # test _expand_parameter_dict(Dict{ParameterRef,IntervalSet}))
    @testset "_expand_parameter_dict (acceptable Form)" begin
        d = (par3 => IntervalSet(0, 1),)
        @test InfiniteOpt._expand_parameter_dict(d) == Dict(d...)
    end
    # test _expand_parameter_dict(Dict{Any,IntervalSet}))
    @testset "_expand_parameter_dict (Array Form)" begin
        d = (pars => IntervalSet(0, 1), par3 => IntervalSet(0, 1))
        @test isa(InfiniteOpt._expand_parameter_dict(d),
                  Dict{GeneralVariableRef, IntervalSet})
    end
    # test _expand_parameter_dict(Dict))
    @testset "_expand_parameter_dict (Fallback)" begin
        d = (pars => 1, par3 => 2)
        @test_throws ErrorException InfiniteOpt._expand_parameter_dict(d)
    end
    # test expansion definition
    @testset "ParameterBounds Expansion" begin
        d = (par3 => IntervalSet(0, 1),)
        @test ParameterBounds(d).intervals == Dict(d...)
        d = (pars => IntervalSet(0, 1), par3 => IntervalSet(0, 1))
        @test ParameterBounds(d).intervals isa Dict{GeneralVariableRef, IntervalSet}
    end
    pb = ParameterBounds((par3 => IntervalSet(0, 1),))
    # test intervals
    @testset "intervals" begin
        @test intervals(pb) == pb.intervals
    end
    # test Base.length
    @testset "Base.length" begin
        @test length(pb) == 1
    end
    # test Base.isempty
    @testset "Base.isempty" begin
        @test !isempty(pb)
        @test isempty(ParameterBounds())
    end
    # test Base.:(==)
    @testset "Base.:(==)" begin
        @test pb == ParameterBounds((par3 => IntervalSet(0, 1),))
        @test pb != ParameterBounds()
    end
    # test Base.copy
    @testset "Base.copy" begin
        @test copy(pb) == pb
    end
    # test Base.getindex
    @testset "Base.getindex" begin
        @test pb[par3] == IntervalSet(0, 1)
    end
    # test Base.setindex!
    @testset "Base.setindex!" begin
        @test (pb[par3] = IntervalSet(0, 2)) == IntervalSet(0, 2)
    end
    # test Base.haskey
    @testset "Base.haskey" begin
        @test haskey(pb, par3)
        @test !haskey(pb, par2)
    end
    # test Base.keys
    @testset "Base.keys" begin
        @test keys(pb) == keys(intervals(pb))
    end
    # test Base.iterate
    @testset "Base.iterate" begin
        @test [p for p in pb] == [par3 => IntervalSet(0, 2)]
    end
    # test Base.merge
    @testset "Base.merge" begin
        new_pb = ParameterBounds()
        @test merge(new_pb, pb) == pb
    end
    # test Base.merge!
    @testset "Base.merge!" begin
        new_pb = ParameterBounds()
        @test merge!(new_pb, pb) isa ParameterBounds
        @test new_pb == pb
    end
    # test Base.filter
    @testset "Base.filter" begin
        @test filter(e -> e[1] != par3, pb) == ParameterBounds()
    end
    # test Base.delete!
    @testset "Base.delete!" begin
        @test delete!(pb, par3) isa ParameterBounds
        @test length(pb) == 0
    end
end

# Test variable datatypes
@testset "Variables" begin
    # useful data
    m = InfiniteModel()
    num = Float64(0)
    sample_info = VariableInfo(true, num, true, num, true, num, true, num, true, true)
    func = (x) -> NaN
    inf_info = VariableInfo{Float64, Float64, Float64, Function}(true, num, true,
                 num, true, num, false, func, true, true)
    pref = GeneralVariableRef(m, 1, IndependentParameterIndex, -1)
    vref = GeneralVariableRef(m, 1, InfiniteVariableIndex)
    # Abstract variables
    @test InfOptVariable <: AbstractVariable
    # Infinite variable
    @test InfiniteVariable <: InfOptVariable
    @test InfiniteVariable(inf_info, IC.VectorTuple(pref), [1], [1], true) isa InfiniteVariable
    # Reduced variable
    @test ReducedVariable <: InfOptVariable
    @test ReducedVariable(vref, Dict(1 => Float64(2)), [1], [1]) isa ReducedVariable
    # Point variable
    @test PointVariable <: InfOptVariable
    @test PointVariable(sample_info, vref, Float64[1]) isa PointVariable
    # Hold variable
    @test HoldVariable <: InfOptVariable
    @test HoldVariable(sample_info, ParameterBounds()) isa HoldVariable
    # VariableData
    @test VariableData <: AbstractDataObject
    @test VariableData(HoldVariable(sample_info, ParameterBounds())) isa VariableData
end

# Test derivative datatypes
@testset "Derivatives" begin
    # useful data
    m = InfiniteModel()
    num = Float64(0)
    sample_info = VariableInfo(true, num, true, num, true, num, true, num, true, true)
    func = (x) -> NaN
    inf_info = VariableInfo{Float64, Float64, Float64, Function}(true, num, true,
                 num, true, num, false, func, true, true)
    pref = GeneralVariableRef(m, 1, IndependentParameterIndex, -1)
    vref = GeneralVariableRef(m, 1, InfiniteVariableIndex)
    dref = GeneralVariableRef(m, 1, DerivativeIndex)
    # derivative
    @test Derivative <: InfOptVariable
    @test Derivative(inf_info, true, vref, pref) isa Derivative
    # Reduced derivative
    @test ReducedVariable(dref, Dict(1 => Float64(2)), [1], [1]) isa ReducedVariable
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
    # abstract types
    @test AbstractMeasureData isa DataType
    # DiscreteMeasureData
    @test DiscreteMeasureData <: AbstractMeasureData
    @test DiscreteMeasureData(pref, ones(2), ones(2), All, w, NaN, NaN, false) isa DiscreteMeasureData
    @test DiscreteMeasureData([pref], ones(2), ones(1, 2), All, w, [NaN], [NaN], false) isa DiscreteMeasureData
    # FunctionalDistributionSet
    @test FunctionalDiscreteMeasureData <: AbstractMeasureData
    @test FunctionalDiscreteMeasureData(pref, w, 2, All, w, NaN, NaN, false) isa FunctionalDiscreteMeasureData
    @test FunctionalDiscreteMeasureData([pref], w, 2, All, w, [NaN], [NaN], false) isa FunctionalDiscreteMeasureData
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
    bounds = ParameterBounds(Dict(pref => IntervalSet(0, 1)))
    # Bounded constraints
    @test BoundedScalarConstraint <: JuMP.AbstractConstraint
    @test BoundedScalarConstraint(zero(AffExpr), MOI.Integer(),
                bounds, bounds).bounds.intervals[pref].lower_bound == 0.0
    # ConstraintData
    @test ConstraintData <: AbstractDataObject
    @test ConstraintData(ScalarConstraint(zero(AffExpr), MOI.Integer()), [1], "",
                         MeasureIndex[], false) isa ConstraintData
end
