# Test collection_sets function
@testset "collection_sets" begin
    set1 = IntervalSet(0, 1)
    set2 = UniDistributionSet(Uniform())
    set = CollectionSet([set1, set2])
    @test_throws ErrorException collection_sets(set1)
    @test collection_sets(set) == [set1, set2]
end

# Test length
@testset "Base.length" begin
    sset1 = IntervalSet(0, 1)
    sset2 = UniDistributionSet(Uniform())
    set1 = CollectionSet([sset1, sset2])
    set2 = MultiDistributionSet(MvNormal(ones(4)))
    @test length(sset1) == 1
    @test length(sset2) == 1
    @test length(set1) == 2
    @test length(set2) == 4
end

# Test supports_in_set functions
@testset "supports_in_set" begin
    # supports_in_set (IntervalSet)
    @testset "IntervalSet" begin
        set = IntervalSet(0, 1)
        @test supports_in_set(0, set)
        @test !supports_in_set(-1, set)
        @test !supports_in_set(2, set)
    end
    # supports_in_set (UnivariateDistribution set)
    @testset "Univariate Distribution" begin
        set = UniDistributionSet(Uniform())
        @test supports_in_set(0, set)
        @test !supports_in_set(-1, set)
        @test !supports_in_set(2, set)
    end
    # supports_in_set (MultivariateDistribution set)
    @testset "Multivariate Distribution" begin
        set = MultiDistributionSet(MvNormal(ones(2)))
        @test supports_in_set(ones(2), set)
        @test supports_in_set(ones(2, 10), set)
        bad_supports = [1 1; 2 2; 3 3];
        @test_throws ErrorException supports_in_set(bad_supports, set)
    end
    # supports_in_set (CollectionSet)
    @testset "CollectionSet" begin
        set1 = IntervalSet(0, 1)
        set2 = UniDistributionSet(Uniform())
        set = CollectionSet([set1, set2])
        bad_supports = [1 1; 2 2; 3 3];
        @test_throws ErrorException supports_in_set(bad_supports, set)
        supports = [1 2; 3 4];
        @test !supports_in_set(supports, set)
        supports = [0 1; 0 1];
        @test supports_in_set(supports, set)
    end
    # supports_in_set (Fallback)
    @testset "Fallback" begin
        @test supports_in_set(0, BadSet())
    end
end

# Test lower bound functions
@testset "Lower Bound" begin
    set1 = IntervalSet(0, 1)
    set2 = UniDistributionSet(Uniform())
    set3 = MultiDistributionSet(MvNormal([0., 0.], [1. 0.; 0. 2.]))
    set4 = CollectionSet([set1, set2])
    set5 = BadSet()
    # JuMP.has_lower_bound
    @testset "JuMP.has_lower_bound" begin
        @test has_lower_bound(set1)
        @test has_lower_bound(set2)
        @test !has_lower_bound(set3)
        @test has_lower_bound(set4)
        @test !has_lower_bound(set5)
    end
    # JuMP.lower_bound
    @testset "JuMP.lower_bound" begin
        @test lower_bound(set1) == 0.0
        @test lower_bound(set2) == 0.0
        @test_throws ErrorException lower_bound(set3)
        @test lower_bound(set4) == [0.0, 0.0]
        @test_throws ErrorException lower_bound(set5)
    end
    # JuMP.set_lower_bound
    @testset "JuMP.set_lower_bound" begin
        @test set_lower_bound(set1, 0.5) == IntervalSet(0.5, 1)
        @test_throws ErrorException set_lower_bound(set2, 2)
        @test_throws ErrorException set_lower_bound(set3, 2)
        @test_throws ErrorException set_lower_bound(set4, [0.1, 0.1])
        @test_throws ErrorException set_lower_bound(set5, 2)
        set6 = CollectionSet([set1, set1])
        @test set_lower_bound(set6, [0.1, 0.1]).sets ==
                                    [IntervalSet(0.1, 1.), IntervalSet(0.1, 1.)]
    end
end

# Test upper bound functions
@testset "Upper Bound" begin
    set1 = IntervalSet(0, 1)
    set2 = UniDistributionSet(Uniform())
    set3 = MultiDistributionSet(MvNormal([0., 0.], [1. 0.; 0. 2.]))
    set4 = CollectionSet([set1, set2])
    set5 = BadSet()
    # JuMP.has_upper_bound
    @testset "JuMP.has_upper_bound" begin
        @test has_upper_bound(set1)
        @test has_upper_bound(set2)
        @test !has_upper_bound(set3)
        @test has_upper_bound(set4)
        @test !has_upper_bound(set5)
    end
    # JuMP.upper_bound
    @testset "JuMP.upper_bound" begin
        @test upper_bound(set1) == 1.0
        @test upper_bound(set2) == 1.0
        @test_throws ErrorException upper_bound(set3)
        @test upper_bound(set4) == [1.0, 1.0]
        @test_throws ErrorException upper_bound(set5)
    end
    # JuMP.set_upper_bound
    @testset "JuMP.set_upper_bound" begin
        @test set_upper_bound(set1, 0.5) == IntervalSet(0, 0.5)
        @test_throws ErrorException set_upper_bound(set2, 2)
        @test_throws ErrorException set_upper_bound(set3, 2)
        @test_throws ErrorException set_upper_bound(set4, [0.1, 0.1])
        @test_throws ErrorException set_upper_bound(set5, 2)
        set6 = CollectionSet([set1, set1])
        @test set_upper_bound(set6, [0.9, 0.9]).sets ==
                                    [IntervalSet(0., 0.9), IntervalSet(0., 0.9)]

    end
end

# Test support generation
@testset "generate_support_values" begin
    @testset "IntervalSet" begin
        set = IntervalSet(0., 1.)
        @test generate_support_values(set, num_supports = 10, sig_digits = 3)[1] isa Vector{<:Number}
        @test generate_support_values(set, num_supports = 10, sig_digits = 3)[2] == UniformGrid
        @test generate_support_values(set, num_supports = 10, sig_digits = 3)[1][2] == 0.111
        @test generate_support_values(set, num_supports = 10, sig_digits = 3)[1][2] != 1/11
        @test length(generate_support_values(set, num_supports = 10, sig_digits = 3)[1]) == 10
        @test generate_support_values(set, Val(MCSample), num_supports = 10)[1] isa Vector{<:Number}
        @test generate_support_values(set, Val(MCSample), num_supports = 10)[2] == MCSample
        @test_throws ErrorException generate_support_values(set, Val(:a))
    end
    @testset "Distribution Sets" begin
        dist1 = Normal(0., 1.)
        dist2 = MvNormal([0.; 0.], [1. 0.; 0. 2.])
        set1 = UniDistributionSet(dist1)
        set2 = MultiDistributionSet(dist2)
        @test generate_support_values(set1, num_supports = 10)[1] isa Vector{<:Number}
        @test generate_support_values(set2, num_supports = 10)[1] isa Array{<:Number, 2}
        @test generate_support_values(set2, num_supports = 10)[2] == WeightedSample
        @test length(generate_support_values(set1, num_supports = 10)[1]) == 10
        @test size(generate_support_values(set2, num_supports = 10)[1]) == (2, 10)
        @test_throws ErrorException generate_support_values(set1, Val(:a))
        @test_throws ErrorException generate_support_values(set2, Val(:a))
    end
    @testset "Matrix Distribution Sets" begin
        dist = MatrixBeta(2, 2, 2)
        set = MultiDistributionSet(dist)
        @test generate_support_values(set, num_supports = 10)[1] isa Array{<:Number, 2}
        @test generate_support_values(set, num_supports = 10)[2] == WeightedSample
        @test size(generate_support_values(set, num_supports = 10)[1]) == (4, 10)
        @test_throws ErrorException generate_support_values(set, Val(:a))
    end
    @testset "_generate_collection_supports" begin
        set1 = IntervalSet(0., 1.)
        set2 = IntervalSet(0., 1.)
        set = CollectionSet([set1, set2])
        @test InfiniteOpt._generate_collection_supports(set, 10, 3) isa Array{<:Number, 2}
        @test InfiniteOpt._generate_collection_supports(set, 10, 3)[2, 2] == 0.111
        @test InfiniteOpt._generate_collection_supports(set, 10, 3)[2, 2] != 1/11
        @test size(InfiniteOpt._generate_collection_supports(set, 10, 3)) == (2, 10)
    end
    @testset "CollectionSet (IntervalSets)" begin
        set1 = IntervalSet(0., 1.)
        set2 = IntervalSet(0., 1.)
        set = CollectionSet([set1, set2])
        @test generate_support_values(set, num_supports = 10, sig_digits = 3)[1] isa Array{<:Number, 2}
        @test generate_support_values(set, num_supports = 10, sig_digits = 3)[2] == UniformGrid
        @test generate_support_values(set, num_supports = 10, sig_digits = 3)[1][2, 2] == 0.111
        @test generate_support_values(set, num_supports = 10, sig_digits = 3)[1][2, 2] != 1/11
        @test size(generate_support_values(set, num_supports = 10, sig_digits = 3)[1]) == (2, 10)
        @test_throws ErrorException generate_support_values(set, Val(:a))
    end
    @testset "CollectionSet (UniDistributionSets)" begin
        set1 = UniDistributionSet(Normal())
        set2 = UniDistributionSet(Normal())
        set = CollectionSet([set1, set2])
        @test generate_support_values(set, num_supports = 10, sig_digits = 3)[1] isa Array{<:Number, 2}
        @test generate_support_values(set, num_supports = 10, sig_digits = 3)[2] == WeightedSample
        @test size(generate_support_values(set, num_supports = 10, sig_digits = 3)[1]) == (2, 10)
        @test_throws ErrorException generate_support_values(set, Val(:a))
    end
    @testset "CollectionSet (InfiniteScalarSets)" begin
        set1 = UniDistributionSet(Normal())
        set2 = IntervalSet(0., 1.)
        set = CollectionSet([set1, set2])
        @test generate_support_values(set, num_supports = 10, sig_digits = 3)[1] isa Array{<:Number, 2}
        @test generate_support_values(set, num_supports = 10, sig_digits = 3)[2] == Mixture
        @test size(generate_support_values(set, num_supports = 10, sig_digits = 3)[1]) == (2, 10)
        @test_throws ErrorException generate_support_values(set, Val(:a))
    end
    @testset "Fallback" begin
        @test_throws ErrorException generate_support_values(BadSet())
    end
    @testset "User Interface" begin
        set = IntervalSet(0., 1.)
        @test generate_supports(set)[1] isa Vector{<:Number}
        @test generate_supports(set, MCSample)[1] isa Vector{<:Number}
    end
end
