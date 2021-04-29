# Test collection_domains function
@testset "collection_domains" begin
    domain1 = IntervalDomain(0, 1)
    domain2 = UniDistributionDomain(Uniform())
    domain = CollectionDomain([domain1, domain2])
    @test_throws ErrorException collection_domains(domain1)
    @test collection_domains(domain) == [domain1, domain2]
end

# Test length
@testset "Base.length" begin
    sdomain1 = IntervalDomain(0, 1)
    sdomain2 = UniDistributionDomain(Uniform())
    domain1 = CollectionDomain([sdomain1, sdomain2])
    domain2 = MultiDistributionDomain(MvNormal(ones(4)))
    @test length(sdomain1) == 1
    @test length(sdomain2) == 1
    @test length(domain1) == 2
    @test length(domain2) == 4
end

# Test supports_in_domain functions
@testset "supports_in_domain" begin
    # supports_in_domain (IntervalDomain)
    @testset "IntervalDomain" begin
        domain = IntervalDomain(0, 1)
        @test supports_in_domain(0, domain)
        @test !supports_in_domain(-1, domain)
        @test !supports_in_domain(2, domain)
    end
    # supports_in_domain (UnivariateDistribution domain)
    @testset "Univariate Distribution" begin
        domain = UniDistributionDomain(Uniform())
        @test supports_in_domain(0, domain)
        @test !supports_in_domain(-1, domain)
        @test !supports_in_domain(2, domain)
    end
    # supports_in_domain (MultivariateDistribution domain)
    @testset "Multivariate Distribution" begin
        domain = MultiDistributionDomain(MvNormal(ones(2)))
        @test supports_in_domain(ones(2), domain)
        @test supports_in_domain(ones(2, 10), domain)
        bad_supports = [1 1; 2 2; 3 3];
        @test_throws ErrorException supports_in_domain(bad_supports, domain)
    end
    # supports_in_domain (CollectionDomain)
    @testset "CollectionDomain" begin
        domain1 = IntervalDomain(0, 1)
        domain2 = UniDistributionDomain(Uniform())
        domain = CollectionDomain([domain1, domain2])
        bad_supports = [1 1; 2 2; 3 3];
        @test_throws ErrorException supports_in_domain(bad_supports, domain)
        supports = [1 2; 3 4];
        @test !supports_in_domain(supports, domain)
        supports = [0 1; 0 1];
        @test supports_in_domain(supports, domain)
    end
    # supports_in_domain (Fallback)
    @testset "Fallback" begin
        @test supports_in_domain(0, BadDomain())
    end
end

# Test lower bound functions
@testset "Lower Bound" begin
    domain1 = IntervalDomain(0, 1)
    domain2 = UniDistributionDomain(Uniform())
    domain3 = MultiDistributionDomain(MvNormal([0., 0.], [1. 0.; 0. 2.]))
    domain4 = CollectionDomain([domain1, domain2])
    domain5 = BadDomain()
    # JuMP.has_lower_bound
    @testset "JuMP.has_lower_bound" begin
        @test has_lower_bound(domain1)
        @test has_lower_bound(domain2)
        @test !has_lower_bound(domain3)
        @test has_lower_bound(domain4)
        @test !has_lower_bound(domain5)
    end
    # JuMP.lower_bound
    @testset "JuMP.lower_bound" begin
        @test lower_bound(domain1) == 0.0
        @test lower_bound(domain2) == 0.0
        @test_throws ErrorException lower_bound(domain3)
        @test lower_bound(domain4) == [0.0, 0.0]
        @test_throws ErrorException lower_bound(domain5)
    end
    # JuMP.set_lower_bound
    @testset "JuMP.set_lower_bound" begin
        @test set_lower_bound(domain1, 0.5) == IntervalDomain(0.5, 1)
        @test_throws ErrorException set_lower_bound(domain2, 2)
        @test_throws ErrorException set_lower_bound(domain3, 2)
        @test_throws ErrorException set_lower_bound(domain4, [0.1, 0.1])
        @test_throws ErrorException set_lower_bound(domain5, 2)
        domain6 = CollectionDomain([domain1, domain1])
        @test set_lower_bound(domain6, [0.1, 0.1]).domains ==
                                    [IntervalDomain(0.1, 1.), IntervalDomain(0.1, 1.)]
    end
end

# Test upper bound functions
@testset "Upper Bound" begin
    domain1 = IntervalDomain(0, 1)
    domain2 = UniDistributionDomain(Uniform())
    domain3 = MultiDistributionDomain(MvNormal([0., 0.], [1. 0.; 0. 2.]))
    domain4 = CollectionDomain([domain1, domain2])
    domain5 = BadDomain()
    # JuMP.has_upper_bound
    @testset "JuMP.has_upper_bound" begin
        @test has_upper_bound(domain1)
        @test has_upper_bound(domain2)
        @test !has_upper_bound(domain3)
        @test has_upper_bound(domain4)
        @test !has_upper_bound(domain5)
    end
    # JuMP.upper_bound
    @testset "JuMP.upper_bound" begin
        @test upper_bound(domain1) == 1.0
        @test upper_bound(domain2) == 1.0
        @test_throws ErrorException upper_bound(domain3)
        @test upper_bound(domain4) == [1.0, 1.0]
        @test_throws ErrorException upper_bound(domain5)
    end
    # JuMP.set_upper_bound
    @testset "JuMP.set_upper_bound" begin
        @test set_upper_bound(domain1, 0.5) == IntervalDomain(0, 0.5)
        @test_throws ErrorException set_upper_bound(domain2, 2)
        @test_throws ErrorException set_upper_bound(domain3, 2)
        @test_throws ErrorException set_upper_bound(domain4, [0.1, 0.1])
        @test_throws ErrorException set_upper_bound(domain5, 2)
        domain6 = CollectionDomain([domain1, domain1])
        @test set_upper_bound(domain6, [0.9, 0.9]).domains ==
                                    [IntervalDomain(0., 0.9), IntervalDomain(0., 0.9)]
    end
end

# Test support labels 
@testset "Support Labels" begin 
    @testset "DataTypes" begin 
        @test AbstractSupportLabel isa DataType 
        @test All <: AbstractSupportLabel
        @test PublicLabel <: AbstractSupportLabel
        @test UserDefined <: PublicLabel
        @test UniformGrid <: PublicLabel
        @test SampleLabel <: PublicLabel
        @test MCSample <: SampleLabel
        @test WeightedSample <: SampleLabel
        @test Mixture <: PublicLabel
        @test UniqueMeasure isa UnionAll
        @test InternalLabel <: AbstractSupportLabel
        @test InfiniteOpt._NoLabel <: AbstractSupportLabel
    end
    @testset "generate_unique_label" begin 
        @test InfiniteOpt.generate_unique_label() <: UniqueMeasure
    end
end

# Test support generation
@testset "generate_support_values" begin
    @testset "IntervalDomain" begin
        domain = IntervalDomain(0., 1.)
        @test generate_support_values(domain, num_supports = 10, sig_digits = 3)[1] isa Vector{<:Number}
        @test generate_support_values(domain, num_supports = 10, sig_digits = 3)[2] == UniformGrid
        @test generate_support_values(domain, num_supports = 10, sig_digits = 3)[1][2] == 0.111
        @test generate_support_values(domain, num_supports = 10, sig_digits = 3)[1][2] != 1/11
        @test length(generate_support_values(domain, num_supports = 10, sig_digits = 3)[1]) == 10
        @test generate_support_values(domain, MCSample, num_supports = 10)[1] isa Vector{<:Number}
        @test generate_support_values(domain, MCSample, num_supports = 10)[2] == MCSample
        @test_throws ErrorException generate_support_values(domain, Val(:a))
        @test generate_support_values(domain, All)[2] == UniformGrid
    end
    @testset "Distribution Domains" begin
        dist1 = Normal(0., 1.)
        dist2 = MvNormal([0.; 0.], [1. 0.; 0. 2.])
        domain1 = UniDistributionDomain(dist1)
        domain2 = MultiDistributionDomain(dist2)
        @test generate_support_values(domain1, num_supports = 10)[1] isa Vector{<:Number}
        @test generate_support_values(domain2, num_supports = 10)[1] isa Array{<:Number, 2}
        @test generate_support_values(domain2, num_supports = 10)[2] == WeightedSample
        @test length(generate_support_values(domain1, num_supports = 10)[1]) == 10
        @test size(generate_support_values(domain2, num_supports = 10)[1]) == (2, 10)
        @test_throws ErrorException generate_support_values(domain1, Val(:a))
        @test_throws ErrorException generate_support_values(domain2, Val(:a))
        @test generate_support_values(domain1, All)[2] == WeightedSample
        @test generate_support_values(domain2, All)[2] == WeightedSample
        @test generate_support_values(domain1, MCSample)[2] == MCSample
    end
    @testset "Matrix Distribution Domains" begin
        dist = MatrixBeta(2, 2, 2)
        domain = MultiDistributionDomain(dist)
        @test generate_support_values(domain, num_supports = 10)[1] isa Array{<:Number, 2}
        @test generate_support_values(domain, num_supports = 10)[2] == WeightedSample
        @test size(generate_support_values(domain, num_supports = 10)[1]) == (4, 10)
        @test_throws ErrorException generate_support_values(domain, Val(:a))
        @test generate_support_values(domain, All)[2] == WeightedSample
    end
    @testset "_generate_collection_supports" begin
        domain1 = IntervalDomain(0., 1.)
        domain2 = IntervalDomain(0., 1.)
        domain = CollectionDomain([domain1, domain2])
        @test InfiniteOpt._generate_collection_supports(domain, 10, 3) isa Array{<:Number, 2}
        @test InfiniteOpt._generate_collection_supports(domain, 10, 3)[2, 2] == 0.111
        @test InfiniteOpt._generate_collection_supports(domain, 10, 3)[2, 2] != 1/11
        @test size(InfiniteOpt._generate_collection_supports(domain, 10, 3)) == (2, 10)
    end
    @testset "CollectionDomain (IntervalDomains)" begin
        domain1 = IntervalDomain(0., 1.)
        domain2 = IntervalDomain(0., 1.)
        domain = CollectionDomain([domain1, domain2])
        @test generate_support_values(domain, num_supports = 10, sig_digits = 3)[1] isa Array{<:Number, 2}
        @test generate_support_values(domain, num_supports = 10, sig_digits = 3)[2] == UniformGrid
        @test generate_support_values(domain, num_supports = 10, sig_digits = 3)[1][2, 2] == 0.111
        @test generate_support_values(domain, num_supports = 10, sig_digits = 3)[1][2, 2] != 1/11
        @test size(generate_support_values(domain, num_supports = 10, sig_digits = 3)[1]) == (2, 10)
        @test_throws ErrorException generate_support_values(domain, Val(:a))
        @test generate_support_values(domain, All)[2] == UniformGrid
        @test generate_support_values(domain, MCSample)[2] == MCSample 
        @test generate_support_values(domain, MCSample)[1] isa Array{<:Number, 2}
    end
    @testset "CollectionDomain (UniDistributionDomains)" begin
        domain1 = UniDistributionDomain(Normal())
        domain2 = UniDistributionDomain(Normal())
        domain = CollectionDomain([domain1, domain2])
        @test generate_support_values(domain, num_supports = 10, sig_digits = 3)[1] isa Array{<:Number, 2}
        @test generate_support_values(domain, num_supports = 10, sig_digits = 3)[2] == WeightedSample
        @test size(generate_support_values(domain, num_supports = 10, sig_digits = 3)[1]) == (2, 10)
        @test_throws ErrorException generate_support_values(domain, Val(:a))
        @test generate_support_values(domain, All)[2] == WeightedSample
    end
    @testset "CollectionDomain (InfiniteScalarDomains)" begin
        domain1 = UniDistributionDomain(Normal())
        domain2 = IntervalDomain(0., 1.)
        domain = CollectionDomain([domain1, domain2])
        @test generate_support_values(domain, num_supports = 10, sig_digits = 3)[1] isa Array{<:Number, 2}
        @test generate_support_values(domain, num_supports = 10, sig_digits = 3)[2] == Mixture
        @test generate_support_values(CollectionDomain([domain2, domain2]), MCSample, num_supports = 10, sig_digits = 3)[2] == MCSample
        @test size(generate_support_values(domain, num_supports = 10, sig_digits = 3)[1]) == (2, 10)
        @test_throws ErrorException generate_support_values(domain, Val(:a))
        @test generate_support_values(domain, All)[2] == Mixture
        @test generate_support_values(domain, MCSample)[2] == MCSample 
        @test generate_support_values(domain, MCSample)[1] isa Array{<:Number, 2}
    end
    @testset "Fallback" begin
        @test_throws ErrorException generate_support_values(BadDomain())
    end
    @testset "User Interface" begin
        domain = IntervalDomain(0., 1.)
        @test generate_supports(domain)[1] isa Vector{<:Number}
        @test generate_supports(domain, MCSample)[1] isa Vector{<:Number}
    end
end
