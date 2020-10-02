@testset "Helper Methods" begin
    # helper method for Lobatto quadrature
    @testset "Lobatto Quadrature" begin
        @test InfiniteOpt._legendre(2) isa Polynomial{Float64}
        @test coeffs(InfiniteOpt._legendre(2)) == [-0.5, 0.0, 1.5]
        @test coeffs(InfiniteOpt._legendre(5)) == [0., 1.875, 0., -8.75, 0., 7.875]
        @test InfiniteOpt._compute_internal_node_basis(1) == [0.]
        @test isapprox(InfiniteOpt._compute_internal_node_basis(2), [-sqrt(5)/5, sqrt(5)/5])
    end

    m = InfiniteModel();
    @infinite_parameter(m, t in [0,1])
    @infinite_parameter(m, x[1:2] in [0,5])
    @infinite_variable(m, T(t))
    oc = OrthogonalCollocation(1, Lobatto)
    bad_oc = OrthogonalCollocation(1, TestOCTechnique)
    tref = dispatch_variable_ref(t)

    # test generate_derivative_supports
    @testset "generate_derivative_supports" begin
        # orthogonal collocation
        @test_throws ErrorException InfiniteOpt.generate_derivative_supports(tref, oc)
        add_supports(t, [0, 0.5, 1])
        @test InfiniteOpt.generate_derivative_supports(tref, oc) == [0.25, 0.75]
        @test_throws ErrorException InfiniteOpt.generate_derivative_supports(tref, bad_oc)
        
        # general fallbacks
        @test_throws ErrorException InfiniteOpt.generate_derivative_supports(tref, TestGenMethod())
        @test InfiniteOpt.generate_derivative_supports(tref, TestMethod()) == Float64[]
    end

    # test add_derivative_supports
    @testset "add_derivative_supports" begin
        set_derivative_method(t, oc)
        @test InfiniteOpt.add_derivative_supports(t) isa Nothing
        @test supports(t, label = All) == [0, 0.25, 0.5, 0.75, 1.]
        @test InfiniteOpt.add_derivative_supports(t) isa Nothing
        @test supports(t, label = All) == [0, 0.25, 0.5, 0.75, 1.]
                
        # fallback
        @test InfiniteOpt.add_derivative_supports(x[1]) isa Nothing
    end
end