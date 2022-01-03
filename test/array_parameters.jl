# Test core DispatchVariableRef extensions
@testset "Core Data Accessers" begin
    # Setup data
    m = InfiniteModel()
    obj_idx = DependentParametersIndex(1)
    idx = DependentParameterIndex(obj_idx, 1)
    domain = CollectionDomain([IntervalDomain(0, 1), IntervalDomain(0, 1)])
    methods = [InfiniteOpt.DefaultDerivativeMethod for i = 1:2]
    params = DependentParameters(domain, Dict{Vector{Float64}, Set{DataType}}(), 5, methods)
    object = MultiParameterData(params, 1, 1:2, ["p1", "p2"])
    pref = DependentParameterRef(m, idx)
    gvref = GeneralVariableRef(m, 1, DependentParameterIndex, 1)
    bad_idx = DependentParameterIndex(DependentParametersIndex(-1), 2)
    bad_pref = DependentParameterRef(m, bad_idx)
    # test dispatch_variable_ref
    @testset "dispatch_variable_ref" begin
        @test isequal(dispatch_variable_ref(m, idx), pref)
        @test isequal(dispatch_variable_ref(gvref), pref)
    end
    # test _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, object) == obj_idx
        @test InfiniteOpt._param_object_indices(m)[end] == obj_idx
    end
    # test _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(pref) === m.dependent_params
        @test InfiniteOpt._data_dictionary(gvref) === m.dependent_params
        @test InfiniteOpt._data_dictionary(m, DependentParameters) === m.dependent_params
    end
    # test _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(pref) == object
        @test InfiniteOpt._data_object(gvref) == object
        @test_throws ErrorException InfiniteOpt._data_object(bad_pref)
    end
    # test _core_variable_object
    @testset "_core_variable_object" begin
        @test InfiniteOpt._core_variable_object(pref) == params
        @test InfiniteOpt._core_variable_object(gvref) == params
    end
    # test _num_parameters
    @testset "_num_parameters" begin
        @test InfiniteOpt._num_parameters(pref) == 2
    end
    # test _delete_data_object
    @testset "_delete_data_object" begin
        @test is_valid(m, pref)
        @test InfiniteOpt._delete_data_object(pref) isa Nothing
        @test !is_valid(m, pref)
    end
end

# Test Definition
@testset "Definition" begin
    # Setup data
    m = InfiniteModel();
    sdomain1 = IntervalDomain(0, 1)
    sdomain2 = UniDistributionDomain(Uniform())
    domain1 = CollectionDomain([sdomain1, sdomain1])
    domain2 = MultiDistributionDomain(MvNormal(Diagonal(map(abs2, ones(2)))))
    domain3 = MultiDistributionDomain(MatrixBeta(2, 2, 2))
    domain4 = CollectionDomain([sdomain1, sdomain2])
    inds1 = IC.indices(ones(2))
    inds2 = IC.indices(ones(2, 2))
    inds3 = IC.indices(JuMPC.SparseAxisArray(Dict((i,) => domain2 for i in 1:2)))
    method = InfiniteOpt.DefaultDerivativeMethod
    # test _make_array_domain
    @testset "_make_array_domain" begin
        # multi-dimensional distributions 
        @test InfiniteOpt._make_array_domain(error, [domain2, domain2], inds1) == domain2
        @test InfiniteOpt._make_array_domain(error, [domain3 for i in 1:4], inds2) == domain3
        @test_throws ErrorException InfiniteOpt._make_array_domain(error, [domain3], inds1)
        @test_throws ErrorException InfiniteOpt._make_array_domain(error, [domain2, domain2], inds3)
        @test_throws ErrorException InfiniteOpt._make_array_domain(error, [domain2, domain3], inds1)
        # collection domains 
        @test InfiniteOpt._make_array_domain(error, [domain4, domain4], inds1) == domain4
        @test_throws ErrorException InfiniteOpt._make_array_domain(error, [domain4, domain1], inds1)
        @test_throws ErrorException InfiniteOpt._make_array_domain(error, [domain4, domain4], inds2)
        @test_throws ErrorException InfiniteOpt._make_array_domain(error, [domain4, domain4], inds3)
        # InfiniteArrayDomains
        @test InfiniteOpt._make_array_domain(error, [BadArrayDomain() for i in 1:2], inds1) == BadArrayDomain()
        @test_throws ErrorException InfiniteOpt._make_array_domain(error, [domain1, domain2], inds1)
        # scalar domains 
        @test InfiniteOpt._make_array_domain(error, [sdomain1, sdomain1], inds1).domains == domain1.domains
        @test InfiniteOpt._make_array_domain(error, [sdomain1, sdomain2], inds1).domains == domain4.domains
        # fallback 
        @test_throws ErrorException InfiniteOpt._make_array_domain(error, domain2, inds1)
    end
    # test _process_supports
    @testset "_process_supports" begin
        # single support
        supps = Dict([0., 0.] => Set([UserDefined]))
        @test InfiniteOpt._process_supports(error, [0, 0], domain1, 2) == supps
        @test_throws ErrorException InfiniteOpt._process_supports(error, [2, 2], domain1, 2)
        # multiple supports 
        supps = Dict([0., 0.] => Set([UserDefined]), [1., 1.] => Set([UserDefined]))
        @test InfiniteOpt._process_supports(error, [[0, 1], [0, 1]], domain1, 2) == supps
        @test_throws ErrorException InfiniteOpt._process_supports(error, [[0, 0], [1]], domain1, 2)
        @test_throws ErrorException InfiniteOpt._process_supports(error, [[0, 2], [0, 2]], domain1, 2)
    end
    # test _process_derivative_methods
    @testset "_process_derivative_methods" begin
        same = [method, method]
        domains = [domain1, domain1]
        @test InfiniteOpt._process_derivative_methods(error, same, domains) == same
        mix = [method, TestMethod()]
        @test InfiniteOpt._process_derivative_methods(error, mix, domains) == mix
        @test InfiniteOpt._process_derivative_methods(error, [], domains) == same
        @test_throws ErrorException InfiniteOpt._process_derivative_methods(error, [1, 2], domains)
    end
    # test _build_parameters
    @testset "_build_parameters" begin
        # test errors
        @test_throws ErrorException InfiniteOpt._build_parameters(error, [domain1, domain1], inds1, bad = 42)
        @test_throws ErrorException InfiniteOpt._build_parameters(error, [domain1, domain2], inds1)
        @test_throws ErrorException InfiniteOpt._build_parameters(error, [domain2, domain2], inds2)
        @test_throws ErrorException InfiniteOpt._build_parameters(error, [domain1, domain1], inds1, supports = [[0, 0], [2]])
        @test_throws ErrorException InfiniteOpt._build_parameters(error, [domain1, domain1], inds1, derivative_method = [2, 2])
        # test vector 
        @test InfiniteOpt._build_parameters(error, [domain1, domain1], inds1).domain == domain1 
        @test InfiniteOpt._build_parameters(error, [domain1, domain1], inds1).supports == Dict{Vector{Float64}, Set{DataType}}()
        @test InfiniteOpt._build_parameters(error, [domain1, domain1],  inds1, supports = [0, 0]).supports == Dict([0., 0.] => Set([UserDefined]))
        @test InfiniteOpt._build_parameters(error, [domain1, domain1], inds1).sig_digits isa Int 
        @test InfiniteOpt._build_parameters(error, [domain1, domain1], inds1).derivative_methods == [method, method]
        # test array 
        @test InfiniteOpt._build_parameters(error, [domain3 for i in 1:4], inds2).domain == domain3 
        @test InfiniteOpt._build_parameters(error, [domain3 for i in 1:4], inds2, supports = [[0, 1] for i in 1:4]).supports == Dict(zeros(4) => Set([UserDefined]), ones(4) => Set([UserDefined]))
        @test length(InfiniteOpt._build_parameters(error, [domain3 for i in 1:4], inds2, num_supports = 2).supports) == 2
        @test InfiniteOpt._build_parameters(error, [domain3 for i in 1:4], inds2, sig_digits = 3).sig_digits == 3
        @test InfiniteOpt._build_parameters(error, [domain3 for i in 1:4], inds2, derivative_method = [TestMethod() for i in 1:4]).derivative_methods isa Vector{TestMethod} 
    end
    # test add_parameters
    @testset "add_parameters" begin
        # test default
        prefs = [GeneralVariableRef(m, 1, DependentParameterIndex, i) for i in 1:2]
        params = InfiniteOpt._build_parameters(error, [domain1, domain1], inds1)
        @test isequal(add_parameters(m, params, ["", ""]), prefs)
        @test InfiniteOpt._param_object_indices(m) == [index(prefs[1]).object_index]
        @test InfiniteOpt._last_param_num(m) == 2
        @test name.(prefs) == ["", ""]
        # test vector build
        prefs = [GeneralVariableRef(m, 2, DependentParameterIndex, i) for i in 1:2]
        params = InfiniteOpt._build_parameters(error, [domain2, domain2], inds1)
        @test isequal(add_parameters(m, params, ["p1", "p2"]), prefs)
        @test InfiniteOpt._param_object_indices(m)[2] == index(prefs[1]).object_index
        @test InfiniteOpt._last_param_num(m) == 4
        @test name.(prefs) == ["p1", "p2"]
        # test array build
        prefs = [GeneralVariableRef(m, 3, DependentParameterIndex, i) for i in 1:4]
        params = InfiniteOpt._build_parameters(error, [domain3 for i in 1:4], inds2)
        @test isequal(add_parameters(m, params, ["p$i" for i in 1:4]), prefs)
        @test InfiniteOpt._param_object_indices(m)[3] == index(prefs[1]).object_index
        @test InfiniteOpt._last_param_num(m) == 8
        @test name.(prefs) == ["p$i" for i in 1:4]
        # test name error 
        @test_throws ErrorException add_parameters(m, params, ["p"])
    end
end

# test macro definition
@testset "Macro Definition" begin
    # Setup data
    m = InfiniteModel();
    dist1 = Uniform()
    dist2 = MvNormal(Diagonal(map(abs2, ones(2))))
    dist3 = MatrixBeta(2, 2, 2)
    sdomain1 = IntervalDomain(0, 1)
    sdomain2 = UniDistributionDomain(dist1)
    domain1 = CollectionDomain([sdomain1, sdomain1])
    domain2 = MultiDistributionDomain(dist2)
    domain3 = MultiDistributionDomain(dist3)
    domain4 = CollectionDomain([sdomain1, sdomain2])
    domain5 = CollectionDomain([sdomain2, sdomain2])
    # test @infinite_parameter
    @testset "@infinite_parameter" begin
        # test macro errors
        @test_macro_throws ErrorException @infinite_parameter(m)
        @test_macro_throws ErrorException @infinite_parameter(m, p[1:2], 42)
        @test_macro_throws ErrorException @infinite_parameter(m, param)
        @test_macro_throws ErrorException @infinite_parameter(m, "bob"[1:2])
        @test_macro_throws ErrorException @infinite_parameter(m, [1:2] in [0, 1], 
                                                              independent = a)
        @test_macro_throws ErrorException @infinite_parameter(m, [1:2] in [0, 1], 
                                                              distribution = dist1)
        @test_macro_throws ErrorException @infinite_parameter(m, [i = 1:2] in [0, 1], 
                                                              num_supports = i)
        # test domain errors
        @test_macro_throws ErrorException @infinite_parameter(m, p[1:2] in 42)
        @test_macro_throws ErrorException @infinite_parameter(m, p[1:2] in dist1)
        @test_macro_throws ErrorException @infinite_parameter(m, p[1:2] ~ [0, 1])
        # test build errors
        @test_macro_throws ErrorException @infinite_parameter(m, p[1:2] ~ dist1, bob = 42)
        @test_macro_throws ErrorException @infinite_parameter(m, p[i = 1:2] in [domain1, domain2][i])
        @test_macro_throws ErrorException @infinite_parameter(m, p[1:2] ~ dist3)
        @test_macro_throws ErrorException @infinite_parameter(m, p[i = 1:3] in domain1)
        @test_macro_throws ErrorException @infinite_parameter(m, p[1:2] in domain1, supports = 4)
        @test_macro_throws ErrorException @infinite_parameter(m, p[i = 1:2] in domain1, supports = [[1, 0], 1][i])
        # test simple explict build
        prefs = [GeneralVariableRef(m, 1, DependentParameterIndex, i) for i in 1:2]
        @test isequal(@infinite_parameter(m, a[1:2] ~ dist1, num_supports = 10), prefs)
        @test name.(prefs) == ["a[1]", "a[2]"]
        @test length(InfiniteOpt._core_variable_object(prefs[1]).supports) == 10
        @test WeightedSample in first(InfiniteOpt._core_variable_object(prefs[1]).supports)[2]
        @test InfiniteOpt._core_variable_object(prefs[1]).domain.domains == domain5.domains
        # test another explicit build
        prefs = [GeneralVariableRef(m, 2, DependentParameterIndex, i) for i in 1:2]
        expected = JuMPC.DenseAxisArray(prefs, 3:4)
        @test isequal(@infinite_parameter(m, b[3:4] in domain2, supports = 0), expected)
        @test name.(prefs) == ["b[3]", "b[4]"]
        @test InfiniteOpt._core_variable_object(prefs[1]).supports == Dict{Vector{Float64}, Set{DataType}}(zeros(2) => Set([UserDefined]))
        # test explicit build with some args
        prefs = [GeneralVariableRef(m, 3, DependentParameterIndex, i) for i in 1:2]
        expected = convert(JuMPC.SparseAxisArray, prefs)
        @test all(isequal.(@infinite_parameter(m, c[1:2] in sdomain1, supports = 0, 
                                  base_name = "z", 
                                  container = SparseAxisArray), expected))
        @test name.(prefs) == ["z[1]", "z[2]"]
        # test explicit build again
        prefs = [GeneralVariableRef(m, 4, DependentParameterIndex, i) for i in 1:2]
        @test isequal(@infinite_parameter(m, d[1:2] in [0, 1], num_supports = 10), prefs)
        @test name.(prefs) == ["d[1]", "d[2]"]
        @test length(InfiniteOpt._core_variable_object(prefs[1]).supports) == 10
        @test UniformGrid in first(InfiniteOpt._core_variable_object(prefs[1]).supports)[2]
        @test InfiniteOpt._core_variable_object(prefs[1]).domain.domains == domain1.domains
        # test test anonymous
        prefs = [GeneralVariableRef(m, 5, DependentParameterIndex, i) for i in 1:4]
        prefs = reshape(prefs, (2, 2))
        @test isequal(@infinite_parameter(m, [1:2, 1:2], distribution = dist3), prefs)
        @test name.(prefs) == ["" ""; "" ""]
        @test isempty(InfiniteOpt._core_variable_object(prefs[1]).supports)
        # test anonymous with domain keyword
        prefs = [GeneralVariableRef(m, 6, DependentParameterIndex, i) for i in 1:2]
        @test isequal(@infinite_parameter(m, [1:2], domain = sdomain2, 
                                  derivative_method = TestMethod()), prefs)
        # test anonymous with dist keyword and base_name
        prefs = [GeneralVariableRef(m, 7, DependentParameterIndex, i) for i in 1:2]
        @test isequal(@infinite_parameter(m, [1:2],
                                  distribution = dist1, base_name = "zz"), prefs)
        @test name.(prefs) == ["zz[1]", "zz[2]"]
    end
end

# test naming stuff
@testset "Naming" begin
    # setup data
    m = InfiniteModel();
    gvrefs = @infinite_parameter(m, a[1:2] in [0, 1])
    prefs = dispatch_variable_ref.(gvrefs)
    bad_idx = DependentParameterIndex(DependentParametersIndex(-1), 2)
    bad_pref = DependentParameterRef(m, bad_idx)
    # test _param_index
    @testset "_param_index" begin
        @test InfiniteOpt._param_index(prefs[1]) == 1
        @test InfiniteOpt._param_index(prefs[2]) == 2
    end
    # test _update_param_name_dict
    @testset "_update_param_name_dict" begin
        m.name_to_param = Dict{String, AbstractInfOptIndex}()
        @test InfiniteOpt._update_param_name_dict(m, m.dependent_params) isa Nothing
        @test m.name_to_param["a[1]"] == DependentParameterIndex(DependentParametersIndex(1), 1)
        @test m.name_to_param["a[2]"] == DependentParameterIndex(DependentParametersIndex(1), 2)
        m.name_to_param = nothing
    end
    # parameter_by_name
    @testset "parameter_by_name" begin
        @test isequal(parameter_by_name(m, "a[1]"), gvrefs[1])
        @test isequal(parameter_by_name(m, "a[2]"), gvrefs[2])
        @test isa(parameter_by_name(m, "a[3]"), Nothing)
        @infinite_parameter(m, b[2:3] in [0, 1], base_name = "a")
        @test_throws ErrorException parameter_by_name(m, "a[2]")
    end
    # test name
    @testset "JuMP.name" begin
        @test name(prefs[1]) == "a[1]"
        @test name(prefs[2]) == "a[2]"
        @test name(gvrefs[1]) == "a[1]"
        @test name(gvrefs[2]) == "a[2]"
        @test name(bad_pref) == ""
    end
    # test set_name
    @testset "JuMP.set_name" begin
        @test set_name(prefs[1], "joe") isa Nothing
        @test name(prefs[1]) == "joe"
        @test set_name(prefs[2], "joe") isa Nothing
        @test name(prefs[2]) == "joe"
        @test set_name(gvrefs[2], "joe2") isa Nothing
        @test name(gvrefs[2]) == "joe2"
    end
    # test is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, prefs[1])
        @test !is_valid(InfiniteModel(), prefs[1])
        @test is_valid(m, gvrefs[1])
    end
end

# test dependency functions
@testset "Parameter Dependencies" begin
    # setup data
    m = InfiniteModel();
    gvrefs = @infinite_parameter(m, a[1:2] in [0, 1])
    prefs = dispatch_variable_ref.(gvrefs)
    data = InfiniteOpt._data_object(first(prefs))
    bad_idx = DependentParameterIndex(DependentParametersIndex(-1), 2)
    bad_pref = DependentParameterRef(m, bad_idx)
    # test _infinite_variable_dependencies
    @testset "_infinite_variable_dependencies" begin
        @test InfiniteOpt._infinite_variable_dependencies(prefs[1]) == data.infinite_var_indices
        @test InfiniteOpt._infinite_variable_dependencies(prefs[2]) == data.infinite_var_indices
        @test InfiniteOpt._infinite_variable_dependencies(gvrefs[1]) == data.infinite_var_indices
        @test_throws ErrorException InfiniteOpt._infinite_variable_dependencies(bad_pref)
    end
    # test _parameter_function_dependencies
    @testset "_parameter_function_dependencies" begin
        @test InfiniteOpt._parameter_function_dependencies(prefs[1]) == data.parameter_func_indices
        @test InfiniteOpt._parameter_function_dependencies(prefs[2]) == data.parameter_func_indices
        @test InfiniteOpt._parameter_function_dependencies(gvrefs[1]) == data.parameter_func_indices
        @test_throws ErrorException InfiniteOpt._parameter_function_dependencies(bad_pref)
    end
    # test _measure_dependencies
    @testset "_measure_dependencies" begin
        @test InfiniteOpt._measure_dependencies(prefs[1]) == data.measure_indices[1]
        @test InfiniteOpt._measure_dependencies(prefs[2]) == data.measure_indices[2]
        @test InfiniteOpt._measure_dependencies(gvrefs[1]) == data.measure_indices[1]
        @test_throws ErrorException InfiniteOpt._measure_dependencies(bad_pref)
    end
    # test _constraint_dependencies
    @testset "_constraint_dependencies" begin
        @test InfiniteOpt._constraint_dependencies(prefs[1]) == data.constraint_indices[1]
        @test InfiniteOpt._constraint_dependencies(prefs[2]) == data.constraint_indices[2]
        @test InfiniteOpt._constraint_dependencies(gvrefs[1]) == data.constraint_indices[1]
        @test_throws ErrorException InfiniteOpt._constraint_dependencies(bad_pref)
    end
    # test _derivative_dependencies
    @testset "_derivative_dependencies" begin
        @test InfiniteOpt._derivative_dependencies(prefs[1]) == data.derivative_indices[1]
        @test InfiniteOpt._derivative_dependencies(prefs[2]) == data.derivative_indices[2]
        @test InfiniteOpt._derivative_dependencies(gvrefs[1]) == data.derivative_indices[1]
        @test_throws ErrorException InfiniteOpt._derivative_dependencies(bad_pref)
    end
    # test used_by_infinite_variable
    @testset "used_by_infinite_variable" begin
        # test not used
        @test !used_by_infinite_variable(prefs[1])
        @test !used_by_infinite_variable(prefs[2])
        @test !used_by_infinite_variable(gvrefs[1])
        # test used
        push!(data.infinite_var_indices, InfiniteVariableIndex(1))
        @test used_by_infinite_variable(prefs[1])
        # undo changes
        empty!(data.infinite_var_indices)
    end
    # test used_by_parameter_function
    @testset "used_by_parameter_function" begin
        # test not used
        @test !used_by_parameter_function(prefs[1])
        @test !used_by_parameter_function(prefs[2])
        @test !used_by_parameter_function(gvrefs[1])
        # test used
        push!(data.parameter_func_indices, ParameterFunctionIndex(1))
        @test used_by_parameter_function(prefs[1])
        # undo changes
        empty!(data.parameter_func_indices)
    end
    # test used_by_derivative
    @testset "used_by_derivative" begin
        # test not used
        @test !used_by_derivative(prefs[1])
        @test !used_by_derivative(prefs[2])
        @test !used_by_derivative(gvrefs[1])
        # test used
        push!(data.derivative_indices[1], DerivativeIndex(1))
        @test used_by_derivative(prefs[1])
        # undo changes
        empty!(data.derivative_indices[1])
    end
    # test used_by_measure
    @testset "used_by_measure" begin
        # test not used
        @test !used_by_measure(prefs[1])
        @test !used_by_measure(prefs[2])
        @test !used_by_measure(gvrefs[1])
        # test used
        push!(data.measure_indices[1], MeasureIndex(1))
        @test used_by_measure(prefs[1])
        # undo changes
        empty!(data.measure_indices[1])
    end
    # test used_by_constraint
    @testset "used_by_constraint" begin
        # test not used
        @test !used_by_constraint(prefs[1])
        @test !used_by_constraint(prefs[2])
        @test !used_by_constraint(gvrefs[1])
        # test used
        push!(data.constraint_indices[1], InfOptConstraintIndex(1))
        @test used_by_constraint(prefs[1])
        # undo changes
        empty!(data.constraint_indices[1])
    end
    # test used_by_objective
    @testset "used_by_objective" begin
        @test !used_by_objective(prefs[1])
        @test !used_by_objective(prefs[2])
        @test !used_by_objective(gvrefs[1])
    end
    # test is_used
    @testset "is_used" begin
        # test not used
        @test !is_used(prefs[1])
        @test !is_used(prefs[2])
        @test !is_used(gvrefs[1])
        # test used
        push!(data.infinite_var_indices, InfiniteVariableIndex(1))
        @test is_used(prefs[1])
    end
end

# test Parameter Object Methods
@testset "Object Methods" begin
    # setup data
    m = InfiniteModel();
    gvrefs = @infinite_parameter(m, a[1:2] in [0, 1])
    prefs = dispatch_variable_ref.(gvrefs)
    data = InfiniteOpt._data_object(first(prefs))
    domain = CollectionDomain([IntervalDomain(0, 2), IntervalDomain(0, 2)])
    methods = [InfiniteOpt.DefaultDerivativeMethod for i = 1:2]
    params = DependentParameters(domain, Dict{Vector{Float64}, Set{DataType}}(), 10, methods)
    bad_idx = DependentParameterIndex(DependentParametersIndex(-1), 2)
    bad_pref = DependentParameterRef(m, bad_idx)
    # test _parameter_number
    @testset "_parameter_number" begin
        @test InfiniteOpt._parameter_number(prefs[1]) == 1
        @test InfiniteOpt._parameter_number(prefs[2]) == 2
        @test InfiniteOpt._parameter_number(gvrefs[1]) == 1
        @test_throws ErrorException InfiniteOpt._parameter_number(bad_pref)
    end
    # test _parameter_numbers
    @testset "_parameter_numbers" begin
        @test InfiniteOpt._parameter_numbers(prefs[1]) == [1]
        @test InfiniteOpt._parameter_numbers(prefs[2]) == [2]
        @test InfiniteOpt._parameter_numbers(gvrefs[1]) == [1]
    end
    # test _object_number
    @testset "_object_number" begin
        @test InfiniteOpt._object_number(prefs[1]) == 1
        @test InfiniteOpt._object_number(prefs[2]) == 1
        @test InfiniteOpt._object_number(gvrefs[1]) == 1
    end
    # test _object_numbers
    @testset "_object_numbers" begin
        @test InfiniteOpt._object_numbers(prefs[1]) == [1]
        @test InfiniteOpt._object_numbers(prefs[2]) == [1]
        @test InfiniteOpt._object_numbers(gvrefs[1]) == [1]
    end
    # test _adaptive_data_update
    @testset "_adaptive_data_update" begin
        # test with same data 
        @test InfiniteOpt._adaptive_data_update(prefs[1], params, data) isa Nothing
        # test with different data 
        domain = MultiDistributionDomain(MvNormal([0, 0], LinearAlgebra.Diagonal(map(abs2, [1, 1]))))
        ps = DependentParameters(domain, Dict{Vector{Float64}, Set{DataType}}(), 10, methods)
        @test InfiniteOpt._adaptive_data_update(prefs[2], ps, data) isa Nothing
        @test InfiniteOpt._core_variable_object(prefs[2]) == ps
    end
    # test _set_core_variable_object
    @testset "_set_core_variable_object" begin
        # test with different type
        @test InfiniteOpt._set_core_variable_object(prefs[1], params) isa Nothing
        @test InfiniteOpt._core_variable_object(prefs[2]) == params
        # test with same type
        @test InfiniteOpt._set_core_variable_object(prefs[2], params) isa Nothing
        @test InfiniteOpt._core_variable_object(prefs[1]) == params
    end
end

# test Derivative Method 
@testset "Derivative Methods" begin 
    # setup the data
    m = InfiniteModel();
    gvrefs = @infinite_parameter(m, a[1:2] in [0, 1], derivative_method = TestMethod())
    prefs = dispatch_variable_ref.(gvrefs)
    push!(InfiniteOpt._constraint_dependencies(prefs[1]), InfOptConstraintIndex(1))
    # test has_generative_supports
    @testset "has_generative_supports" begin
        @test !has_generative_supports(prefs[1])
        @test !has_generative_supports(gvrefs[2])
    end
    # test add_generative_supports
    @testset "add_generative_supports" begin
        @test add_generative_supports(prefs[1]) isa Nothing
        @test add_generative_supports(gvrefs[2]) isa Nothing
    end
    # test _set_has_generative_supports
    @testset "_set_has_generative_supports" begin
        @test InfiniteOpt._set_has_generative_supports(prefs[1], true) isa Nothing
        @test InfiniteOpt._set_has_generative_supports(gvrefs[2], false) isa Nothing
    end
    # test has_derivative_constraints
    @testset "has_derivative_constraints" begin
        @test !has_derivative_constraints(prefs[1])
        InfiniteOpt._data_object(prefs[1]).has_deriv_constrs[1] = true
        @test has_derivative_constraints(gvrefs[1])
    end
    # test _set_has_derivative_constraints
    @testset "_set_has_derivative_constraints" begin
        @test InfiniteOpt._set_has_derivative_constraints(prefs[1], true) isa Nothing
        @test has_derivative_constraints(prefs[1])
        @test InfiniteOpt._set_has_derivative_constraints(gvrefs[1], false) isa Nothing
        @test !has_derivative_constraints(prefs[1])
    end
    # test _derivative_methods
    @testset "_derivative_methods" begin
        @test InfiniteOpt._derivative_methods(prefs[1]) isa Vector{TestMethod}
    end
    # test derivative_method
    @testset "derivative_method" begin
        @test derivative_method(prefs[1]) isa TestMethod
        @test derivative_method(gvrefs[2]) isa TestMethod
    end
    # test _adaptive_method_update
    @testset "_adaptive_method_update" begin
        # same type
        params = InfiniteOpt._core_variable_object(prefs[1])
        @test InfiniteOpt._adaptive_method_update(prefs[1], params, TestMethod()) isa Nothing
        @test derivative_method(prefs[1]) isa TestMethod
        # different type
        params = InfiniteOpt._core_variable_object(prefs[1])
        method = InfiniteOpt.DefaultDerivativeMethod
        @test InfiniteOpt._adaptive_method_update(prefs[2], params, method) isa Nothing
        @test InfiniteOpt._derivative_methods(prefs[1]) == [TestMethod(), method]
    end
    # test set_derivative_method
    @testset "set_derivative_method" begin
        # same type
        method = InfiniteOpt.DefaultDerivativeMethod
        @test set_derivative_method(prefs[1], method) isa Nothing
        @test derivative_method(prefs[1]) == method
        @test InfiniteOpt._derivative_methods(prefs[1]) == [method, method]
        # different type
        @test set_derivative_method(prefs[2], TestMethod()) isa Nothing
        @test InfiniteOpt._derivative_methods(prefs[1]) == [method, TestMethod()]
        # test error 
        @test_throws ErrorException set_derivative_method(gvrefs[2], TestGenMethod())
    end
    # test set_all_derivative_methods
    @testset "set_all_derivative_methods" begin 
        @infinite_parameter(m, pref in [0, 1])
        @test set_all_derivative_methods(m, TestMethod()) isa Nothing
        @test InfiniteOpt._derivative_methods(prefs[1]) == [TestMethod(), TestMethod()]
        @test derivative_method(pref) == TestMethod()
        @test_throws ErrorException set_all_derivative_methods(m, TestGenMethod())
    end
end

# test Infinite Domain Methods
@testset "Infinite Domain Methods" begin
    # setup data
    m = InfiniteModel();
    gvrefs1 = @infinite_parameter(m, a[1:2] in [0, 1], num_supports = 2)
    prefs1 = dispatch_variable_ref.(gvrefs1)
    gvrefs2 = @infinite_parameter(m, b[1:2, 1:2] ~ MatrixBeta(2, 2, 2))
    prefs2 = dispatch_variable_ref.(gvrefs2)
    push!(InfiniteOpt._constraint_dependencies(prefs1[1]), InfOptConstraintIndex(1))
    # test _parameter_domain (raw domain)
    @testset "_parameter_domain (Raw Domain)" begin
        @test InfiniteOpt._parameter_domain(prefs1[1]) isa CollectionDomain
        @test InfiniteOpt._parameter_domain(prefs2[1, 1]) isa MultiDistributionDomain
    end
    # test _parameter_domain (CollectionDomain)
    @testset "_parameter_domain (CollectionDomain)" begin
        domain = InfiniteOpt._parameter_domain(prefs1[1])
        @test InfiniteOpt._parameter_domain(domain, prefs1[1]) == IntervalDomain(0, 1)
    end
    # test _parameter_domain (Fallback)
    @testset "_parameter_domain (Fallback)" begin
        domain = InfiniteOpt._parameter_domain(prefs2[1])
        @test_throws ErrorException InfiniteOpt._parameter_domain(domain, prefs2[1])
    end
    # test infinite_domain (single)
    @testset "infinite_domain (Single)" begin
        @test infinite_domain(prefs1[1]) == IntervalDomain(0, 1)
        @test_throws ErrorException infinite_domain(prefs2[1])
        @test infinite_domain(gvrefs1[1]) == IntervalDomain(0, 1)
    end
    # test _check_complete_param_array
    @testset "_check_complete_param_array" begin
        @test InfiniteOpt._check_complete_param_array(prefs1) isa Nothing
        @test_throws ErrorException InfiniteOpt._check_complete_param_array(prefs2[:, 1])
    end
    # test infinite_domain (Array)
    @testset "infinite_domain (Array)" begin
        @test infinite_domain(prefs1) isa CollectionDomain
        @test_throws ErrorException infinite_domain(prefs2[:, 1])
        @test infinite_domain(prefs2) isa MultiDistributionDomain
        @test infinite_domain(gvrefs1) isa CollectionDomain
    end
    # test _update_parameter_domain
    @testset "_update_parameter_domain" begin
        old_domain = infinite_domain(prefs1)
        new_domain = CollectionDomain([IntervalDomain(0, 1), IntervalDomain(0, 2)])
        @test InfiniteOpt._update_parameter_domain(prefs1[1], new_domain) isa Nothing
        @test num_supports(prefs1) == 0
        @test !optimizer_model_ready(m)
        @test infinite_domain(prefs1) == new_domain
        @test InfiniteOpt._update_parameter_domain(prefs1[1], old_domain) isa Nothing
    end
    # test set_infinite_domain (Single)
    @testset "set_infinite_domain (Single)" begin
        # test normal
        @test set_infinite_domain(prefs1[2], IntervalDomain(-1, 1)) isa Nothing
        @test infinite_domain(prefs1[2]) == IntervalDomain(-1, 1)
        @test set_infinite_domain(gvrefs1[2], IntervalDomain(0, 1)) isa Nothing
        @test set_infinite_domain(prefs1[1], UniDistributionDomain(Normal())) isa Nothing
        @test infinite_domain(prefs1[1]) isa UniDistributionDomain
        @test set_infinite_domain(prefs1[1], IntervalDomain(0, 1)) isa Nothing
        # test errors
        @test_throws ErrorException set_infinite_domain(prefs2[1], IntervalDomain(0, 2))
        push!(InfiniteOpt._measure_dependencies(prefs1[1]), MeasureIndex(1))
        @test_throws ErrorException set_infinite_domain(prefs1[1], IntervalDomain(0, 2))
        empty!(InfiniteOpt._measure_dependencies(prefs1[1]))
    end
    # test set_infinite_domain (Array)
    @testset "set_infinite_domain (Array)" begin
        # test normal
        old_domain = infinite_domain(prefs1)
        new_domain = CollectionDomain([IntervalDomain(0, 1), IntervalDomain(0, 2)])
        @test set_infinite_domain(prefs1, new_domain) isa Nothing
        @test infinite_domain(prefs1) == new_domain
        @test set_infinite_domain(gvrefs1, old_domain) isa Nothing
        # test errors
        @test_throws ErrorException set_infinite_domain(prefs2[:, 1], new_domain)
        push!(InfiniteOpt._measure_dependencies(prefs1[1]), MeasureIndex(1))
        @test_throws ErrorException set_infinite_domain(prefs1, new_domain)
        empty!(InfiniteOpt._measure_dependencies(prefs1[1]))
    end
    # test JuMP.has_lower_bound
    @testset "JuMP.has_lower_bound" begin
        @test has_lower_bound(prefs1[1])
        @test has_lower_bound(gvrefs1[2])
        @test !has_lower_bound(prefs2[2])
    end
    # test JuMP.lower_bound
    @testset "JuMP.lower_bound" begin
        @test lower_bound(prefs1[1]) == 0
        @test lower_bound(gvrefs1[2]) == 0
        @test_throws ErrorException lower_bound(prefs2[2])
    end
    # test JuMP.set_lower_bound
    @testset "JuMP.set_lower_bound" begin
        @test set_lower_bound(prefs1[1], -1) isa Nothing
        @test lower_bound(prefs1[1]) == -1
        @test set_lower_bound(gvrefs1[2], 0) isa Nothing
        @test lower_bound(prefs1[2]) == 0
        @test_throws ErrorException set_lower_bound(prefs2[2], -3)
    end
    # test JuMP.has_upper_bound
    @testset "JuMP.has_upper_bound" begin
        @test has_upper_bound(prefs1[1])
        @test has_upper_bound(gvrefs1[2])
        @test !has_upper_bound(prefs2[2])
    end
    # test JuMP.upper_bound
    @testset "JuMP.upper_bound" begin
        @test upper_bound(prefs1[1]) == 1
        @test upper_bound(gvrefs1[2]) == 1
        @test_throws ErrorException upper_bound(prefs2[2])
    end
    # test JuMP.set_upper_bound
    @testset "JuMP.set_upper_bound" begin
        @test set_upper_bound(prefs1[1], 42) isa Nothing
        @test upper_bound(prefs1[1]) == 42
        @test set_upper_bound(gvrefs1[2], 0) isa Nothing
        @test upper_bound(prefs1[2]) == 0
        @test_throws ErrorException set_upper_bound(prefs2[2], 42)
    end
end

# test Support Methods
@testset "Support Methods" begin
    # Setup data
    m = InfiniteModel();
    gvrefs1 = @infinite_parameter(m, a[1:2] in [0, 1], num_supports = 2)
    prefs1 = dispatch_variable_ref.(gvrefs1)
    gvrefs2 = @infinite_parameter(m, b[1:2, 1:2] ~ MatrixBeta(2, 2, 2))
    prefs2 = dispatch_variable_ref.(gvrefs2)
    gvrefs3 = @infinite_parameter(m, c[1:2] in [0, 1])
    gvrefs4 = @infinite_parameter(m, d[1:2] in [0, 1], supports = 0)
    push!(InfiniteOpt._constraint_dependencies(prefs1[1]), InfOptConstraintIndex(1))
    # test has_internal_supports 
    @testset "has_internal_supports" begin
        @test !has_internal_supports(prefs1[1])
        @test !has_internal_supports(gvrefs4[2])
    end
    # test _set_has_internal_supports 
    @testset "_set_has_internal_supports" begin
        @test InfiniteOpt._set_has_internal_supports(prefs1[1], true) isa Nothing
        @test has_internal_supports(prefs1[1])
        @test InfiniteOpt._set_has_internal_supports(gvrefs1[1], false) isa Nothing
        @test !has_internal_supports(prefs1[1])
    end
    # test _parameter_supports
    @testset "_parameter_supports" begin
        @test sort(collect(keys(InfiniteOpt._parameter_supports(prefs1[1])))) == [[0., 0.], [1., 1.]]
        @test isempty(InfiniteOpt._parameter_supports(prefs2[2]))
    end
    # test significant_digits
    @testset "significant_digits" begin
        @test significant_digits(prefs1[1]) == InfiniteOpt.DefaultSigDigits
        @test significant_digits(gvrefs1[2]) == InfiniteOpt.DefaultSigDigits
    end
    # test num_supports (Single)
    @testset "num_supports (Single)" begin
        # test default
        @test num_supports(prefs1[1]) == 2
        @test num_supports(gvrefs1[2]) == 2
        @test num_supports(prefs2[2]) == 0
        # test label
        @test num_supports(prefs1[1], label = MCSample) == 0
        @test num_supports(prefs1[1], label = UniformGrid) == 2
        @test num_supports(gvrefs1[2], label = MCSample) == 0
        @test num_supports(prefs1[1], label = All) == 2
        @test num_supports(prefs1[1], label = SampleLabel) == 0
    end
    # test num_supports (Array)
    @testset "num_supports (Array)" begin
        # test default
        @test num_supports(prefs1) == 2
        @test num_supports(gvrefs1) == 2
        @test num_supports(prefs2) == 0
        # test label
        @test num_supports(prefs1, label = MCSample) == 0
        @test num_supports(prefs1, label = UniformGrid) == 2
        @test num_supports(gvrefs1, label = MCSample) == 0
        # test error
        @test_throws ErrorException num_supports(prefs2[:, 1])
    end
    # test has_supports (Single)
    @testset "has_supports (Single)" begin
        @test has_supports(prefs1[1])
        @test has_supports(gvrefs1[2])
        @test !has_supports(prefs2[2])
    end
    # test has_supports (Array)
    @testset "has_supports (Array)" begin
        @test has_supports(prefs1)
        @test has_supports(gvrefs1)
        @test !has_supports(prefs2)
        # test error
        @test_throws ErrorException has_supports(prefs2[:, 1])
    end
    # test supports (Single)
    @testset "supports (Single)" begin
        # test default
        @test sort!(supports(prefs1[1])) == Float64[0, 1]
        @test sort!(supports(gvrefs1[2])) == Float64[0, 1]
        @test supports(prefs2[2]) == Float64[]
        # test label
        @test supports(prefs1[1], label = MCSample) == Float64[]
        @test sort!(supports(prefs1[1], label = UniformGrid)) == Float64[0, 1]
        @test supports(gvrefs1[2], label = MCSample) == Float64[]
        @test supports(gvrefs4[1], label = All) == [0]
    end
    # test supports (AbstractArray)
    @testset "supports (AbstractArray)" begin
        # test default
        new_prefs1 = convert(JuMPC.SparseAxisArray, prefs1)
        new_gvrefs1 = convert(JuMPC.SparseAxisArray, gvrefs1)
        expected = convert.(JuMPC.SparseAxisArray, [zeros(2), ones(2)])
        @test sort!(supports(new_prefs1)) == expected
        @test sort!(supports(new_gvrefs1)) == expected
        @test supports(prefs2) == []
        # test label
        @test supports(new_prefs1, label = MCSample) == []
        @test sort!(supports(new_prefs1, label = UniformGrid)) == expected
        @test supports(new_gvrefs1, label = MCSample) == []
        @test supports(gvrefs4, label = All) == zeros(2, 1)
        @test supports(gvrefs4, label = AbstractSupportLabel) == zeros(2, 1)
    end
    # test supports (Vector)
    @testset "supports (Vector)" begin
        # test default
        @test sortcols(supports(prefs1)) == [0 1; 0 1]
        @test sortcols(supports(gvrefs1)) == [0 1; 0 1]
        @test supports(gvrefs3) == zeros(2, 0)
        @test supports(gvrefs4) == zeros(2, 1)
        # test label
        @test supports(prefs1, label = MCSample) == zeros(2, 0)
        @test sortcols(supports(prefs1, label = UniformGrid)) == [0 1; 0 1]
        @test supports(gvrefs1, label = MCSample) == zeros(2, 0)
        @test supports(gvrefs1, label = InternalLabel) == zeros(2, 0)
        # test error
        @test_throws ErrorException supports(prefs2[:, 1])
    end
    # test _update_parameter_supports
    @testset "_update_parameter_supports" begin
        old_supports = supports(prefs1)
        @test InfiniteOpt._update_parameter_supports(prefs1, ones(Int, 2, 3),
                                                     UserDefined) isa Nothing
        @test supports(prefs1) == ones(Float64, 2, 1)
        @test !optimizer_model_ready(m)
        @test InfiniteOpt._update_parameter_supports(prefs1, old_supports,
                                                     UniformGrid) isa Nothing
    end
    # test _make_support_matrix
    @testset "_make_support_matrix" begin
        # test errors
        supps = [ones(3, 2)]
        inds = IC.indices(prefs2)
        @test_throws ErrorException InfiniteOpt._make_support_matrix(supps, inds)
        supps = [ones(2, 2), ones(3, 2)]
        @test_throws ErrorException InfiniteOpt._make_support_matrix(supps, inds)
        # test normal
        supps = [zeros(2, 2), ones(2, 2)]
        @test InfiniteOpt._make_support_matrix(supps, inds) == hcat(zeros(4), ones(4))
    end
    # test set_supports (Single)
    @testset "set_supports (Single)" begin
        @test_throws ErrorException set_supports(prefs1[1], [0, 1])
        @test_throws ErrorException set_supports(gvrefs1[1], [0, 1], label = Mixture)
    end
    # test set_supports (Array)
    @testset "set_supports (Array)" begin
        # test errors
        @test_throws ErrorException set_supports(prefs2[:, 1], ones(2, 2))
        @test_throws ErrorException set_supports(prefs1, ones(2, 2))
        @test_throws ErrorException set_supports(prefs1, ones(2, 2) * 2, force = true)
        # test default
        old_supports = supports(prefs1)
        @test set_supports(prefs1, ones(2, 3), force = true) isa Nothing
        @test supports(prefs1) == ones(Float64, 2, 1)
        @test collect(values(InfiniteOpt._parameter_supports(prefs1[1]))) == [Set([UserDefined])]
        # test keywords
        @test set_supports(gvrefs1, old_supports, force = true,
                           label = UniformGrid) isa Nothing
        @test sortcols(supports(prefs1)) == sortcols(old_supports)
        @test collect(values(InfiniteOpt._parameter_supports(prefs1[1]))) == [Set([UniformGrid]) for i = 1:2]
        # test setting with internal supports 
        @test set_supports(gvrefs1, old_supports, force = true,
                            label = InternalLabel) isa Nothing
        @test sortcols(supports(prefs1, label = All)) == sortcols(old_supports)
        @test has_internal_supports(prefs1[1])
        @test set_supports(gvrefs1, old_supports, force = true,
                           label = UniformGrid) isa Nothing
        @test !has_internal_supports(prefs1[1])
    end
    # test set_supports (Vector of Supports)
    @testset "set_supports (Vector of Supports)" begin
        # default
        supps = [ones(2, 2) for i in 1:3]
        @test set_supports(prefs2, supps) isa Nothing
        @test supports(prefs2) == [ones(Float64, 2, 2)]
        @test collect(values(InfiniteOpt._parameter_supports(prefs2[1]))) == [Set([UserDefined])]
        # test keywords
        # supps = [Float64[] for i in CartesianIndices(prefs2)]
        # @test set_supports(gvrefs2, supps, force = true, label = Mixture) isa Nothing
        # @test supports(prefs2) == [ones(Float64, 0) for i in CartesianIndices(prefs2)]
        empty!(InfiniteOpt._parameter_supports(prefs2[1]))
    end
    # test add_supports (Single)
    @testset "add_supports (Single)" begin
        @test_throws ErrorException add_supports(prefs1[1], [0, 1])
        @test_throws ErrorException add_supports(gvrefs1[1], [0, 1], label = Mixture)
    end
    # test add_supports (Vector)
    @testset "add_supports (Vector)" begin
        # test errors
        @test_throws ErrorException add_supports(prefs2[:, 1], ones(2, 2))
        @test_throws ErrorException add_supports(prefs1, ones(2, 2) * 2)
        # test default
        @test add_supports(prefs1, 0.1 * ones(2, 2)) isa Nothing
        @test sortcols(supports(prefs1)) == Float64[0 0.1 1; 0 0.1 1]
        @test InfiniteOpt._parameter_supports(prefs1[1])[0.1 * ones(2)] == Set([UserDefined])
        # test keywords
        @test add_supports(gvrefs1, zeros(2, 1), check = false,
                           label = InternalLabel) isa Nothing
        @test sortcols(supports(prefs1)) == Float64[0 0.1 1; 0 0.1 1]
        @test InfiniteOpt._parameter_supports(prefs1[1])[zeros(2)] == Set([UniformGrid, InternalLabel])
        # test warning 
        InfiniteOpt._data_object(prefs1[1]).has_deriv_constrs[1] = true
        warn = "Support/method changes will invalidate existing derivative evaluation " *
               "constraints that have been added to the InfiniteModel. Thus, " *
               "these are being deleted."
        @test_logs (:warn, warn) add_supports(prefs1, 0.2 * ones(2, 2)) isa Nothing
        InfiniteOpt._data_object(prefs1[1]).has_deriv_constrs[1] = false
    end
    # test add_supports (AbstractArray)
    @testset "add_supports (AbstractArray)" begin
        # default
        supps = [ones(2, 2) for i = 1:3]
        @test add_supports(prefs2, supps) isa Nothing
        @test supports(prefs2) == [ones(Float64, 2, 2)]
        @test collect(values(InfiniteOpt._parameter_supports(prefs2[1]))) == [Set([UserDefined])]
        # test keywords
        supps = [ones(2, 2) * 0.5]
        @test add_supports(gvrefs2, supps, check = false, label = MCSample) isa Nothing
        @test sort!(supports(prefs2), by = first) == [0.5 * ones(2, 2), ones(2, 2)]
        @test Set([MCSample]) in values(InfiniteOpt._parameter_supports(prefs2[1]))
    end
    # test delete_supports (Single)
    @testset "delete_supports (Single)" begin
        @test_throws ErrorException delete_supports(prefs1[1])
        @test_throws ErrorException delete_supports(gvrefs1[1])
    end
    # test delete_supports (AbstractArray)
    @testset "delete_supports (AbstractArray)" begin
        # test errors
        @test_throws ErrorException delete_supports(prefs2[:, 1])
        push!(InfiniteOpt._measure_dependencies(prefs1[1]), MeasureIndex(1))
        @test_throws ErrorException delete_supports(prefs1)
        empty!(InfiniteOpt._measure_dependencies(prefs1[1]))
        # test label deletion 
        @test delete_supports(prefs1, label = InternalLabel) isa Nothing
        @test !has_internal_supports(prefs1[1])
        # normal
        @test delete_supports(prefs1) isa Nothing
        @test supports(prefs1) == zeros(Float64, 2, 0)
        @test delete_supports(gvrefs2) isa Nothing
        @test supports(prefs2) == []
    end
end

# test Support Filling
@testset "Support Filling" begin
    # Setup data
    m = InfiniteModel();
    gvrefs1 = @infinite_parameter(m, a[1:2] in [0, 1], num_supports = 2)
    prefs1 = dispatch_variable_ref.(gvrefs1)
    gvrefs2 = @infinite_parameter(m, b[1:2, 1:2] ~ MatrixBeta(2, 2, 2))
    prefs2 = dispatch_variable_ref.(gvrefs2)
    gvref = @infinite_parameter(m, c in [0, 1])
    pref = dispatch_variable_ref(gvref)
    # test generate_and_add_supports!
    @testset "generate_and_add_supports!" begin
        # old_supports = supports(prefs1)
        domain = infinite_domain(prefs1)
        @test generate_and_add_supports!(prefs1, domain, num_supports = 2) isa Nothing
        @test sortcols(supports(prefs1)) == Float64[0 1; 0 1]
        expected = [Set([UniformGrid]) for i = 1:2]
        @test collect(values(InfiniteOpt._parameter_supports(prefs1[1]))) == expected
        # @test set_supports(prefs1, old_supports, force = true) isa Nothing
    end
    # test fill_in_supports! (Single)
    @testset "fill_in_supports! (Single)" begin
        @test_throws ErrorException fill_in_supports!(prefs1[1])
        @test_throws ErrorException fill_in_supports!(gvrefs1[1])
    end
    # test fill_in_supports! (AbstractArray)
    @testset "fill_in_supports! (AbstractArray)" begin
        # test error
        @test_throws ErrorException fill_in_supports!(prefs2[1, :])
        # test default
        @test fill_in_supports!(prefs1) isa Nothing
        @test num_supports(prefs1) == 8
        @test fill_in_supports!(gvrefs2) isa Nothing
        @test num_supports(prefs2) == 10
        # test keywords
        @test fill_in_supports!(prefs1, num_supports = 2) isa Nothing
        @test num_supports(prefs1) == 8
        @test fill_in_supports!(gvrefs2, modify = false,
                                num_supports = 20) isa Nothing
        @test num_supports(prefs2) == 10
        # delete additions
        @test delete_supports(prefs1) isa Nothing
        @test delete_supports(prefs2) isa Nothing
    end
    # test fill_in_supports! (InfiniteModel)
    @testset "fill_in_supports! (InfiniteModel)" begin
        # test default
        @test fill_in_supports!(m) isa Nothing
        @test num_supports(prefs1) == 10
        @test num_supports(prefs2) == 10
        @test num_supports(pref) == 10
        # test with keywords
        @test fill_in_supports!(m, num_supports = 12, modify = false) isa Nothing
        @test num_supports(prefs1) == 10
        @test num_supports(prefs2) == 10
        @test num_supports(pref) == 10
    end
end

# Test parameter counting and listing methods
@testset "General Queries" begin
    # Setup data
    m = InfiniteModel();
    prefs1 = @infinite_parameter(m, a[1:2] in [0, 1], num_supports = 2)
    prefs2 = @infinite_parameter(m, b[1:2, 1:2] ~ MatrixBeta(2, 2, 2))
    pref = @infinite_parameter(m, c in [0, 1])
    fpref = @finite_parameter(m, d == 10)
    # test num_parameters (Default)
    @testset "num_parameters (Default)" begin
        @test num_parameters(m) == 8
    end
    # test num_parameters (Specific Scalar)
    @testset "num_parameters (Specific Scalar)" begin
        @test num_parameters(m, IndependentParameter) == 1
        @test num_parameters(m, FiniteParameter) == 1
    end
    # test num_parameters (ScalarParameter)
    @testset "num_parameters (ScalarParameter)" begin
        @test num_parameters(m, ScalarParameter) == 2
    end
    # test num_parameters (DependentParameters)
    @testset "num_parameters (DependentParameters)" begin
        @test num_parameters(m, DependentParameters) == 6
    end
    # test num_parameters (InfiniteParameter)
    @testset "num_parameters (InfiniteParameter)" begin
        @test num_parameters(m, InfiniteParameter) == 7
    end
    # test all_parameters (Default)
    @testset "all_parameters (Default)" begin
        @test isequal(all_parameters(m), [pref; prefs1; [prefs2...]; fpref])
    end
    # test all_parameters (Specific Scalar)
    @testset "all_parameters (Specific Scalar)" begin
        @test isequal(all_parameters(m, IndependentParameter), [pref])
        @test isequal(all_parameters(m, FiniteParameter), [fpref])
    end
    # test all_parameters (ScalarParameter)
    @testset "all_parameters (ScalarParameter)" begin
        @test isequal(all_parameters(m, ScalarParameter), [pref, fpref])
    end
    # test all_parameters (DependentParameters)
    @testset "all_parameters (DependentParameters)" begin
        @test isequal(all_parameters(m, DependentParameters), [prefs1; [prefs2...]])
    end
    # test all_parameters (InfiniteParameter)
    @testset "all_parameters (InfiniteParameter)" begin
        @test isequal(all_parameters(m, InfiniteParameter), [prefs1; [prefs2...]; pref])
    end
end
