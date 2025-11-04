# test Base Extensions
@testset "Base Extensions" begin
    # Setup data
    m = InfiniteModel();
    gvref = GeneralVariableRef(m, 1, PointVariableIndex)
    gvref2 = GeneralVariableRef(m, 2, PointVariableIndex)
    ivref = InfiniteVariableRef(m, InfiniteVariableIndex(1))
    hvref = FiniteVariableRef(m, FiniteVariableIndex(1))
    # test Base.copy
    @testset "Base.copy" begin
        @test isequal(copy(gvref), gvref)
        @test isequal(copy(ivref), ivref)
    end
    # test Base.isequal
    @testset "Base.isequal" begin
        @test isequal(gvref, gvref)
        @test isequal(ivref, ivref)
        @test !isequal(gvref, ivref)
        @test !isequal(hvref, ivref)
        @test !isequal(gvref, gvref2)
    end
    # test Base.broadcastable
    @testset "Base.broadcastable" begin
        @test Base.broadcastable(gvref) isa Base.RefValue{GeneralVariableRef}
        @test Base.broadcastable(ivref) isa Base.RefValue{InfiniteVariableRef}
    end
    # test Base.length
    @testset "Base.length" begin
        @test length(gvref) == 1
        @test length(ivref) == 1
    end
    # test JuMP.isequal_canonical
    @testset "JuMP.isequal_canonical" begin
        @test JuMP.isequal_canonical(gvref, gvref)
        @test JuMP.isequal_canonical(ivref, ivref)
    end
    # test _remove_name_index
    @testset "_remove_name_index" begin
        pref = @infinite_parameter(m, test, domain = IntervalDomain(0, 1))
        pref2 = @infinite_parameter(m, θ, domain = IntervalDomain(0, 1))
        prefs = @infinite_parameter(m, x[1:2], domain = IntervalDomain(0, 1))
        @test InfiniteOpt._remove_name_index(pref) == "test"
        @test InfiniteOpt._remove_name_index(prefs[1]) == "x"
        @test InfiniteOpt._remove_name_index(pref2) == "θ"
    end
    # test Basic getters
    @testset "Attribute Accessors" begin
        @test InfiniteOpt._index_type(gvref) == PointVariableIndex
        @test InfiniteOpt._raw_index(gvref) == 1
        @test InfiniteOpt._param_index(gvref) == -1
    end
end

# test Reference Accessers
@testset "Reference Accessers" begin
    # Setup data
    m = InfiniteModel();
    gvref1 = GeneralVariableRef(m, 1, PointVariableIndex)
    gvref2 = GeneralVariableRef(m, 1, DependentParameterIndex, 2)
    ivref = InfiniteVariableRef(m, InfiniteVariableIndex(1))
    # test JuMP.index (GeneralVariableRef)
    @testset "JuMP.index (GeneralVariableRef)" begin
        @test index(gvref1) == PointVariableIndex(1)
        @test index(gvref2) == DependentParameterIndex(DependentParametersIndex(1), 2)
    end
    # test JuMP.index (DispatchVariableRef)
    @testset "JuMP.index (DispatchVariableRef)" begin
        @test index(ivref) == InfiniteVariableIndex(1)
    end
    # test JuMP.owner_model (GeneralVariableRef)
    @testset "JuMP.owner_model (GeneralVariableRef)" begin
        @test owner_model(gvref1) == m
        @test owner_model(gvref2) == m
    end
    # test JuMP.owner_model (DispatchVariableRef)
    @testset "JuMP.owner_model (DispatchVariableRef)" begin
        @test owner_model(ivref) == m
    end
end

# test Dispatch Variabel Methods
@testset "dispatch_variable_ref" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test the test type
    @testset "TestVariableRef" begin
        @test isequal(dispatch_variable_ref(m, idx), dvref)
    end
    # test GeneralVariableRef input
    @testset "GeneralVariableRef" begin
        @test isequal(dispatch_variable_ref(gvref), dvref)
    end
end

# test Naming Methods
@testset "Naming Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test JuMP.name (Fallback)
    @testset "JuMP.name (Fallback)" begin
        @test_throws ArgumentError name(dvref)
    end
    # test JuMP.name (GeneralVariableRef)
    @testset "JuMP.name (GeneralVariableRef)" begin
        @test_throws ArgumentError name(gvref)
    end
    # test JuMP.set_name (Fallback)
    @testset "JuMP.set_name (Fallback)" begin
        @test_throws ArgumentError set_name(dvref, "new")
    end
    # test JuMP.set_name (GeneralVariableRef)
    @testset "JuMP.set_name (GeneralVariableRef)" begin
        @test_throws ArgumentError set_name(gvref, "new")
    end
end

# test Core data methods
@testset "Core Data Methods" begin
    # Setup data
    m = InfiniteModel();
    gvref = GeneralVariableRef(m, 1, TestIndex)
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    # test _data_dictionary with GeneralVariableRef
    @testset "_data_dictionary (GeneralVariableRef)" begin
        @test_throws MethodError InfiniteOpt._data_dictionary(gvref)
    end
    # test _data_object with GeneralVariableRef
    @testset "_data_object (GeneralVariableRef)" begin
        @test_throws MethodError InfiniteOpt._data_object(gvref)
    end
    # test _delete_data_object with DispatchVariableRef
    @testset "_delete_data_object (DispatchVariableRef)" begin
        @test_throws MethodError InfiniteOpt._delete_data_object(dvref)
    end
end

# test Validity Methods
@testset "Validity Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test JuMP.is_valid (Fallback)
    @testset "JuMP.is_valid (Fallback)" begin
        @test_throws MethodError is_valid(m, dvref)
    end
    # test JuMP.is_valid (GeneralVariableRef)
    @testset "JuMP.is_valid (GeneralVariableRef)" begin
        @test_throws MethodError is_valid(m, gvref)
    end
end

# test Core Object Methods
@testset "Core Object Methods" begin
    # Setup data
    m = InfiniteModel();
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test _core_variable_object (GeneralVariableRef)
    @testset "_core_variable_object (GeneralVariableRef)" begin
        @test_throws MethodError core_object(gvref)
    end
end

# test Dependency Methods
@testset "Dependency Methods" begin
    # Setup data
    m = InfiniteModel();
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # Loop over the various functions
    for f in (:_infinite_variable_dependencies, :_semi_infinite_variable_dependencies,
              :_point_variable_dependencies, :_measure_dependencies,
              :_constraint_dependencies, :_derivative_dependencies, 
              :_derivative_constraint_dependencies, 
              :_parameter_function_dependencies, :_generative_measures)
        @test_throws MethodError InfiniteOpt.eval(f)(gvref)
    end
end

# test Usage Methods
@testset "Usage Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    for f in (:used_by_infinite_variable, :used_by_semi_infinite_variable,
              :used_by_point_variable, :used_by_measure, :used_by_constraint,
              :used_by_objective, :used_by_derivative, :is_used, 
              :has_derivative_constraints, :used_by_parameter_function)
        @test_throws ArgumentError eval(f)(dvref)
        @test_throws ArgumentError eval(f)(gvref)
    end
end

# test Delete Methods
@testset "Delete Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test JuMP.delete (Fallback)
    @testset "JuMP.delete (Fallback)" begin
        @test_throws ArgumentError JuMP.delete(m, dvref)
    end
    # test JuMP.delete (GeneralVariableRef)
    @testset "JuMP.delete (GeneralVariableRef)" begin
        @test_throws ArgumentError JuMP.delete(m, gvref)
    end
    # test JuMP.delete (GeneralVariableRefs)
    @testset "JuMP.delete (GeneralVariableRefs)" begin
        @test_throws ArgumentError JuMP.delete(m, [gvref])
    end
end

# test calling the general variable reference as a function
@testset "Functional Calls" begin
    # Setup data
    m = InfiniteModel();
    gvref = GeneralVariableRef(m, 1, TestIndex2)
    # test _functional_reference_call (Fallback)
    @testset "_functional_reference_call (Fallback)" begin
        @test_throws ErrorException InfiniteOpt._functional_reference_call(gvref, TestIndex2, 42)
    end
    # test GeneralVariableRef as a function
    @testset "GeneralVariableRef(args...)" begin
        @test_throws ErrorException gvref(42, a = 4)
    end
end

# test Parameter Methods
@testset "Parameter Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    gvrefs = [GeneralVariableRef(m, 1, TestIndex),
              GeneralVariableRef(m, 2, TestIndex)]
    sdomain = IntervalDomain(0, 1)
    adomain = CollectionDomain([sdomain, sdomain])
    # test parameter_group_int_index
    @testset "parameter_group_int_index" begin
        @test_throws MethodError InfiniteOpt.parameter_group_int_index(gvref)
    end
    # test 1 argument methods and fallbacks 
    for f in (:significant_digits, :parameter_value, :derivative_method, 
              :has_generative_supports, :has_internal_supports,
              :add_generative_supports, :raw_function, :generative_support_info)
        @test_throws ArgumentError eval(f)(dvref)
        @test_throws ArgumentError eval(f)(gvref)
    end
    # test infinite_domain (Fallback)
    @testset "infinite_domain (Fallback)" begin
        @test_throws ArgumentError infinite_domain(dvref)
    end
    # test infinite_domain (GeneralVariableRef)
    @testset "infinite_domain (GeneralVariableRef)" begin
        @test_throws ArgumentError infinite_domain(gvref)
    end
    # test infinite_domain (GeneralVariableRef Array)
    @testset "infinite_domain (GeneralVariableRef Array)" begin
        @test_throws ArgumentError infinite_domain(gvrefs)
    end
    # test set_infinite_domain (Fallback)
    @testset "set_infinite_domain (Fallback)" begin
        @test_throws ArgumentError set_infinite_domain(dvref, sdomain)
    end
    # test set_infinite_domain (GeneralVariableRef)
    @testset "set_infinite_domain (GeneralVariableRef)" begin
        @test_throws ArgumentError set_infinite_domain(gvref, sdomain)
    end
    # test set_infinite_domain (GeneralVariableRef Array)
    @testset "set_infinite_domain (GeneralVariableRef Array)" begin
        @test_throws ArgumentError set_infinite_domain(gvrefs, adomain)
    end
    # test num_supports (Fallback)
    @testset "num_supports (Fallback)" begin
        @test_throws ArgumentError num_supports(dvref, label = All)
    end
    # test num_supports (GeneralVariableRef)
    @testset "num_supports (GeneralVariableRef)" begin
        @test_throws ArgumentError num_supports(gvref, label = All)
    end
    # test num_supports (GeneralVariableRef Array)
    @testset "num_supports (GeneralVariableRef Array)" begin
        @test_throws ArgumentError num_supports(gvrefs, label = All)
    end
    # test has_supports (Fallback)
    @testset "has_supports (Fallback)" begin
        @test_throws ArgumentError has_supports(dvref)
    end
    # test has_supports (GeneralVariableRef)
    @testset "has_supports (GeneralVariableRef)" begin
        @test_throws ArgumentError has_supports(gvref)
    end
    # test has_supports (GeneralVariableRef Array)
    @testset "has_supports (GeneralVariableRef Array)" begin
        @test_throws ArgumentError has_supports(gvrefs)
    end
    # test supports (Fallback)
    @testset "supports (Fallback)" begin
        @test_throws ArgumentError supports(dvref, label = All)
    end
    # test supports (GeneralVariableRef)
    @testset "supports (GeneralVariableRef)" begin
        @test_throws ArgumentError supports(gvref, label = All)
    end
    # test supports (GeneralVariableRef Array)
    @testset "supports (GeneralVariableRef Array)" begin
        @test_throws ArgumentError supports(gvrefs, label = All)
    end
    # test set_supports (Fallback)
    @testset "set_supports (Fallback)" begin
        @test_throws ArgumentError set_supports(dvref, 2, label = All, force = true)
    end
    # test set_supports (GeneralVariableRef)
    @testset "set_supports (GeneralVariableRef)" begin
        @test_throws ArgumentError set_supports(gvref, 2, label = All, force = true)
    end
    # test set_supports (GeneralVariableRef Array)
    @testset "set_supports (GeneralVariableRef Array)" begin
        @test_throws ArgumentError set_supports(gvrefs, ones(1, 1), label = All,
                                                force = true)
    end
    # test add_supports (Fallback)
    @testset "add_supports (Fallback)" begin
        @test_throws ArgumentError add_supports(dvref, 2, label = All, check = true)
    end
    # test add_supports (GeneralVariableRef)
    @testset "add_supports (GeneralVariableRef)" begin
        @test_throws ArgumentError add_supports(gvref, 2, label = All, check = true)
    end
    # test add_supports (GeneralVariableRef Array)
    @testset "add_supports (GeneralVariableRef Array)" begin
        @test_throws ArgumentError add_supports(gvrefs, ones(1, 1), label = All,
                                                check = true)
    end
    # test delete_supports (Fallback)
    @testset "delete_supports (Fallback)" begin
        @test_throws ArgumentError delete_supports(dvref)
    end
    # test delete_supports (GeneralVariableRef)
    @testset "delete_supports (GeneralVariableRef)" begin
        @test_throws ArgumentError delete_supports(gvref)
    end
    # test delete_supports (GeneralVariableRef Array)
    @testset "delete_supports (GeneralVariableRef Array)" begin
        @test_throws ArgumentError delete_supports(gvrefs)
    end
    # test fill_in_supports! (Fallback)
    @testset "fill_in_supports! (Fallback)" begin
        @test_throws ArgumentError fill_in_supports!(dvref, num_supports = 1)
    end
    # test fill_in_supports! (GeneralVariableRef)
    @testset "fill_in_supports! (GeneralVariableRef)" begin
        @test_throws ArgumentError fill_in_supports!(gvref, num_supports = 1)
    end
    # test fill_in_supports! (GeneralVariableRef Array)
    @testset "fill_in_supports! (GeneralVariableRef Array)" begin
        @test_throws ArgumentError fill_in_supports!(gvrefs, num_supports = 1)
    end
    # test JuMP.set_parameter_value (Fallback)
    @testset "JuMP.set_value (Fallback)" begin
        @test_throws ArgumentError set_parameter_value(dvref, 2)
    end
    # test JuMP.set_parameter_value (GeneralVariableRef)
    @testset "JuMP.set_value (GeneralVariableRef)" begin
        @test_throws ArgumentError set_parameter_value(gvref, 2)
    end
    # test set_derivative_method (Fallback)
    @testset "set_derivative_method (Fallback)" begin
        @test_throws ArgumentError set_derivative_method(dvref, FiniteDifference())
    end
    # test set_derivative_method (GeneralVariableRef)
    @testset "set_derivative_method (GeneralVariableRef)" begin
        @test_throws ArgumentError set_derivative_method(gvref, FiniteDifference())
    end
    # test call_function (Fallback)
    @testset "call_function (Fallback)" begin
        @test_throws ArgumentError call_function(dvref, 1)
    end
    # test call_function (GeneralVariableRef)
    @testset "call_function (GeneralVariableRef)" begin
        @test_throws ArgumentError call_function(gvref, 1)
    end
    # test parameter status setters
    for f in (:_set_has_generative_supports, :_set_has_internal_supports,
              :_set_has_derivative_constraints)
        @test_throws ArgumentError InfiniteOpt.eval(f)(dvref, true)
        @test_throws ArgumentError InfiniteOpt.eval(f)(gvref, true)
    end
end

# test Variable Methods
@testset "Variable Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test 1 argument methods 
    for f in (:raw_parameter_refs, :parameter_refs, :parameter_list,
              :infinite_variable_ref, :eval_support, :raw_parameter_values,
              :parameter_values)
        @test_throws ArgumentError eval(f)(dvref)
        @test_throws ArgumentError eval(f)(gvref)
    end
end

# test Measure Methods
@testset "Measure Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # loop of the methods 
    for f in (:measure_function, :measure_data, :is_analytic, :expand)
        @test_throws ArgumentError eval(f)(dvref)
        @test_throws ArgumentError eval(f)(gvref)
    end
end

# test Derivative Methods
@testset "Derivative Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # loop through methods 
    for f in (:derivative_argument, :operator_parameter, :evaluate, 
              :derivative_constraints, :delete_derivative_constraints,
              :derivative_order)
        @test_throws ArgumentError eval(f)(dvref)
        @test_throws ArgumentError eval(f)(gvref)
    end
end

# test variable info methods 
@testset "Variable Info Methods" begin 
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    pref = GeneralVariableRef(m, -1, IndependentParameterIndex)
    # test 1 argument methods 
    for f in (:has_lower_bound, :has_upper_bound, :lower_bound, :upper_bound, 
              :fix_value, :start_value, :set_binary, :set_integer, 
              :LowerBoundRef, :UpperBoundRef, :FixRef, :BinaryRef, :IntegerRef, 
              :delete_lower_bound, :delete_upper_bound, 
              :unfix, :unset_binary, :unset_integer)
        @test_throws ArgumentError eval(f)(dvref)
        @test_throws ArgumentError eval(f)(gvref)
    end
    for f in (:is_fixed, :is_binary, :is_integer)
        @test !eval(f)(dvref)
        @test !eval(f)(gvref)
    end
    # test setting methods
    for f in (:set_lower_bound, :set_upper_bound, :set_start_value)
        @test_throws ArgumentError eval(f)(dvref, 42)
        @test_throws ArgumentError eval(f)(gvref, 42)
    end
    # test JuMP.fix (Fallback)
    @testset "JuMP.fix (Fallback)" begin
        @test_throws ArgumentError fix(dvref, 42, force = true)
    end
    # test JuMP.fix (GeneralVariableRef)
    @testset "JuMP.fix (GeneralVariableRef)" begin
        @test_throws ArgumentError fix(gvref, 42, force = true)
    end
    # test constant_over_collocation (Fallback)
    @testset "constant_over_collocation (Fallback)" begin
        @test_throws ArgumentError constant_over_collocation(dvref, pref)
    end
    # test constant_over_collocation (GeneralVariableRef)
    @testset "constant_over_collocation (GeneralVariableRef)" begin
        @test_throws ArgumentError constant_over_collocation(gvref, pref)
    end
end
