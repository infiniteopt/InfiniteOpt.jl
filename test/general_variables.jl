# test Base Extensions
@testset "Base Extensions" begin
    # Setup data
    m = InfiniteModel();
    gvref = GeneralVariableRef(m, 1, PointVariableIndex)
    gvref2 = GeneralVariableRef(m, 2, PointVariableIndex)
    ivref = InfiniteVariableRef(m, InfiniteVariableIndex(1))
    hvref = HoldVariableRef(m, HoldVariableIndex(1))
    # test Base.copy
    @testset "Base.copy" begin
        @test copy(gvref) == gvref
        @test copy(ivref) == ivref
    end
    # test Base.:(==)
    @testset "Base.:(==)" begin
        @test gvref == gvref
        @test ivref == ivref
        @test !(gvref == ivref)
        @test !(hvref == ivref)
        @test !(gvref == gvref2)
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
    # test JuMP.variable_type
    @testset "JuMP.variable_type" begin
        @test JuMP.variable_type(m) == GeneralVariableRef
    end
    # test _remove_name_index
    @testset "_remove_name_index" begin
        pref = @independent_parameter(m, test, set = IntervalSet(0, 1))
        pref2 = @independent_parameter(m, θ, set = IntervalSet(0, 1))
        prefs = @independent_parameter(m, x[1:2], set = IntervalSet(0, 1))
        @test InfiniteOpt._remove_name_index(pref) == "test"
        @test InfiniteOpt._remove_name_index(prefs[1]) == "x"
        @test InfiniteOpt._remove_name_index(pref2) == "θ"
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
        @test dispatch_variable_ref(m, idx) == dvref
    end
    # test GeneralVariableRef input
    @testset "GeneralVariableRef" begin
        @test dispatch_variable_ref(gvref) == dvref
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
    idx = DependentParameterIndex(DependentParametersIndex(1), 1)
    pref = DependentParameterRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    gvref2 = GeneralVariableRef(m, 1, DependentParameterIndex, 1)
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
    # test _delete_data_object with DependentParameterRef
    @testset "_delete_data_object (DependentParameterRef)" begin
        @test_throws KeyError InfiniteOpt._delete_data_object(pref)
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
        @test_throws MethodError InfiniteOpt._core_variable_object(gvref)
    end
end

# test Dependency Methods
@testset "Dependency Methods" begin
    # Setup data
    m = InfiniteModel();
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test _infinite_variable_dependencies (GeneralVariableRef)
    @testset "_infinite_variable_dependencies (GeneralVariableRef)" begin
        @test_throws MethodError InfiniteOpt._infinite_variable_dependencies(gvref)
    end
    # test _reduced_variable_dependencies (GeneralVariableRef)
    @testset "_reduced_variable_dependencies (GeneralVariableRef)" begin
        @test_throws MethodError InfiniteOpt._reduced_variable_dependencies(gvref)
    end
    # test _point_variable_dependencies (GeneralVariableRef)
    @testset "_point_variable_dependencies (GeneralVariableRef)" begin
        @test_throws MethodError InfiniteOpt._point_variable_dependencies(gvref)
    end
    # test _measure_dependencies (GeneralVariableRef)
    @testset "_measure_dependenciess (GeneralVariableRef)" begin
        @test_throws MethodError InfiniteOpt._measure_dependencies(gvref)
    end
    # test _constraint_dependencies (GeneralVariableRef)
    @testset "_constraint_dependencies (GeneralVariableRef)" begin
        @test_throws MethodError InfiniteOpt._constraint_dependencies(gvref)
    end
end

# test Usage Methods
@testset "Usage Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test used_by_infinite_variable (Fallback)
    @testset "used_by_infinite_variable (Fallback)" begin
        @test_throws ArgumentError used_by_infinite_variable(dvref)
    end
    # test used_by_infinite_variable (GeneralVariableRef)
    @testset "used_by_infinite_variable (GeneralVariableRef)" begin
        @test_throws ArgumentError used_by_infinite_variable(gvref)
    end
    # test used_by_reduced_variable (Fallback)
    @testset "used_by_reduced_variable (Fallback)" begin
        @test_throws ArgumentError used_by_reduced_variable(dvref)
    end
    # test used_by_reduced_variable (GeneralVariableRef)
    @testset "used_by_reduced_variable (GeneralVariableRef)" begin
        @test_throws ArgumentError used_by_reduced_variable(gvref)
    end
    # test used_by_point_variable (Fallback)
    @testset "used_by_point_variable (Fallback)" begin
        @test_throws ArgumentError used_by_point_variable(dvref)
    end
    # test used_by_point_variable (GeneralVariableRef)
    @testset "used_by_point_variable (GeneralVariableRef)" begin
        @test_throws ArgumentError used_by_point_variable(gvref)
    end
    # test used_by_measure (Fallback)
    @testset "used_by_measure (Fallback)" begin
        @test_throws ArgumentError used_by_measure(dvref)
    end
    # test used_by_measure (GeneralVariableRef)
    @testset "used_by_measure (GeneralVariableRef)" begin
        @test_throws ArgumentError used_by_measure(gvref)
    end
    # test used_by_objective (Fallback)
    @testset "used_by_objective (Fallback)" begin
        @test_throws ArgumentError used_by_objective(dvref)
    end
    # test used_by_objective (GeneralVariableRef)
    @testset "used_by_objective (GeneralVariableRef)" begin
        @test_throws ArgumentError used_by_objective(gvref)
    end
    # test used_by_constraint (Fallback)
    @testset "used_by_constraint (Fallback)" begin
        @test_throws ArgumentError used_by_constraint(dvref)
    end
    # test used_by_constraint (GeneralVariableRef)
    @testset "used_by_constraint (GeneralVariableRef)" begin
        @test_throws ArgumentError used_by_constraint(gvref)
    end
    # test is_used (Fallback)
    @testset "is_used (Fallback)" begin
        @test_throws ArgumentError is_used(dvref)
    end
    # test is_used (GeneralVariableRef)
    @testset "is_used (GeneralVariableRef)" begin
        @test_throws ArgumentError is_used(gvref)
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
    sparse_gvrefs = convert(JuMPC.SparseAxisArray, gvrefs)
    dense_gvrefs = JuMPC.DenseAxisArray(gvrefs, 3:4)
    sset = IntervalSet(0, 1)
    aset = CollectionSet([sset, sset])
    # test _parameter_number
    @testset "_parameter_number" begin
        @test_throws MethodError InfiniteOpt._parameter_number(gvref)
    end
    # test _object_number
    @testset "_object_number" begin
        @test_throws MethodError InfiniteOpt._object_number(gvref)
    end
    # test infinite_set (Fallback)
    @testset "infinite_set (Fallback)" begin
        @test_throws ArgumentError infinite_set(dvref)
    end
    # test infinite_set (GeneralVariableRef)
    @testset "infinite_set (GeneralVariableRef)" begin
        @test_throws ArgumentError infinite_set(gvref)
    end
    # test infinite_set (GeneralVariableRef Array)
    @testset "infinite_set (GeneralVariableRef Array)" begin
        @test_throws ArgumentError infinite_set(gvrefs)
        @test_throws ArgumentError infinite_set(dense_gvrefs)
        @test_throws ArgumentError infinite_set(sparse_gvrefs)
    end
    # test set_infinite_set (Fallback)
    @testset "set_infinite_set (Fallback)" begin
        @test_throws ArgumentError set_infinite_set(dvref, sset)
    end
    # test set_infinite_set (GeneralVariableRef)
    @testset "set_infinite_set (GeneralVariableRef)" begin
        @test_throws ArgumentError set_infinite_set(gvref, sset)
    end
    # test set_infinite_set (GeneralVariableRef Array)
    @testset "set_infinite_set (GeneralVariableRef Array)" begin
        @test_throws ArgumentError set_infinite_set(gvrefs, aset)
        @test_throws ArgumentError set_infinite_set(dense_gvrefs, aset)
        @test_throws ArgumentError set_infinite_set(sparse_gvrefs, aset)
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
        @test_throws ArgumentError num_supports(dense_gvrefs, label = All)
        @test_throws ArgumentError num_supports(sparse_gvrefs, label = All)
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
        @test_throws ArgumentError has_supports(dense_gvrefs)
        @test_throws ArgumentError has_supports(sparse_gvrefs)
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
        @test_throws ArgumentError supports(dense_gvrefs, label = All)
        @test_throws ArgumentError supports(sparse_gvrefs, label = All)
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
        @test_throws ArgumentError set_supports(dense_gvrefs, ones(1, 1),
                                                label = All, force = true)
        @test_throws ArgumentError set_supports(sparse_gvrefs, ones(1, 1),
                                                label = All, force = true)
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
        @test_throws ArgumentError add_supports(dense_gvrefs, ones(1, 1),
                                                label = All, check = true)
        @test_throws ArgumentError add_supports(sparse_gvrefs, ones(1, 1),
                                                label = All, check = true)
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
        @test_throws ArgumentError delete_supports(dense_gvrefs)
        @test_throws ArgumentError delete_supports(sparse_gvrefs)
    end
    # test fill_in_supports! (Fallback)
    @testset "fill_in_supports! (Fallback)" begin
        @test_throws ArgumentError fill_in_supports!(dvref, num_supports = 1,
                                                     sig_figs = 1)
    end
    # test fill_in_supports! (GeneralVariableRef)
    @testset "fill_in_supports! (GeneralVariableRef)" begin
        @test_throws ArgumentError fill_in_supports!(gvref, num_supports = 1,
                                                     sig_figs = 1)
    end
    # test fill_in_supports! (GeneralVariableRef Array)
    @testset "fill_in_supports! (GeneralVariableRef Array)" begin
        @test_throws ArgumentError fill_in_supports!(gvrefs, num_supports = 1,
                                                     sig_figs = 1)
        @test_throws ArgumentError fill_in_supports!(dense_gvrefs,
                                                 num_supports = 1, sig_figs = 1)
        @test_throws ArgumentError fill_in_supports!(sparse_gvrefs,
                                                 num_supports = 1, sig_figs = 1)
    end
    # TODO insert more as needed
end

# test Variable Methods
@testset "Variable Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test ...
    # TODO insert as needed
end

# test Lower Bound Methods
@testset "Lower Bound Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test JuMP.has_lower_bound (Fallback)
    @testset "JuMP.has_lower_bound (Fallback)" begin
        @test_throws ArgumentError has_lower_bound(dvref)
    end
    # test JuMP.has_lower_bound (GeneralVariableRef)
    @testset "JuMP.has_lower_bound (GeneralVariableRef)" begin
        @test_throws ArgumentError has_lower_bound(gvref)
    end
    # test JuMP.lower_bound (Fallback)
    @testset "JuMP.lower_bound (Fallback)" begin
        @test_throws ArgumentError lower_bound(dvref)
    end
    # test JuMP.lower_bound (GeneralVariableRef)
    @testset "JuMP.lower_bound (GeneralVariableRef)" begin
        @test_throws ArgumentError lower_bound(gvref)
    end
    # test JuMP.set_lower_bound (Fallback)
    @testset "JuMP.set_lower_bound (Fallback)" begin
        @test_throws ArgumentError set_lower_bound(dvref, 42)
    end
    # test JuMP.set_lower_bound (GeneralVariableRef)
    @testset "JuMP.set_lower_bound (GeneralVariableRef)" begin
        @test_throws ArgumentError set_lower_bound(gvref, 42)
    end
    # test JuMP.LowerBoundRef (Fallback)
    @testset "JuMP.LowerBoundRef (Fallback)" begin
        @test_throws ArgumentError LowerBoundRef(dvref)
    end
    # test JuMP.LowerBoundRef (GeneralVariableRef)
    @testset "JuMP.LowerBoundRef (GeneralVariableRef)" begin
        @test_throws ArgumentError LowerBoundRef(gvref)
    end
    # test JuMP.delete_lower_bound (Fallback)
    @testset "JuMP.delete_lower_bound (Fallback)" begin
        @test_throws ArgumentError delete_lower_bound(dvref)
    end
    # test JuMP.delete_lower_bound (GeneralVariableRef)
    @testset "JuMP.delete_lower_bound (GeneralVariableRef)" begin
        @test_throws ArgumentError delete_lower_bound(gvref)
    end
end

# test Upper Bound Methods
@testset "Upper Bound Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test JuMP.has_upper_bound (Fallback)
    @testset "JuMP.has_upper_bound (Fallback)" begin
        @test_throws ArgumentError has_upper_bound(dvref)
    end
    # test JuMP.has_upper_bound (GeneralVariableRef)
    @testset "JuMP.has_upper_bound (GeneralVariableRef)" begin
        @test_throws ArgumentError has_upper_bound(gvref)
    end
    # test JuMP.upper_bound (Fallback)
    @testset "JuMP.upper_bound (Fallback)" begin
        @test_throws ArgumentError upper_bound(dvref)
    end
    # test JuMP.upper_bound (GeneralVariableRef)
    @testset "JuMP.upper_bound (GeneralVariableRef)" begin
        @test_throws ArgumentError upper_bound(gvref)
    end
    # test JuMP.set_upper_bound (Fallback)
    @testset "JuMP.set_upper_bound (Fallback)" begin
        @test_throws ArgumentError set_upper_bound(dvref, 42)
    end
    # test JuMP.set_upper_bound (GeneralVariableRef)
    @testset "JuMP.set_upper_bound (GeneralVariableRef)" begin
        @test_throws ArgumentError set_upper_bound(gvref, 42)
    end
    # test JuMP.UpperBoundRef (Fallback)
    @testset "JuMP.UpperBoundRef (Fallback)" begin
        @test_throws ArgumentError UpperBoundRef(dvref)
    end
    # test JuMP.UpperBoundRef (GeneralVariableRef)
    @testset "JuMP.UpperBoundRef (GeneralVariableRef)" begin
        @test_throws ArgumentError UpperBoundRef(gvref)
    end
    # test JuMP.delete_upper_bound (Fallback)
    @testset "JuMP.delete_upper_bound (Fallback)" begin
        @test_throws ArgumentError delete_upper_bound(dvref)
    end
    # test JuMP.delete_upper_bound (GeneralVariableRef)
    @testset "JuMP.delete_upper_bound (GeneralVariableRef)" begin
        @test_throws ArgumentError delete_upper_bound(gvref)
    end
end

# test Fixing Methods
@testset "Fixing Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test JuMP.is_fixed (Fallback)
    @testset "JuMP.is_fixed (Fallback)" begin
        @test_throws ArgumentError is_fixed(dvref)
    end
    # test JuMP.is_fixed (GeneralVariableRef)
    @testset "JuMP.is_fixed (GeneralVariableRef)" begin
        @test_throws ArgumentError is_fixed(gvref)
    end
    # test JuMP.fix_value (Fallback)
    @testset "JuMP.fix_value (Fallback)" begin
        @test_throws ArgumentError fix_value(dvref)
    end
    # test JuMP.fix_value (GeneralVariableRef)
    @testset "JuMP.fix_value (GeneralVariableRef)" begin
        @test_throws ArgumentError fix_value(gvref)
    end
    # test JuMP.fix (Fallback)
    @testset "JuMP.fix (Fallback)" begin
        @test_throws ArgumentError fix(dvref, 42, force = true)
    end
    # test JuMP.fix (GeneralVariableRef)
    @testset "JuMP.fix (GeneralVariableRef)" begin
        @test_throws ArgumentError fix(gvref, 42, force = true)
    end
    # test JuMP.FixRef (Fallback)
    @testset "JuMP.FixRef (Fallback)" begin
        @test_throws ArgumentError FixRef(dvref)
    end
    # test JuMP.FixRef (GeneralVariableRef)
    @testset "JuMP.FixRef (GeneralVariableRef)" begin
        @test_throws ArgumentError FixRef(gvref)
    end
    # test JuMP.unfix (Fallback)
    @testset "JuMP.unfix (Fallback)" begin
        @test_throws ArgumentError unfix(dvref)
    end
    # test JuMP.unfix (GeneralVariableRef)
    @testset "JuMP.unfix (GeneralVariableRef)" begin
        @test_throws ArgumentError unfix(gvref)
    end
end

# test Start Value Methods
@testset "Start Value Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test JuMP.start_value (Fallback)
    @testset "JuMP.start_value (Fallback)" begin
        @test_throws ArgumentError start_value(dvref)
    end
    # test JuMP.start_value (GeneralVariableRef)
    @testset "JuMP.start_value (GeneralVariableRef)" begin
        @test_throws ArgumentError start_value(gvref)
    end
    # test JuMP.set_start_value (Fallback)
    @testset "JuMP.set_start_value(Fallback)" begin
        @test_throws ArgumentError set_start_value(dvref, 42)
    end
    # test JuMP.set_start_value (GeneralVariableRef)
    @testset "JuMP.set_start_value (GeneralVariableRef)" begin
        @test_throws ArgumentError set_start_value(gvref, 42)
    end
end

# test Binary Methods
@testset "Binary Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test JuMP.is_binary (Fallback)
    @testset "JuMP.is_binary (Fallback)" begin
        @test_throws ArgumentError is_binary(dvref)
    end
    # test JuMP.is_binary (GeneralVariableRef)
    @testset "JuMP.is_binary (GeneralVariableRef)" begin
        @test_throws ArgumentError is_binary(gvref)
    end
    # test JuMP.set_binary (Fallback)
    @testset "JuMP.set_binary (Fallback)" begin
        @test_throws ArgumentError set_binary(dvref)
    end
    # test JuMP.set_binary (GeneralVariableRef)
    @testset "JuMP.set_binary (GeneralVariableRef)" begin
        @test_throws ArgumentError set_binary(gvref)
    end
    # test JuMP.BinaryRef (Fallback)
    @testset "JuMP.BinaryRef (Fallback)" begin
        @test_throws ArgumentError BinaryRef(dvref)
    end
    # test JuMP.BinaryRef (GeneralVariableRef)
    @testset "JuMP.BinaryRef (GeneralVariableRef)" begin
        @test_throws ArgumentError BinaryRef(gvref)
    end
    # test JuMP.unset_binary (Fallback)
    @testset "JuMP.unset_binary (Fallback)" begin
        @test_throws ArgumentError unset_binary(dvref)
    end
    # test JuMP.unset_binary (GeneralVariableRef)
    @testset "JuMP.unset_binary (GeneralVariableRef)" begin
        @test_throws ArgumentError unset_binary(gvref)
    end
end

# test Integer Methods
@testset "Integer Methods" begin
    # Setup data
    m = InfiniteModel();
    idx = TestIndex(1)
    dvref = TestVariableRef(m, idx)
    gvref = GeneralVariableRef(m, 1, TestIndex)
    # test JuMP.is_integer (Fallback)
    @testset "JuMP.is_integer (Fallback)" begin
        @test_throws ArgumentError is_integer(dvref)
    end
    # test JuMP.is_integer (GeneralVariableRef)
    @testset "JuMP.is_integer (GeneralVariableRef)" begin
        @test_throws ArgumentError is_integer(gvref)
    end
    # test JuMP.set_integer (Fallback)
    @testset "JuMP.set_integer (Fallback)" begin
        @test_throws ArgumentError set_integer(dvref)
    end
    # test JuMP.set_integer (GeneralVariableRef)
    @testset "JuMP.set_integer (GeneralVariableRef)" begin
        @test_throws ArgumentError set_integer(gvref)
    end
    # test JuMP.IntegerRef (Fallback)
    @testset "JuMP.IntegerRef (Fallback)" begin
        @test_throws ArgumentError IntegerRef(dvref)
    end
    # test JuMP.IntegerRef (GeneralVariableRef)
    @testset "JuMP.IntegerRef (GeneralVariableRef)" begin
        @test_throws ArgumentError IntegerRef(gvref)
    end
    # test JuMP.unset_integer (Fallback)
    @testset "JuMP.unset_integer (Fallback)" begin
        @test_throws ArgumentError unset_integer(dvref)
    end
    # test JuMP.unset_integer (GeneralVariableRef)
    @testset "JuMP.unset_integer (GeneralVariableRef)" begin
        @test_throws ArgumentError unset_integer(gvref)
    end
end
