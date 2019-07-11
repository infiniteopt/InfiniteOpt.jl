# Test extensions to basic Base methods
@testset "Base Extensions" begin
    m = InfiniteModel()
    m2 = InfiniteModel()
    ivref = InfiniteVariableRef(m, 1)
    pvref = PointVariableRef(m, 2)
    gvref = GlobalVariableRef(m, 3)
    pref = ParameterRef(m, 1)
    # variable compare
    @testset "(==)" begin
        @test ivref == ivref
        @test pvref == pvref
        @test gvref == gvref
        @test ivref == InfiniteVariableRef(m, 1)
        @test pvref == PointVariableRef(m, 2)
        @test gvref == GlobalVariableRef(m, 3)
        @test !(ivref == InfiniteVariableRef(m, 2))
        @test !(ivref == InfiniteVariableRef(m2, 1))
        @test !(ivref != InfiniteVariableRef(m, 1))
        @test !(pref == ivref)
    end
    # copy(v)
    @testset "copy(v)" begin
        @test copy(ivref) == ivref
        @test copy(pvref) == pvref
        @test copy(gvref) == gvref
    end
    # copy(v, m)
    @testset "copy(v, m)" begin
        @test copy(ivref, m2) == InfiniteVariableRef(m2, 1)
        @test copy(pvref, m2) == PointVariableRef(m2, 2)
        @test copy(gvref, m2) == GlobalVariableRef(m2, 3)
    end
    # broadcastable
    @testset "broadcastable" begin
        @test isa(Base.broadcastable(ivref), Base.RefValue{InfiniteVariableRef})
        @test isa(Base.broadcastable(pvref), Base.RefValue{PointVariableRef})
        @test isa(Base.broadcastable(gvref), Base.RefValue{GlobalVariableRef})
    end
end

# Test core JuMP methods
@testset "Core JuMP Extensions" begin
    m = InfiniteModel()
    m2 = InfiniteModel()
    ivref = InfiniteVariableRef(m, 1)
    pvref = PointVariableRef(m, 2)
    gvref = GlobalVariableRef(m, 3)
    pref = ParameterRef(m, 1)
    # isequal_canonical
    @testset "JuMP.isequal_canonical" begin
        @test isequal_canonical(ivref, ivref)
        @test isequal_canonical(pvref, pvref)
        @test isequal_canonical(gvref, gvref)
        @test !isequal_canonical(ivref, InfiniteVariableRef(m2, 1))
        @test !isequal_canonical(ivref, InfiniteVariableRef(m, 2))
    end
    # variable_type(m)
    @testset "JuMP.variable_type(m)" begin
        @test variable_type(m) == GeneralVariableRef
    end
    # variable_type(m, t)
    @testset "JuMP.variable_type(m, t)" begin
        @test variable_type(m, Infinite) == InfiniteVariableRef
        @test variable_type(m, Point) == PointVariableRef
        @test variable_type(m, Global) == GlobalVariableRef
        @test variable_type(m, Parameter) == ParameterRef
        @test_throws ErrorException variable_type(m, :bad)
    end
end

# Test precursor functions needed for add_parameter
@testset "Basic Reference Queries" begin
    m = InfiniteModel()
    ivref = InfiniteVariableRef(m, 1)
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    m.vars[1] = InfiniteVariable(info, (pref, ))
    # JuMP.index
    @testset "JuMP.index" begin
        @test index(ivref) == 1
    end
    # JuMP.owner_model
    @testset "JuMP.owner_model" begin
        @test owner_model(ivref) == m
    end
    # JuMP.is_valid
    @testset "JuMP.is_valid" begin
        @test is_valid(m, ivref)
        @test !is_valid(InfiniteModel(), ivref)
        @test !is_valid(m, InfiniteVariableRef(m, 5))
    end
end

# Test name methods
@testset "Infinite Variable Name" begin
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test1")
    pref2 = add_parameter(m, param, "test2")
    var = InfiniteVariable(info, (pref, pref2))
    m.vars[1] = var
    m.var_to_name[1] = "var"
    vref = InfiniteVariableRef(m, 1)
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "var"
    end
    # parameter_refs
    @testset "parameter_refs" begin
        @test parameter_refs(vref) == (pref, pref2)
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        vref2 = InfiniteVariableRef(m, 2)
        m.vars[2] = var
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new(test1, test2)"
        @test isa(set_name(vref2, ""), Nothing)
        @test name(vref2) == "noname(test1, test2)"
    end
    # parameter_by_name
    @testset "JuMP.variable_by_name" begin
        @test variable_by_name(m, "new(test1, test2)") == vref
        @test isa(variable_by_name(m, "test(test1, test2)"), Nothing)
        m.vars[2] = var
        m.var_to_name[2] = "new(test1, test2)"
        m.name_to_var = nothing
        @test_throws ErrorException variable_by_name(m, "new(test1, test2)")
    end
    # _root_name
    @testset "_root_name" begin
        @test InfiniteOpt._root_name(vref) == "new"
    end
end

# TODO test variable info methods here

# Test variable definition methods
@testset "Infinite Variable Definition" begin
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    pref2 = add_parameter(m, param, "θ")
    prefs = @infinite_parameter(m, x[1:2], set = IntervalSet(0, 1))
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    # _check_parameter_tuple
    @testset "_check_parameter_tuple" begin
        @test isa(InfiniteOpt._check_parameter_tuple(error, (pref, prefs)),
                  Nothing)
        @test_throws ErrorException InfiniteOpt._check_parameter_tuple(error,
                                                               (pref, prefs, 2))
    end
    # _make_formatted_tuple
    @testset "_make_formatted_tuple" begin
        @test isa(InfiniteOpt._make_formatted_tuple((pref, prefs)), Tuple)
        @test isa(InfiniteOpt._make_formatted_tuple((pref, prefs))[2],
                  JuMP.Containers.SparseAxisArray)
        @test isa(InfiniteOpt._make_formatted_tuple((pref, prefs))[1],
                  ParameterRef)
    end
    # _check_tuple_groups
    @testset "_check_tuple_groups" begin
        tuple = InfiniteOpt._make_formatted_tuple((pref, prefs))
        @test isa(InfiniteOpt._check_tuple_groups(error, tuple), Nothing)
        tuple = InfiniteOpt._make_formatted_tuple(([pref; pref2], prefs))
        @test_throws ErrorException InfiniteOpt._check_tuple_groups(error,
                                                                    tuple)
        @test_throws ErrorException InfiniteOpt._check_tuple_groups(error,
                                                                   (pref, pref))
    end
    # build_variable
    @testset "JuMP.build_variable" begin
        @test_throws ErrorException build_variable(error, info, Infinite,
                                                   bob = 42)
        @test_throws ErrorException build_variable(error, info, :bad)
        @test_throws ErrorException build_variable(error, info, Point,
                                                   parameter_refs = pref)
        @test_throws ErrorException build_variable(error, info, Infinite)
        @test_throws ErrorException build_variable(error, info, Infinite,
                                                   parameter_refs = (pref, 2))
        @test_throws ErrorException build_variable(error, info, Infinite,
                                                  parameter_refs = (pref, pref),
                                                  error = error)
        expected = InfiniteVariable(info, (pref,))
        @test build_variable(error, info, Infinite,
                             parameter_refs = pref).info == expected.info
        @test build_variable(error, info, Infinite,
                parameter_refs = pref).parameter_refs == expected.parameter_refs
        @test build_variable(error, info, Infinite,
                 parameter_refs = (pref, pref2)).parameter_refs == (pref, pref2)
        tuple = InfiniteOpt._make_formatted_tuple((pref, prefs))
        @test build_variable(error, info, Infinite,
                         parameter_refs = (pref, prefs)).parameter_refs == tuple
        tuple = InfiniteOpt._make_formatted_tuple((prefs,))
        @test build_variable(error, info, Infinite,
                             parameter_refs = prefs).parameter_refs == tuple
    end
    # _update_param_var_mapping
    @testset "_update_param_var_mapping" begin
        m2 = InfiniteModel()
        param = InfOptParameter(IntervalSet(0, 1), Number[], false)
        pref3 = add_parameter(m2, param, "test")
        prefs2 = @infinite_parameter(m2, x[1:2], set = IntervalSet(0, 1))
        ivref = InfiniteVariableRef(m2, 1)
        ivref2 = InfiniteVariableRef(m2, 2)
        tuple = (pref3, prefs2)
        tuple = InfiniteOpt._make_formatted_tuple(tuple)
        @test isa(InfiniteOpt._update_param_var_mapping(ivref, tuple), Nothing)
        @test m2.param_to_vars[1] == [1]
        @test m2.param_to_vars[2] == [1]
        @test m2.param_to_vars[3] == [1]
        @test isa(InfiniteOpt._update_param_var_mapping(ivref2, tuple), Nothing)
        @test m2.param_to_vars[1] == [1, 2]
        @test m2.param_to_vars[2] == [1, 2]
        @test m2.param_to_vars[3] == [1, 2]
    end
    # _check_parameters_valid
    @testset "_check_parameters_valid" begin
        tuple = (pref, prefs, copy(pref2, InfiniteModel()))
        tuple = InfiniteOpt._make_formatted_tuple(tuple)
        @test_throws ErrorException InfiniteOpt._check_parameters_valid(m, tuple)
        @test isa(InfiniteOpt._check_parameters_valid(m, (pref, pref2)), Nothing)
    end
    # add_variable
    # TODO test info updates
    @testset "JuMP.add_variable" begin
        m2 = InfiniteModel()
        param = InfOptParameter(IntervalSet(0, 1), Number[], false)
        pref3 = add_parameter(m2, param, "test")
        v = build_variable(error, info, Infinite,
                           parameter_refs = pref3)
        @test_throws ErrorException add_variable(m, v)
        v = build_variable(error, info, Infinite, parameter_refs = pref)
        @test add_variable(m, v, "name") == InfiniteVariableRef(m, 2)
        @test haskey(m.vars, 2)
        @test m.param_to_vars[1] == [2]
        @test m.var_to_name[2] == "name(test)"
    end
end

# Test name methods
@testset "Point Variable Name" begin
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test1")
    pref2 = add_parameter(m, param, "test2")
    ivar = InfiniteVariable(info, (pref, pref2))
    ivref = add_variable(m, ivar, "ivar")
    var = PointVariable(info, ivref, (0.5, 0.5))
    m.vars[2] = var
    m.var_to_name[2] = "var"
    vref = PointVariableRef(m, 2)
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "var"
    end
    # infinite_variable_ref
    @testset "infinite_variable_ref" begin
        @test infinite_variable_ref(vref) == ivref
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        vref2 = PointVariableRef(m, 3)
        m.vars[3] = var
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
        @test isa(set_name(vref2, ""), Nothing)
        @test name(vref2) == "ivar(0.5, 0.5)"
    end
    # parameter_by_name
    @testset "JuMP.variable_by_name" begin
        @test variable_by_name(m, "new") == vref
        @test isa(variable_by_name(m, "test"), Nothing)
        m.vars[3] = var
        m.var_to_name[3] = "new"
        m.name_to_var = nothing
        @test_throws ErrorException variable_by_name(m, "new")
    end
end

# Test variable definition methods
@testset "Point Variable Definition" begin
    m = InfiniteModel()
    param = InfOptParameter(IntervalSet(0, 1), Number[], false)
    pref = add_parameter(m, param, "test")
    pref2 = add_parameter(m, param, "θ")
    prefs = @infinite_parameter(m, x[1:2], set = IntervalSet(0, 1))
    info = VariableInfo(false, 0., false, 0., false, 0., false, 0., false, false)
    ivar = InfiniteVariable(info, (pref, pref2))
    ivref = add_variable(m, ivar, "ivar")
    ivar2 = build_variable(error, info, Infinite, parameter_refs = (pref, prefs))
    ivref2 = add_variable(m, ivar2, "ivar2")
    # _check_tuple_shape
    @testset "_check_tuple_shape" begin
        @test isa(InfiniteOpt._check_tuple_shape(error, ivref, (0.5, 0.5)),
                  Nothing)
        tuple = InfiniteOpt._make_formatted_tuple((0.5, [0.5, 0.5]))
        @test isa(InfiniteOpt._check_tuple_shape(error, ivref2, tuple), Nothing)
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, ivref,
                                                                   (0.5,))
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, ivref,
                                                                   (0.5, [0.5]))
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, ivref2,
                                                                   (0.5, 0.5))
        tuple = InfiniteOpt._make_formatted_tuple((0.5, [0.5, 0.5, 0.5]))
        @test_throws ErrorException InfiniteOpt._check_tuple_shape(error, ivref2,
                                                                   tuple)
    end
    # _check_tuple_values
    @testset "_check_tuple_values" begin
        @test isa(InfiniteOpt._check_tuple_values(error, ivref, (0.5, 0.5)),
                  Nothing)
        tuple = InfiniteOpt._make_formatted_tuple((0, [0.5, 1]))
        @test isa(InfiniteOpt._check_tuple_values(error, ivref2, tuple), Nothing)
        @test_throws ErrorException InfiniteOpt._check_tuple_values(error, ivref,
                                                                    (0, 2))
        tuple = InfiniteOpt._make_formatted_tuple((0, [2, 1]))
        @test_throws ErrorException InfiniteOpt._check_tuple_values(error, ivref2,
                                                                    tuple)
    end
    # _update_point_info
    @testset "_update_point_info" begin
        new_info = VariableInfo(true, 0., true, 0., false, 0., true, 0., true,
                                false)
        InfiniteOpt._update_variable_info(ivref, new_info)
        @test InfiniteOpt._update_point_info(info, ivref) == new_info
        new_info = VariableInfo(false, 0., false, 0., true, 0., true, 0., false,
                                true)
        InfiniteOpt._update_variable_info(ivref, new_info)
        @test InfiniteOpt._update_point_info(info, ivref) == new_info
        curr_info = VariableInfo(true, 0., true, 0., false, 0., true, 0., true,
                                 false)
        @test InfiniteOpt._update_point_info(curr_info, ivref) == curr_info
    end
    # build_variable
    @testset "JuMP.build_variable" begin
        @test_throws ErrorException build_variable(error, info, Infinite,
                                                  infinite_variable_ref = ivref)
        @test_throws ErrorException build_variable(error, info, Infinite,
                                                   parameter_values = 3)
        @test_throws ErrorException build_variable(error, info, Point)
        @test_throws ErrorException build_variable(error, info, Point,
                                                  infinite_variable_ref = ivref)
        @test_throws ErrorException build_variable(error, info, Point,
                                                   parameter_values = 3)
        @test build_variable(error, info, Point, infinite_variable_ref = ivref,
                   parameter_values = (0.5, 0.5)).infinite_variable_ref == ivref
        @test build_variable(error, info, Point, infinite_variable_ref = ivref,
                   parameter_values = (0.5, 0.5)).parameter_values == (0.5, 0.5)
        @test_throws ErrorException build_variable(error, info, Point,
                                                  infinite_variable_ref = ivref,
                                                  parameter_values = (0.5, 2))
        @test build_variable(error, info, Point, infinite_variable_ref = ivref2,
               parameter_values = (0.5, [0, 0])).infinite_variable_ref == ivref2
        tuple = tuple = InfiniteOpt._make_formatted_tuple((0.5, [0, 0]))
        @test build_variable(error, info, Point, infinite_variable_ref = ivref2,
                     parameter_values = (0.5, [0, 0])).parameter_values == tuple
        @test_throws ErrorException build_variable(error, info, Point,
                                            infinite_variable_ref = ivref2,
                                            parameter_values = (0.5, [0, 0, 0]))
    end
    # _update_param_supports
    @testset "_update_param_supports" begin
        @test isa(InfiniteOpt._update_param_supports(ivref, (0.5, 1)), Nothing)
        @test supports(pref) == [0.5]
        @test supports(pref2) == [1]
        tuple = InfiniteOpt._make_formatted_tuple((0.5, [0, 1]))
        @test isa(InfiniteOpt._update_param_supports(ivref2, tuple), Nothing)
        @test supports(pref) == [0.5]
        @test supports(prefs[1]) == [0]
        @test supports(prefs[2]) == [1]
    end
    # add_variable
    # TODO test info updates
    @testset "JuMP.add_variable" begin
        m2 = InfiniteModel()
        pref3 = add_parameter(m2, param, "test")
        ivar3 = InfiniteVariable(info, (pref3,))
        ivref3 = add_variable(m2, ivar3, "ivar")
        v = build_variable(error, info, Point, infinite_variable_ref = ivref3,
                           parameter_values = 0.5)
        @test_throws ErrorException add_variable(m, v)
        v = build_variable(error, info, Point, infinite_variable_ref = ivref,
                           parameter_values = (0, 1))
        @test add_variable(m, v, "name") == PointVariableRef(m, 4)
        @test haskey(m.vars, 4)
        @test supports(pref) == [0, 0.5]
        @test supports(pref2) == [1]
        @test m.var_to_name[4] == "name"
    end
end

# Test name methods
@testset "Global Variable Name" begin
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    var = GlobalVariable(info)
    m.vars[1] = var
    m.var_to_name[1] = "test"
    vref = GlobalVariableRef(m, 1)
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(vref) == "test"
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(vref, "new"), Nothing)
        @test name(vref) == "new"
    end
    # parameter_by_name
    @testset "JuMP.variable_by_name" begin
        @test variable_by_name(m, "new") == vref
        @test isa(variable_by_name(m, "test2"), Nothing)
        m.vars[2] = var
        m.var_to_name[2] = "new"
        m.name_to_var = nothing
        @test_throws ErrorException variable_by_name(m, "new")
    end
end

# Test variable definition methods
@testset "Global Variable Definition" begin
    m = InfiniteModel()
    info = VariableInfo(false, 0, false, 0, false, 0, false, 0, false, false)
    # build_variable
    @testset "JuMP.build_variable" begin
        expected = GlobalVariable(info)
        @test build_variable(error, info, Global) == expected
    end
    # add_variable
    # TODO test info updates
    @testset "JuMP.add_variable" begin
        v = build_variable(error, info, Global)
        @test add_variable(m, v, "name") == GlobalVariableRef(m, 1)
        @test haskey(m.vars, 1)
        @test m.var_to_name[1] == "name"
    end
end
