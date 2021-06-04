# Test basic extensions
@testset "Basics" begin
    # initialize the model, references, and other info
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 10])
    con = ScalarConstraint(zero(AffExpr), MOI.EqualTo(0.0))
    object = ConstraintData(con, Int[], "c1", MeasureIndex[], false)
    idx = InfOptConstraintIndex(1)
    cref = InfOptConstraintRef(m, idx)
    bad_cref = InfOptConstraintRef(m, InfOptConstraintIndex(-1))
    # owner_model
    @testset "JuMP.owner_model" begin
      @test owner_model(cref) === m
    end
    # index
    @testset "JuMP.index" begin
      @test JuMP.index(cref) == InfOptConstraintIndex(1)
    end
    # test Base.:(==) of constraint references
    @testset "Base.:(==) References" begin
        @test !(cref == bad_cref)
        @test cref == InfOptConstraintRef(m, idx)
        @test cref != InfOptConstraintRef(InfiniteModel(), idx)
    end
    # test broadcastable
    @testset "Base.broadcastable Reference" begin
        @test isa(Base.broadcastable(cref), Base.RefValue)
    end
    # _add_data_object
    @testset "_add_data_object" begin
        @test InfiniteOpt._add_data_object(m, object) == idx
    end
    # _data_dictionary
    @testset "_data_dictionary" begin
        @test InfiniteOpt._data_dictionary(cref) === m.constraints
    end
    # is_valid
    @testset "JuMP.is_valid" begin
      @test is_valid(m, cref)
    end
    # _data_object
    @testset "_data_object" begin
        @test InfiniteOpt._data_object(cref) === object
        @test_throws ErrorException InfiniteOpt._data_object(bad_cref)
    end
    # _core_constraint_object
    @testset "_core_constraint_object" begin
        @test InfiniteOpt._core_constraint_object(cref) === con
    end
    # _set_core_constraint_object
    @testset "_set_core_constraint_object" begin
        @test InfiniteOpt._set_core_constraint_object(cref, con) isa Nothing
        con2 = ScalarConstraint(zero(AffExpr), MOI.LessThan(0.0))
        @test InfiniteOpt._set_core_constraint_object(cref, con2) isa Nothing
        @test InfiniteOpt._core_constraint_object(cref) === con2
        @test InfiniteOpt._set_core_constraint_object(cref, con) isa Nothing
        @test InfiniteOpt._core_constraint_object(cref) === con
    end
    # _object_numbers
    @testset "_object_numbers" begin
        @test InfiniteOpt._object_numbers(cref) == []
    end
    # _measure_dependencies
    @testset "_measure_dependencies" begin
        @test InfiniteOpt._measure_dependencies(cref) == MeasureIndex[]
    end
    @testset "_is_info_constraint" begin
        @test !InfiniteOpt._is_info_constraint(cref)
    end
    # constraint_object
    @testset "JuMP.constraint_object" begin
      @test constraint_object(cref) == con
    end
    # JuMP.name
    @testset "JuMP.name" begin
        @test name(cref) == "c1"
        @test name(bad_cref) == ""
    end
    # JuMP.set_name
    @testset "JuMP.set_name" begin
        @test isa(set_name(cref, "new"), Nothing)
        @test_throws ErrorException set_name(bad_cref, "")
        @test name(cref) == "new"
    end
    # test has_domain_restrictions
    @testset "has_domain_restrictions" begin
        @test !has_domain_restrictions(cref)
    end
    # test domain_restrictions
    @testset "domain_restrictions" begin
        @test domain_restrictions(cref) == DomainRestrictions()
    end
    # _make_constraint_ref
    @testset "_make_constraint_ref" begin
        @test InfiniteOpt._make_constraint_ref(m, idx) == cref
    end
    # constraint_by_name
    @testset "JuMP.constraint_by_name" begin
        # test normal
        @test constraint_by_name(m, "new") == cref
        @test isa(constraint_by_name(m, "test2"), Nothing)
        # prepare constraint with duplicate name
        idx2 = InfiniteOpt._add_data_object(m, object)
        cref2 = InfOptConstraintRef(m, idx2)
        @test isa(set_name(cref2, "new"), Nothing)
        # test for duplciate name error
        @test_throws ErrorException constraint_by_name(m, "new")
    end
    # _delete_data_object
    @testset "_delete_data_object" begin
        @test InfiniteOpt._delete_data_object(cref) isa Nothing
        @test length(InfiniteOpt._data_dictionary(cref)) == 1
        @test !is_valid(m, cref)
    end
end

# Test definition methods
@testset "Definition" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @infinite_parameter(m, pars[1:2] in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, inf2, Infinite(par, pars))
    @variable(m, pt, Point(inf, 0.5))
    @variable(m, x)
    @variable(m, rv, SemiInfinite(inf2, 0.5, pars))
    data = TestData(par, 1, 1)
    meas = @measure(inf + par - x, data, name = "test")
    mindex = MeasureIndex(1)
    dinf = @deriv(inf, par)
    @finite_parameter(m, fin == 42)
    # test _check_restrictions
    @testset "_check_restrictions" begin
        # test normal
        rs = DomainRestrictions(par => [0, 1])
        @test InfiniteOpt._check_restrictions(rs) isa Nothing
        rs = DomainRestrictions(pars => [0, 1])
        @test InfiniteOpt._check_restrictions(rs) isa Nothing
        rs = DomainRestrictions(pars => 0)
        @test InfiniteOpt._check_restrictions(rs) isa Nothing
        rs = DomainRestrictions(pars[1] => [0, 1])
        @test InfiniteOpt._check_restrictions(rs) isa Nothing
        # test errors
        rs = DomainRestrictions(fin => [0, 1])
        @test_throws ErrorException InfiniteOpt._check_restrictions(rs)
        rs = DomainRestrictions(par => [-1, 1])
        @test_throws ErrorException InfiniteOpt._check_restrictions(rs)
        rs = DomainRestrictions(par => [0, 2])
        @test_throws ErrorException InfiniteOpt._check_restrictions(rs)
        rs = DomainRestrictions(pars[1] => 0)
        @test_throws ErrorException InfiniteOpt._check_restrictions(rs)
        rs = DomainRestrictions(pars[1] => 0, pars[2] => [0, 3])
        @test_throws ErrorException InfiniteOpt._check_restrictions(rs)
    end
    # test build_constraint for DomainRestrictedConstraint
    @testset "build_constraint (DomainRestrictedConstraint)" begin 
        c = JuMP.ScalarConstraint(inf + x, MOI.EqualTo(0.0))
        rs = DomainRestrictions(par => [0, 1])
        @test build_constraint(error, inf + x, MOI.EqualTo(0.0), rs).constraint isa ScalarConstraint
        @test build_constraint(error, inf + x, MOI.EqualTo(0.0), rs).restrictions == rs
    end
    # test _validate_restrictions 
    @testset "_validate_restrictions" begin 
        # test normal
        rs = DomainRestrictions(par => [0, 1])
        @test InfiniteOpt._validate_restrictions(m, rs) isa Nothing
        @test num_supports(par) == 1
        # test error
        @infinite_parameter(InfiniteModel(), par2 ~ Normal())
        rs = DomainRestrictions(par2 => [0, 1])
        err = VariableNotOwned{GeneralVariableRef} 
        @test_throws err InfiniteOpt._validate_restrictions(m, rs)
        # test support addition
        rs = DomainRestrictions(par => 0)
        @test InfiniteOpt._validate_restrictions(m, rs) isa Nothing
        @test supports(par, label = UserDefined) == [0, 0.5]
        rs = DomainRestrictions(pars => 0)
        @test InfiniteOpt._validate_restrictions(m, rs) isa Nothing
        @test supports(pars) == zeros(2, 1)
    end
    # test _update_var_constr_mapping
    @testset "_update_var_constr_mapping" begin
        # add dumby constraints
        con = ScalarConstraint(zero(AffExpr), MOI.EqualTo(0.0))
        object1 = ConstraintData(con, Int[], "c", MeasureIndex[], false)
        object2 = ConstraintData(con, Int[], "c", MeasureIndex[], false)
        idx1 = InfOptConstraintIndex(1)
        idx2 = InfOptConstraintIndex(2)
        cref1 = InfOptConstraintRef(m, idx1)
        cref2 = InfOptConstraintRef(m, idx2)
        @test InfiniteOpt._add_data_object(m, object1) == idx1
        @test InfiniteOpt._add_data_object(m, object2) == idx2
        # test initial use of variable
        @test isa(InfiniteOpt._update_var_constr_mapping([inf, par, meas, rv], cref1),
                  Nothing)
        @test InfiniteOpt._constraint_dependencies(inf) == [idx1]
        @test InfiniteOpt._constraint_dependencies(par) == [idx1]
        @test InfiniteOpt._constraint_dependencies(meas) == [idx1]
        @test InfiniteOpt._constraint_dependencies(rv) == [idx1]
        @test InfiniteOpt._measure_dependencies(cref1) == [mindex]
        # test secondary use of variable
        @test isa(InfiniteOpt._update_var_constr_mapping([inf, par, meas, rv], cref2),
                  Nothing)
        @test InfiniteOpt._constraint_dependencies(inf) == [idx1, idx2]
        @test InfiniteOpt._constraint_dependencies(par) == [idx1, idx2]
        @test InfiniteOpt._constraint_dependencies(meas) == [idx1, idx2]
        @test InfiniteOpt._constraint_dependencies(rv) == [idx1, idx2]
        @test InfiniteOpt._measure_dependencies(cref2) == [mindex]
        # Undo changes
        empty!(InfiniteOpt._constraint_dependencies(inf))
        empty!(InfiniteOpt._constraint_dependencies(par))
        empty!(InfiniteOpt._constraint_dependencies(meas))
        empty!(InfiniteOpt._constraint_dependencies(rv))
        empty!(InfiniteOpt._measure_dependencies(cref1))
        empty!(InfiniteOpt._measure_dependencies(cref2))
        empty!(m.constraints)
    end
    # test add_constraint
    @testset "JuMP.add_constraint" begin
        # test rejects bad variables
        @variable(InfiniteModel(), z)
        con = ScalarConstraint(z, MOI.EqualTo(42.0))
        err = VariableNotOwned{GeneralVariableRef}
        @test_throws err add_constraint(m, con, "a")
        # test restricted constraint
        rs = DomainRestrictions(par => [0, 1])
        con = ScalarConstraint(inf + pt, MOI.EqualTo(42.0))
        con = DomainRestrictedConstraint(con, rs)
        idx = InfOptConstraintIndex(1)
        cref = InfOptConstraintRef(m, idx)
        @test add_constraint(m, con, "a") == cref
        @test name(cref) == "a"
        @test domain_restrictions(cref) == rs
        @test constraint_object(cref).set == MOI.EqualTo(42.0)
        # test infinite constraint
        con = ScalarConstraint(inf + pt, MOI.EqualTo(42.0))
        idx = InfOptConstraintIndex(2)
        cref = InfOptConstraintRef(m, idx)
        @test add_constraint(m, con, "b") == cref
        @test name(cref) == "b"
        @test domain_restrictions(cref) == DomainRestrictions()
        @test constraint_object(cref).set == MOI.EqualTo(42.0)
        # test finite constraint
        con = ScalarConstraint(x + pt, MOI.EqualTo(42.0))
        idx = InfOptConstraintIndex(3)
        cref = InfOptConstraintRef(m, idx)
        @test add_constraint(m, con, "d") == cref
        @test name(cref) == "d"
        @test !InfiniteOpt._is_info_constraint(cref)
        @test !optimizer_model_ready(m)
        @test used_by_constraint(pt)
        # test vector constraint
        con = VectorConstraint([inf + pt, 2inf], MOI.Zeros(2))
        idx = InfOptConstraintIndex(4)
        cref = InfOptConstraintRef(m, idx)
        @test add_constraint(m, con, "e") == cref
        @test constraint_object(cref) isa VectorConstraint
        @test !has_domain_restrictions(cref)
    end
    # test macro
    @testset "JuMP.@constraint" begin
        # test scalar constraint
        idx = InfOptConstraintIndex(5)
        cref = InfOptConstraintRef(m, idx)
        @test @constraint(m, f, x + pt -2 <= 2) == cref
        @test InfiniteOpt._core_constraint_object(cref) isa ScalarConstraint
        # test restricted scalar constraint
        rs = DomainRestrictions(par => [0, 1])
        idx = InfOptConstraintIndex(6)
        cref = InfOptConstraintRef(m, idx)
        @test @constraint(m, g, inf + meas - dinf <= 2, rs) == cref
        @test InfiniteOpt._object_numbers(cref) == [1]
        @test InfiniteOpt._core_constraint_object(cref) isa ScalarConstraint
        @test used_by_constraint(dinf)
    end
    # test deprecation
    @testset "@BDconstraint" begin 
        @test_macro_throws ErrorException @BDconstraint(m, bd(par == 0), x == 0)
    end
end

# Test parameter reference methods
@testset "Parameter References" begin
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @infinite_parameter(m, y in [0, 1])
    @infinite_parameter(m, x[1:3] in [0, 1])
    @variable(m, z)
    @constraint(m, c1, z >= 0)
    @constraint(m, c2, z + t + x[1] >= 0)
    # test parameter_refs
    @testset "parameter_refs" begin
        @test parameter_refs(c1) == ()
        @test parameter_refs(c2) == (t, x)
    end
end

# Test domain restrictions methods
@testset "Domain Restrictions" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 10])
    @infinite_parameter(m, pars[1:2] in [0, 10])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0.5))
    @variable(m, x >= 0, Int)
    @constraint(m, c1, inf + x == 0, DomainRestrictions(par => [0, 1]))
    @constraint(m, c2, x * pt + x == 2)
    # test has_domain_restrictions
    @testset "has_domain_restrictions" begin
        @test has_domain_restrictions(c1)
        @test !has_domain_restrictions(c2)
        @test_deprecated has_parameter_bounds(c1)
    end
    # test domain_restrictions
    @testset "domain_restrictions" begin
        @test domain_restrictions(c1) == DomainRestrictions(par => [0, 1])
        @test domain_restrictions(c2) == DomainRestrictions()
        @test_deprecated parameter_bounds(c1)
    end
    # test set_domain_restrictions
    @testset "set_domain_restrictions" begin
        # test force error
        rs = DomainRestrictions(par => [0, 1])
        @test_throws ErrorException set_domain_restrictions(c1, rs)
        # test normal
        @test set_domain_restrictions(c2, rs) isa Nothing
        @test domain_restrictions(c2) == rs
        @test !optimizer_model_ready(m)
        # test test error with restrictions
        rs2 = DomainRestrictions(par => [-1, 1])
        @test_throws ErrorException set_domain_restrictions(c1, rs2, force = true)
        # test deprecation
        empty!(domain_restrictions(c1).intervals)
        @test_deprecated set_parameter_bounds(c1, rs)
        @test_macro_throws ErrorException @set_parameter_bounds(c1, rs)
    end
    # test _update_restrictions
    @testset "_update_restrictions" begin
        # test normal
        rs1 = DomainRestrictions(pars[1] => [0, 1], pars[2] => [-1, 2])
        rs2 = DomainRestrictions(par => 0, pars[1] => [0.5, 1.5],
                                 pars[2] => [0, 1])
        @test isa(InfiniteOpt._update_restrictions(rs1, rs2), Nothing)
        @test length(rs1) == 3
        @test rs1[par] == IntervalDomain(0, 0)
        @test rs1[pars[1]] == IntervalDomain(0.5, 1)
        @test rs1[pars[2]] == IntervalDomain(0, 1)
        # test error
        new_rs = DomainRestrictions(par => 1)
        @test_throws ErrorException InfiniteOpt._update_restrictions(rs1, new_rs)
        new_rs = DomainRestrictions(par => -1)
        @test_throws ErrorException InfiniteOpt._update_restrictions(rs1, new_rs)
    end
    # test add_parameter_restrictions
    @testset "add_domain_restrictions" begin
        # test already has restrictions
        rs = DomainRestrictions(par => [0, 1])
        @test add_domain_restrictions(c1, copy(rs)) isa Nothing
        @test domain_restrictions(c1) == rs
        # test doesn't have restrictions
        @constraint(m, c3, inf + pt == 0)
        @test add_domain_restrictions(c3, copy(rs)) isa Nothing
        @test domain_restrictions(c3) == rs
        # test deprecation 
        @test_deprecated add_parameter_bounds(c1, copy(rs))
        @test_macro_throws ErrorException @add_parameter_bounds(c1, rs)
    end
    # test delete_domain_restrictions
    @testset "delete_domain_restrictions" begin
        # test partial deletion
        rs = DomainRestrictions(pars[1] => [0, 1])
        @test add_domain_restrictions(c1, rs) isa Nothing
        @test delete_domain_restrictions(c1) isa Nothing
        @test !has_domain_restrictions(c1)
        # test already gone
        @test delete_domain_restrictions(c1) isa Nothing
        @test !has_domain_restrictions(c1)
        # test deprecation 
        @test_deprecated delete_parameter_bounds(c1)
    end
end

# Test coefficient queries and modifications
@testset "Coefficient Methods" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf, Infinite(par))
    @variable(m, pt, Point(inf, 0.5))
    @variable(m, x >= 0., Int)
    @constraint(m, c1, inf + x == 0., DomainRestrictions(par => [0, 1]))
    @constraint(m, c2, x * pt + x == 2.)
    @constraint(m, c3, [inf, x] in MOI.Zeros(2))
    # test set_normalized_coefficient
    @testset "JuMP.set_normalized_coefficient" begin
        @test isa(set_normalized_coefficient(c1, x, 2), Nothing)
        @test jump_function(constraint_object(c1)) == inf + 2x
        @test isa(set_normalized_coefficient(c2, inf, 2), Nothing)
        @test jump_function(constraint_object(c2)) == x * pt + x + 2inf
        @test_throws ErrorException set_normalized_coefficient(c3, x, 42)
    end
    # test normalized_coefficient
    @testset "JuMP.normalized_coefficient" begin
        @test normalized_coefficient(c1, x) == 2.
        @test normalized_coefficient(c1, pt) == 0.
        @test normalized_coefficient(c2, inf) == 2.
        @test normalized_coefficient(LowerBoundRef(x), x) == 1.
        @test normalized_coefficient(LowerBoundRef(x), pt) == 0.
        @test_throws ErrorException normalized_coefficient(c3, x)
    end
    # test _set_set_value
    @testset "_set_set_value" begin
        @test InfiniteOpt._set_set_value(MOI.EqualTo(1.), 2) == MOI.EqualTo(2.)
        @test InfiniteOpt._set_set_value(MOI.GreaterThan(1.), 2) == MOI.GreaterThan(2.)
        @test InfiniteOpt._set_set_value(MOI.LessThan(1.), 2) == MOI.LessThan(2.)
    end
    # test set_normalized_rhs
    @testset "JuMP.set_normalized_rhs" begin
        @test isa(set_normalized_rhs(c1, 42.), Nothing)
        @test constraint_object(c1).set == MOI.EqualTo(42.0)
        @test isa(set_normalized_rhs(c2, 42.), Nothing)
        @test constraint_object(c2).set == MOI.EqualTo(42.0)
        @test isa(set_normalized_rhs(LowerBoundRef(x), 42.0), Nothing)
        @test constraint_object(LowerBoundRef(x)).set == MOI.GreaterThan(42.0)
        @test_throws ErrorException set_normalized_rhs(c3, 42)
    end
    # test normalized_rhs
    @testset "JuMP.normalized_rhs" begin
        @test normalized_rhs(c1) == 42.
        @test normalized_rhs(c2) == 42.
        @test normalized_rhs(LowerBoundRef(x)) == 42.
        @test_throws ErrorException normalized_rhs(c3)
    end
    # test add_to_function_constant
    @testset "JuMP.add_to_function_constant" begin
        @test isa(add_to_function_constant(c1, 4.), Nothing)
        @test constraint_object(c1).set == MOI.EqualTo(38.0)
        @test isa(add_to_function_constant(c2, 4.), Nothing)
        @test constraint_object(c2).set == MOI.EqualTo(38.0)
        @test isa(add_to_function_constant(LowerBoundRef(x), 4.0), Nothing)
        @test constraint_object(LowerBoundRef(x)).set == MOI.GreaterThan(38.0)
        @test_throws ErrorException add_to_function_constant(c3, 42)
    end
end

# Test num_constraints extensions
@testset "JuMP.num_constraints Extensions" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf >= 0, Infinite(par))
    @variable(m, pt, Point(inf, 0.5))
    @variable(m, 0 <= x <= 1, Int)
    # test search function type and set type
    @testset "Function and Set" begin
        @test num_constraints(m, GeneralVariableRef, MOI.LessThan) == 1
        @test num_constraints(m, GeneralVariableRef, MOI.GreaterThan) == 3
    end
    # test search function type
    @testset "Function" begin
        @test num_constraints(m, GeneralVariableRef) == 5
    end
    # test search set type
    @testset "Set" begin
        @test num_constraints(m, MOI.LessThan) == 1
        @test num_constraints(m, MOI.GreaterThan) == 3
    end
    # test search total
    @testset "Total" begin
        @test num_constraints(m) == 5
    end
end

# Test all_constraints extensions
@testset "JuMP.all_constraints Extensions" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf >= 0, Infinite(par))
    @variable(m, pt, Point(inf, 0.5))
    @variable(m, 0 <= x <= 1, Int)
    # test search function type and set type
    @testset "Function and Set" begin
        list = [InfOptConstraintRef(m, InfOptConstraintIndex(4))]
        @test all_constraints(m, GeneralVariableRef, MOI.LessThan) == list
        list = [InfOptConstraintRef(m, InfOptConstraintIndex(1)),
                InfOptConstraintRef(m, InfOptConstraintIndex(2)),
                InfOptConstraintRef(m, InfOptConstraintIndex(3))]
        @test all_constraints(m, GeneralVariableRef, MOI.GreaterThan) == list
    end
    # test search function type
    @testset "Function" begin
        list = [InfOptConstraintRef(m, InfOptConstraintIndex(1)),
                InfOptConstraintRef(m, InfOptConstraintIndex(2)),
                InfOptConstraintRef(m, InfOptConstraintIndex(3)),
                InfOptConstraintRef(m, InfOptConstraintIndex(4)),
                InfOptConstraintRef(m, InfOptConstraintIndex(5))]
        @test all_constraints(m, GeneralVariableRef) == list
    end
    # test search set type
    @testset "Set" begin
        list = [InfOptConstraintRef(m, InfOptConstraintIndex(4))]
        @test all_constraints(m, MOI.LessThan) == list
        list = [InfOptConstraintRef(m, InfOptConstraintIndex(1)),
                InfOptConstraintRef(m, InfOptConstraintIndex(2)),
                InfOptConstraintRef(m, InfOptConstraintIndex(3))]
        @test all_constraints(m, MOI.GreaterThan) == list
    end
    # test search total
    @testset "Total" begin
        list = [InfOptConstraintRef(m, InfOptConstraintIndex(1)),
                InfOptConstraintRef(m, InfOptConstraintIndex(2)),
                InfOptConstraintRef(m, InfOptConstraintIndex(3)),
                InfOptConstraintRef(m, InfOptConstraintIndex(4)),
                InfOptConstraintRef(m, InfOptConstraintIndex(5))]
        @test all_constraints(m) == list
    end
end

# Test list_of_constraint_types
@testset "JuMP.list_of_constraint_types" begin
    # initialize model and references
    m = InfiniteModel()
    @infinite_parameter(m, par in [0, 1])
    @variable(m, inf >= 0., Infinite(par))
    @variable(m, pt, Point(inf, 0.5))
    @variable(m, 0. <= x <= 1., Int)
    expected = [(GeneralVariableRef, MOI.GreaterThan{Float64}),
                (GeneralVariableRef, MOI.LessThan{Float64}),
                (GeneralVariableRef, MOI.Integer)]
    @test isempty(setdiff(list_of_constraint_types(m), expected))
end
