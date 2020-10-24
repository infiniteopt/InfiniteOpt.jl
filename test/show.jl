using JuMP: REPLMode, IJuliaMode

# Test string creation
@testset "String Creators" begin
    # initialize model and attributes
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par1 <= 1)
    @infinite_parameter(m, pars[1:2] in MvNormal([1, 1], 1))
    @infinite_parameter(m, pars2[1:2] in [0, 2])
    @infinite_parameter(m, pars3[1:2] in [0, 1], independent = true)
    @infinite_variable(m, x(par1))
    @infinite_variable(m, z(pars))
    @infinite_variable(m, inf(pars, par1, pars3))
    @hold_variable(m, y)
    d1 = @deriv(x, par1)
    d2 = @deriv(inf, pars[1], par1)
    d3 = @deriv(z, pars[2])
    set_name(d3, "d3")
    @objective(m, Min, 2 + y)
    @constraint(m, c1, x + y - 2 <= 0)
    ac1 = @constraint(m, x + y - 2 <= 0)
    @constraint(m, c2, y^2 - 3 == 0)
    ac2 = @constraint(m, y^2 - 3 == 0)
    @BDconstraint(m, c3(par1 in [0, 0.5]), x == 5)
    ac3 = @BDconstraint(m, (par1 in [0, 0.5]), x == 5)
    # test _infopt_math_symbol (REPL)
    @testset "_infopt_math_symbol (REPL)" begin
        if Sys.iswindows()
            @test InfiniteOpt._infopt_math_symbol(REPLMode, :intersect) == "and"
            @test InfiniteOpt._infopt_math_symbol(REPLMode, :partial) == "d"
            @test InfiniteOpt._infopt_math_symbol(REPLMode, :integral) == "integral"
            @test InfiniteOpt._infopt_math_symbol(REPLMode, :expect) == "E"
        else
            @test InfiniteOpt._infopt_math_symbol(REPLMode, :intersect) == "∩"
            @test InfiniteOpt._infopt_math_symbol(REPLMode, :partial) == "∂"
            @test InfiniteOpt._infopt_math_symbol(REPLMode, :integral) == "∫"
            @test InfiniteOpt._infopt_math_symbol(REPLMode, :expect) == "𝔼"
        end
        @test InfiniteOpt._infopt_math_symbol(REPLMode, :times) == "*"
        @test InfiniteOpt._infopt_math_symbol(REPLMode, :prop) == "~"
    end
    # test _infopt_math_symbol (IJulia)
    @testset "_infopt_math_symbol (IJulia)" begin
        @test InfiniteOpt._infopt_math_symbol(IJuliaMode, :intersect) == "\\cap"
        @test InfiniteOpt._infopt_math_symbol(IJuliaMode, :eq) == "="
        @test InfiniteOpt._infopt_math_symbol(IJuliaMode, :prop) == "\\sim"
        @test InfiniteOpt._infopt_math_symbol(IJuliaMode, :partial) == "\\partial"
        @test InfiniteOpt._infopt_math_symbol(IJuliaMode, :open_rng) == "\\left["
        @test InfiniteOpt._infopt_math_symbol(IJuliaMode, :close_rng) == "\\right]"
        @test InfiniteOpt._infopt_math_symbol(IJuliaMode, :integral) == "\\int"
        @test InfiniteOpt._infopt_math_symbol(IJuliaMode, :expect) == "\\mathbb{E}"
    end
    # test _plural
    @testset "_plural" begin
        @test InfiniteOpt._plural(1) == ""
        @test InfiniteOpt._plural(2) == "s"
    end
    # test set_string (IntervalSet)
    @testset "set_string (IntervalSet)" begin
        # test simple case
        set = IntervalSet(0, 1)
        @test InfiniteOpt.set_string(REPLMode, set) == "[0, 1]"
        @test InfiniteOpt.set_string(IJuliaMode, set) == "[0, 1]"
        # test rounding case
        set = IntervalSet(-0, 1)
        @test InfiniteOpt.set_string(REPLMode, set) == "[0, 1]"
        @test InfiniteOpt.set_string(IJuliaMode, set) == "[0, 1]"
        # test decimal case
        set = IntervalSet(0.1, 1.3)
        @test InfiniteOpt.set_string(REPLMode, set) == "[0.1, 1.3]"
        @test InfiniteOpt.set_string(IJuliaMode, set) == "[0.1, 1.3]"
    end
    # test set_string (DistributionSet)
    @testset "set_string (DistributionSet)" begin
        # test univariate set
        set = UniDistributionSet(Uniform())
        @test InfiniteOpt.set_string(REPLMode, set) == "Uniform{Float64}(a=0.0, b=1.0)"
        @test InfiniteOpt.set_string(IJuliaMode, set) == "Uniform{Float64}(a=0.0, b=1.0)"
        # test mulivariate set
        set = MultiDistributionSet(MvNormal([1], 1))
        str = "IsoNormal(\ndim: 1\nμ: [1.0]\nΣ: [1.0]\n)\n"
        @test InfiniteOpt.set_string(REPLMode, set) == str
        @test InfiniteOpt.set_string(IJuliaMode, set) == str
    end
    # test set_string (CollectionSet)
    @testset "set_string (CollectionSet)" begin
        set = CollectionSet([IntervalSet(0, 1), IntervalSet(0, 0.1)])
        str = "CollectionSet with 2 sets:\n [0, 1]\n [0, 0.1]"
        @test InfiniteOpt.set_string(REPLMode, set) == str
        @test InfiniteOpt.set_string(IJuliaMode, set) == str
    end
    # test set_string (Fallback)
    @testset "set_string (Fallback)" begin
        set = BadSet()
        @test InfiniteOpt.set_string(REPLMode, set) == "BadSet()"
        @test InfiniteOpt.set_string(IJuliaMode, set) == "BadSet()"
    end
    # test in_set_string (IntervalSet)
    @testset "JuMP.in_set_string (Interval)" begin
        # test simple case
        set = IntervalSet(0, 1)
        str = JuMP._math_symbol(REPLMode, :in) * " [0, 1]"
        @test in_set_string(REPLMode, set) == str
        str = JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]"
        @test in_set_string(IJuliaMode, set) == str
        # test rounding case
        set = IntervalSet(-0, 1)
        str = JuMP._math_symbol(REPLMode, :in) * " [0, 1]"
        @test in_set_string(REPLMode, set) == str
        str = JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]"
        @test in_set_string(IJuliaMode, set) == str
        # test decimal case
        set = IntervalSet(0.1, 1.3)
        str = JuMP._math_symbol(REPLMode, :in) * " [0.1, 1.3]"
        @test in_set_string(REPLMode, set) == str
        str = JuMP._math_symbol(IJuliaMode, :in) * " [0.1, 1.3]"
        @test in_set_string(IJuliaMode, set) == str
        # test finite case
        set = IntervalSet(0.1, 0.1)
        str = JuMP._math_symbol(REPLMode, :eq) * " 0.1"
        @test in_set_string(REPLMode, set) == str
        str = JuMP._math_symbol(IJuliaMode, :eq) * " 0.1"
        @test in_set_string(IJuliaMode, set) == str
    end
    # test in_set_string (Distribution)
    @testset "JuMP.in_set_string (Distribution)" begin
        # test univariate set
        set = UniDistributionSet(Uniform())
        str = InfiniteOpt._infopt_math_symbol(REPLMode, :prop) * " Uniform"
        @test in_set_string(REPLMode, set) == str
        str = InfiniteOpt._infopt_math_symbol(IJuliaMode, :prop) * " Uniform"
        @test in_set_string(IJuliaMode, set) == str
        # test mulivariate set
        set = MultiDistributionSet(MvNormal([1], 1))
        str = InfiniteOpt._infopt_math_symbol(REPLMode, :prop) * " MvNormal(dim: (1))"
        @test in_set_string(REPLMode, set) == str
        str = InfiniteOpt._infopt_math_symbol(IJuliaMode, :prop) * " MvNormal(dim: (1))"
        @test in_set_string(IJuliaMode, set) == str
        # test matrix set
        set = MultiDistributionSet(MatrixBeta(2, 2, 2))
        str = InfiniteOpt._infopt_math_symbol(REPLMode, :prop) * " MatrixBeta(dims: (2, 2))"
        @test in_set_string(REPLMode, set) == str
        str = InfiniteOpt._infopt_math_symbol(IJuliaMode, :prop) * " MatrixBeta(dims: (2, 2))"
        @test in_set_string(IJuliaMode, set) == str
    end
    # test in_set_string (Fallback)
    @testset "JuMP.in_set_string (Fallback)" begin
        set = BadSet()
        in1 = JuMP._math_symbol(REPLMode, :in)
        in2 = JuMP._math_symbol(IJuliaMode, :in)
        @test in_set_string(REPLMode, set) == in1 * " BadSet()"
        @test in_set_string(IJuliaMode, set) == in2 * " BadSet()"
    end
    # test in_set_string (IntervalSet with Bounds)
    @testset "JuMP.in_set_string (IntervalSet with Bounds)" begin
        # test in bounds
        bounds = ParameterBounds((par1 => IntervalSet(0, 0),))
        set = IntervalSet(0, 1)
        str = JuMP._math_symbol(REPLMode, :eq) * " 0"
        @test in_set_string(REPLMode, par1, set, bounds) == str
        str = JuMP._math_symbol(IJuliaMode, :eq) * " 0"
        @test in_set_string(IJuliaMode, par1, set, bounds) == str
        # test not in bounds
        str = JuMP._math_symbol(REPLMode, :in) * " [0, 1]"
        @test in_set_string(REPLMode, pars[1], set, bounds) == str
        str = JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]"
        @test in_set_string(IJuliaMode, pars[1], set, bounds) == str
    end
    # test in_set_string (InfiniteScalarSet with Bounds)
    @testset "JuMP.in_set_string (InfiniteScalarSet with Bounds)" begin
        # test in bounds
        bounds = ParameterBounds((par1 => IntervalSet(0, 0),))
        set = UniDistributionSet(Uniform())
        str = JuMP._math_symbol(REPLMode, :eq) * " 0"
        @test in_set_string(REPLMode, par1, set, bounds) == str
        str = JuMP._math_symbol(IJuliaMode, :eq) * " 0"
        @test in_set_string(IJuliaMode, par1, set, bounds) == str
        # test in bounds and not equality
        bounds = ParameterBounds((par1 => IntervalSet(0, 1),))
        set = UniDistributionSet(Uniform())
        str = InfiniteOpt._infopt_math_symbol(REPLMode, :prop) * " Uniform " *
              InfiniteOpt._infopt_math_symbol(REPLMode, :intersect) * " [0, 1]"
        @test in_set_string(REPLMode, par1, set, bounds) == str
        str = InfiniteOpt._infopt_math_symbol(IJuliaMode, :prop) * " Uniform " *
              InfiniteOpt._infopt_math_symbol(IJuliaMode, :intersect) * " [0, 1]"
        @test in_set_string(IJuliaMode, par1, set, bounds) == str
        # test not in bounds
        str = InfiniteOpt._infopt_math_symbol(REPLMode, :prop) * " Uniform"
        @test in_set_string(REPLMode, pars[1], set, bounds) == str
        str = InfiniteOpt._infopt_math_symbol(IJuliaMode, :prop) * " Uniform"
        @test in_set_string(IJuliaMode, pars[1], set, bounds) == str
    end
    # test measure_data_string with 1-D DiscreteMeasureData/FunctionalDiscreteMeasureData
    @testset "measure_data_string (1-D)" begin
        # test with bounds
        data = FunctionalDiscreteMeasureData(par1, ones, 0, All, default_weight, 0, 1, false)
        str = "par1 " * JuMP._math_symbol(REPLMode, :in) * " [0, 1]"
        @test InfiniteOpt.measure_data_string(REPLMode, data) == str
        str = "par1 " * JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]"
        @test InfiniteOpt.measure_data_string(IJuliaMode, data) == str
        # test without bounds
        data = FunctionalDiscreteMeasureData(par1, ones, 0, All, default_weight, NaN, NaN, false)
        @test InfiniteOpt.measure_data_string(REPLMode, data) == "par1"
        @test InfiniteOpt.measure_data_string(IJuliaMode, data) == "par1"
    end
    # test measure_data_string with Multi-D DiscreteMeasureData/FunctionalDiscreteMeasureData
    @testset "measure_data_string (Multi-D)" begin
        # test with homogenous bounds
        data = FunctionalDiscreteMeasureData(pars2, ones, 0, All, default_weight, [0, 0], [1, 1], false)
        str = "pars2 " * JuMP._math_symbol(REPLMode, :in) * " [0, 1]^2"
        @test InfiniteOpt.measure_data_string(REPLMode, data) == str
        str = "pars2 " * JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]^2"
        @test InfiniteOpt.measure_data_string(IJuliaMode, data) == str
        # test heterogeneous bounds
        data = FunctionalDiscreteMeasureData(pars2, ones, 0, All, default_weight, [0, 0], [0.5, 1], false)
        str = "pars2[1] " * JuMP._math_symbol(REPLMode, :in) * " [0, 0.5], " *
              "pars2[2] " * JuMP._math_symbol(REPLMode, :in) * " [0, 1]"
        @test InfiniteOpt.measure_data_string(REPLMode, data) == str
        str = "pars2_{1} " * JuMP._math_symbol(IJuliaMode, :in) * " [0, 0.5], " *
              "pars2_{2} " * JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]"
        @test InfiniteOpt.measure_data_string(IJuliaMode, data) == str
        # test no bounds with homogenous names
        data = FunctionalDiscreteMeasureData(pars2, ones, 0, All, default_weight, [NaN, NaN], [NaN, NaN], false)
        @test InfiniteOpt.measure_data_string(REPLMode, data) == "pars2"
        @test InfiniteOpt.measure_data_string(IJuliaMode, data) == "pars2"
        # test no bounds with heterogeneous names
        data = FunctionalDiscreteMeasureData([par1, pars[1]], ones, 0, All, default_weight, [NaN, NaN], [NaN, NaN], false)
        @test InfiniteOpt.measure_data_string(REPLMode, data) == "[par1, pars[1]]"
        @test InfiniteOpt.measure_data_string(IJuliaMode, data) == "[par1, pars[1]]"
    end
    # test _get_root_parameter_name
    @testset "_get_root_parameter_name" begin
        # test single
        data = TestData(par1, 0, 1)
        @test InfiniteOpt._get_root_parameter_name(data) == "par1"
        @test InfiniteOpt._get_root_parameter_name(data) == "par1"
        # test with homogenous names
        data = TestData(pars2, 0, 1)
        @test InfiniteOpt._get_root_parameter_name(data) == "pars2"
        @test InfiniteOpt._get_root_parameter_name(data) == "pars2"
        # test with heterogeneous names
        data = TestData([par1, pars[1]], 0, 1)
        @test InfiniteOpt._get_root_parameter_name(data) == "[par1, pars[1]]"
        @test InfiniteOpt._get_root_parameter_name(data) == "[par1, pars[1]]"
    end
    # test measure_data_string (Fallback)
    @testset "measure_data_string (Fallback)" begin
        # test single
        data = TestData(par1, 0, 1)
        @test InfiniteOpt.measure_data_string(REPLMode, data) == "par1"
        @test InfiniteOpt.measure_data_string(IJuliaMode, data) == "par1"
        # test with homogenous names
        data = TestData(pars2, 0, 1)
        @test InfiniteOpt.measure_data_string(REPLMode, data) == "pars2"
        @test InfiniteOpt.measure_data_string(IJuliaMode, data) == "pars2"
        # test with heterogeneous names
        data = TestData([par1, pars[1]], 0, 1)
        @test InfiniteOpt.measure_data_string(REPLMode, data) == "[par1, pars[1]]"
        @test InfiniteOpt.measure_data_string(IJuliaMode, data) == "[par1, pars[1]]"
    end
    # test variabel_string (MeasureRef)
    @testset "variable_string (MeasureRef)" begin
        # test non measure toolbox measures
        data = FunctionalDiscreteMeasureData(pars2, ones, 0, All, default_weight, [0, 0], [1, 1], false)
        meas = dispatch_variable_ref(measure(y, data, name = "test"))
        str = "test{pars2 " * JuMP._math_symbol(REPLMode, :in) * " [0, 1]^2}[y]"
        @test InfiniteOpt.variable_string(REPLMode, meas) == str
        str = "\\text{test}_{pars2 " * JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]^2}\\left[y\\right]"
        @test InfiniteOpt.variable_string(IJuliaMode, meas) == str
        # test measure toolbox special cases for integrals and expectations
        @infinite_parameter(m, t in [0, 1])
        meas = dispatch_variable_ref(expect(y, t))
        str = "\\mathbb{E}_{t}\\left[y\\right]"
        @test InfiniteOpt.variable_string(IJuliaMode, meas) == str
        str = InfiniteOpt._infopt_math_symbol(REPLMode, :expect) * "{t}[y]"
        @test InfiniteOpt.variable_string(REPLMode, meas) == str
        meas = dispatch_variable_ref(integral(y, t))
        str = "\\int_{t " * JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]}ydt"
        @test InfiniteOpt.variable_string(IJuliaMode, meas) == str
        int = InfiniteOpt._infopt_math_symbol(REPLMode, :integral)
        str = int * "{t " * JuMP._math_symbol(REPLMode, :in) * " [0, 1]}[y]"
        @test InfiniteOpt.variable_string(REPLMode, meas) == str
        meas = dispatch_variable_ref(integral(y, par1))
        str = "\\int_{par1 " * JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]}yd(par1)"
        @test InfiniteOpt.variable_string(IJuliaMode, meas) == str
    end
    # test _get_base_name
    @testset "_get_base_name (REPL)" begin
        @test InfiniteOpt._get_base_name(REPLMode, pars[1]) == "pars[1]"
        bad_ref = HoldVariableRef(m, HoldVariableIndex(-1))
        @test InfiniteOpt._get_base_name(REPLMode, bad_ref) == "noname"
    end
    # test _get_base_name
    @testset "_get_base_name (IJulia)" begin
        @test InfiniteOpt._get_base_name(IJuliaMode, pars[1]) == "pars_{1}"
        bad_ref = HoldVariableRef(m, HoldVariableIndex(-1))
        @test InfiniteOpt._get_base_name(IJuliaMode, bad_ref) == "noname"
    end
    # test variable_string (InfiniteVariableRef)
    @testset "variable_string (InfiniteVariableRef)" begin
        # test bad ref
        bad_ref = InfiniteVariableRef(m, InfiniteVariableIndex(-1))
        @test InfiniteOpt.variable_string(REPLMode, bad_ref) == "noname"
        @test InfiniteOpt.variable_string(IJuliaMode, bad_ref) == "noname"
        # test complex normal one
        dvref = dispatch_variable_ref(inf)
        @test InfiniteOpt.variable_string(REPLMode, dvref) == "inf(pars, par1, pars3)"
        @test InfiniteOpt.variable_string(IJuliaMode, dvref) == "inf(pars, par1, pars3)"
        # test mixed types
        @infinite_variable(m, inf2([par1, pars3[2]]))
        dvref = dispatch_variable_ref(inf2)
        @test InfiniteOpt.variable_string(REPLMode, dvref) == "inf2([par1, pars3[2]])"
        @test InfiniteOpt.variable_string(IJuliaMode, dvref) == "inf2([par1, pars3[2]])"
    end
    # test variable_string (DerivativeRef)
    @testset "variable_string (DerivativeRef)" begin
        # test bad ref
        bad_ref = DerivativeRef(m, DerivativeIndex(-1))
        @test InfiniteOpt.variable_string(REPLMode, bad_ref) == "noname"
        @test InfiniteOpt.variable_string(IJuliaMode, bad_ref) == "noname"
        # test normal one
        dref = dispatch_variable_ref(d1)
        d_re = InfiniteOpt._infopt_math_symbol(REPLMode, :partial)
        @test InfiniteOpt.variable_string(REPLMode, dref) == "$d_re/$(d_re)par1[x(par1)]"
        @test InfiniteOpt.variable_string(IJuliaMode, dref) == "\\frac{\\partial}{\\partial par1}\\left[x(par1)\\right]"
        # test nested one 
        dref = dispatch_variable_ref(d2)
        expected = "$d_re/$(d_re)par1[$d_re/$(d_re)pars[1][inf(pars, par1, pars3)]]"
        @test InfiniteOpt.variable_string(REPLMode, dref) == expected
        expected = "\\frac{\\partial}{\\partial par1}\\left[\\frac{\\partial}{\\partial pars_{1}}\\left[inf(pars, par1, pars3)\\right]\\right]"
        @test InfiniteOpt.variable_string(IJuliaMode, dref) == expected
        # test has explicit name
        dref = dispatch_variable_ref(d3)
        @test InfiniteOpt.variable_string(REPLMode, dref) == "d3(pars)"
        @test InfiniteOpt.variable_string(IJuliaMode, dref) == "d3(pars)"
    end
    # _make_str_value (Number)
    @testset "_make_str_value (Number)" begin
        @test InfiniteOpt._make_str_value(1.0) == "1"
    end
    # _make_str_value (Array)
    @testset "_make_str_value (Array)" begin
        # test single
        @test InfiniteOpt._make_str_value([1.1]) == "1.1"
        # test short array
        values = [1., 2., 3.]
        @test InfiniteOpt._make_str_value(values) == "[1, 2, 3]"
        # test long array
        values = [1., 2., 3., 4., 5., 6.]
        @test InfiniteOpt._make_str_value(values) == "[1, ..., 6]"
    end
    # test variable_string (PointVariableRef)
    @testset "variable_string (PointVariableRef)" begin
        # test bad ref
        bad_ref = PointVariableRef(m, PointVariableIndex(-1))
        @test InfiniteOpt.variable_string(REPLMode, bad_ref) == "noname"
        @test InfiniteOpt.variable_string(IJuliaMode, bad_ref) == "noname"
        # test short one
        dvref = dispatch_variable_ref(@point_variable(m, x(0)))
        @test InfiniteOpt.variable_string(REPLMode, dvref) == "x(0)"
        @test InfiniteOpt.variable_string(IJuliaMode, dvref) == "x(0)"
        # test complex one
        dvref = dispatch_variable_ref(@point_variable(m, inf([0, 0], 0, [0, 0])))
        @test InfiniteOpt.variable_string(REPLMode, dvref) == "inf([0, 0], 0, [0, 0])"
        @test InfiniteOpt.variable_string(IJuliaMode, dvref) == "inf([0, 0], 0, [0, 0])"
        # test with alias name
        dvref = dispatch_variable_ref(@point_variable(m, z([0, 0]), z0))
        @test InfiniteOpt.variable_string(REPLMode, dvref) == "z0"
        @test InfiniteOpt.variable_string(IJuliaMode, dvref) == "z0"
        # test with derivative 
        dvref = dispatch_variable_ref(@point_variable(m, d1(0)))
        d_re = InfiniteOpt._infopt_math_symbol(REPLMode, :partial)
        @test InfiniteOpt.variable_string(REPLMode, dvref) == "$(d_re)/$(d_re)par1[x(par1)](0)"
        # test named derivative 
        dvref = dispatch_variable_ref(@point_variable(m, d3([0, 0])))
        @test InfiniteOpt.variable_string(REPLMode, dvref) == "d3([0, 0])"
    end
    # test variable_string (ReducedVariableRef)
    @testset "variable_string (ReducedVariableRef)" begin
        # test bad ref
        bad_ref = ReducedVariableRef(m, ReducedVariableIndex(-1))
        @test InfiniteOpt.variable_string(REPLMode, bad_ref) == "noname"
        @test InfiniteOpt.variable_string(IJuliaMode, bad_ref) == "noname"
        # test short one
        eval_supps = Dict{Int, Float64}(1 => 0, 2 => 0, 4=>0)
        var = build_variable(error, inf, eval_supps, check = false)
        dvref = dispatch_variable_ref(add_variable(m, var))
        @test InfiniteOpt.variable_string(REPLMode, dvref) == "inf([0, 0], par1, [0, pars3[2]])"
        @test InfiniteOpt.variable_string(IJuliaMode, dvref) == "inf([0, 0], par1, [0, pars3[2]])"
        # test named derivative 
        eval_supps = Dict{Int, Float64}(1 => 0)
        var = build_variable(error, d3, eval_supps, check = false)
        dvref = dispatch_variable_ref(add_variable(m, var))
        @test InfiniteOpt.variable_string(REPLMode, dvref) == "d3([0, pars[2]])"
        # unnamed derivative
        eval_supps = Dict{Int, Float64}(1 => 0)
        var = build_variable(error, d1, eval_supps, check = false)
        dvref = dispatch_variable_ref(add_variable(m, var))
        d_re = InfiniteOpt._infopt_math_symbol(REPLMode, :partial)
        @test InfiniteOpt.variable_string(REPLMode, dvref) == "$(d_re)/$(d_re)par1[x(par1)](0)"
    end
    # test variable_string (Fallback)
    @testset "variable_string (Fallback)" begin
        @test InfiniteOpt.variable_string(REPLMode, y) == "y"
        @test InfiniteOpt.variable_string(IJuliaMode, y) == "y"
    end
    # test the function_string extensions
    @testset "JuMP.function_string" begin
        @test JuMP.function_string(REPLMode, dispatch_variable_ref(y)) == "y"
        @test JuMP.function_string(IJuliaMode, dispatch_variable_ref(inf)) == "inf(pars, par1, pars3)"
        @test JuMP.function_string(REPLMode, y) == "y"
        @test JuMP.function_string(IJuliaMode, inf) == "inf(pars, par1, pars3)"
        d_re = InfiniteOpt._infopt_math_symbol(REPLMode, :partial)
        @test JuMP.function_string(REPLMode, d1) == "$d_re/$(d_re)par1[x(par1)]"
    end
    # test bound_string
    @testset "bound_string" begin
        # test with single bound
        bounds = ParameterBounds(Dict(par1 => IntervalSet(0.5, 0.7)))
        str = "par1 " * JuMP._math_symbol(REPLMode, :in) * " [0.5, 0.7]"
        @test InfiniteOpt.bound_string(REPLMode, bounds) == str
        str = "par1 " *  JuMP._math_symbol(IJuliaMode, :in) * " [0.5, 0.7]"
        @test InfiniteOpt.bound_string(IJuliaMode, bounds) == str
    end
    # test constraint_string (Finite constraint)
    @testset "JuMP.constraint_string (Finite)" begin
        # test named
        str = "c2 : y² " * JuMP._math_symbol(REPLMode, :eq) * " 3.0"
        @test constraint_string(REPLMode, c2) == str
        str =  "c2 : \$ y^2 " * JuMP._math_symbol(IJuliaMode, :eq) * " 3.0 \$"
        @test constraint_string(IJuliaMode, c2) == str
        # test unnamed
        str = "y² " * JuMP._math_symbol(REPLMode, :eq) * " 3.0"
        @test constraint_string(REPLMode, ac2) == str
        str =  "\$ y^2 " * JuMP._math_symbol(IJuliaMode, :eq) * " 3.0 \$"
        @test constraint_string(IJuliaMode, ac2) == str
        # test named in math mode
        str = "c2 : y² " * JuMP._math_symbol(REPLMode, :eq) * " 3.0"
        @test constraint_string(REPLMode, c2, in_math_mode = true) == str
        str =  "y^2 " * JuMP._math_symbol(IJuliaMode, :eq) * " 3.0"
        @test constraint_string(IJuliaMode, c2, in_math_mode = true) == str
    end
    # test _param_domain_string (IndependentParameter)
    @testset "_param_domain_string (IndependentParameter)" begin
        bounds = ParameterBounds((par1 => IntervalSet(0, 0),))
        idx = index(par1)
        str = "par1 " * JuMP._math_symbol(REPLMode, :eq) * " 0"
        @test InfiniteOpt._param_domain_string(REPLMode, m, idx, bounds) == str
        str = "par1 " * JuMP._math_symbol(IJuliaMode, :eq) * " 0"
        @test InfiniteOpt._param_domain_string(IJuliaMode, m, idx, bounds) == str
    end
    # test _param_domain_string (DependentParameters)
    @testset "_param_domain_string (DependentParameters)" begin
        # Collection set
        bounds = ParameterBounds((pars2[1] => IntervalSet(0, 1),))
        idx = index(pars2[1]).object_index
        str = "pars2[1] " * JuMP._math_symbol(REPLMode, :in) * " [0, 1], " *
              "pars2[2] " * JuMP._math_symbol(REPLMode, :in) * " [0, 2]"
        @test InfiniteOpt._param_domain_string(REPLMode, m, idx, bounds) == str
        str = "pars2_{1} " * JuMP._math_symbol(IJuliaMode, :in) * " [0, 1], " *
              "pars2_{2} " * JuMP._math_symbol(IJuliaMode, :in) * " [0, 2]"
        @test InfiniteOpt._param_domain_string(IJuliaMode, m, idx, bounds) == str
        # other set with equalities
        bounds = ParameterBounds((pars[1] => IntervalSet(0, 0),
                                  pars[2] => IntervalSet(1, 1)))
        idx = index(pars[1]).object_index
        str = InfiniteOpt.bound_string(REPLMode, bounds)
        str2 = string(split(str, ", ")[2], ", ", split(str, ", ")[1])
        @test InfiniteOpt._param_domain_string(REPLMode, m, idx, bounds) in [str, str2]
        str = InfiniteOpt.bound_string(IJuliaMode, bounds)
        str2 = string(split(str, ", ")[2], ", ", split(str, ", ")[1])
        @test InfiniteOpt._param_domain_string(IJuliaMode, m, idx, bounds) in [str, str2]
        # other set without equalities and including in the bounds
        bounds = ParameterBounds((pars[1] => IntervalSet(0, 1),))
        str = "pars " * InfiniteOpt._infopt_math_symbol(REPLMode, :prop) *
              " MvNormal(dim: (2)) " *
              InfiniteOpt._infopt_math_symbol(REPLMode, :intersect) *
              " (pars[1] " * JuMP._math_symbol(REPLMode, :in) * " [0, 1])"
        @test InfiniteOpt._param_domain_string(REPLMode, m, idx, bounds) == str
        str = "pars " * InfiniteOpt._infopt_math_symbol(IJuliaMode, :prop) *
              " MvNormal(dim: (2)) " *
              InfiniteOpt._infopt_math_symbol(IJuliaMode, :intersect) *
              " (pars_{1} " * JuMP._math_symbol(IJuliaMode, :in) * " [0, 1])"
        @test InfiniteOpt._param_domain_string(IJuliaMode, m, idx, bounds) == str
        # other set without equalities and not included in bounds
        bounds = ParameterBounds((par1 => IntervalSet(0, 1),))
        str = "pars " * InfiniteOpt._infopt_math_symbol(REPLMode, :prop) *
              " MvNormal(dim: (2))"
        @test InfiniteOpt._param_domain_string(REPLMode, m, idx, bounds) == str
        str = "pars " * InfiniteOpt._infopt_math_symbol(IJuliaMode, :prop) *
              " MvNormal(dim: (2))"
        @test InfiniteOpt._param_domain_string(IJuliaMode, m, idx, bounds) == str
    end
    # test constraint_string (infinite constraint)
    @testset "JuMP.constraint_string (Infinite)" begin
        # test c1 with name
        str = "c1 : x(par1) + y " * JuMP._math_symbol(REPLMode, :leq) * " 2.0, " *
              JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
              JuMP._math_symbol(REPLMode, :in) * " [0, 1]"
        @test constraint_string(REPLMode, c1) == str
        str = "c1 : \$ x(par1) + y " * JuMP._math_symbol(IJuliaMode, :leq) * " 2.0, " *
              JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
              JuMP._math_symbol(IJuliaMode, :in) * " [0, 1] \$"
        @test constraint_string(IJuliaMode, c1) == str
        # test c1 without name
        str = "x(par1) + y " * JuMP._math_symbol(REPLMode, :leq) * " 2.0, " *
              JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
              JuMP._math_symbol(REPLMode, :in) * " [0, 1]"
        @test constraint_string(REPLMode, ac1) == str
        str = "\$ x(par1) + y " * JuMP._math_symbol(IJuliaMode, :leq) * " 2.0, " *
              JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
              JuMP._math_symbol(IJuliaMode, :in) * " [0, 1] \$"
        @test constraint_string(IJuliaMode, ac1) == str
        # test c1 with name and in_math_mode
        str = "c1 : x(par1) + y " * JuMP._math_symbol(REPLMode, :leq) * " 2.0, " *
              JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
              JuMP._math_symbol(REPLMode, :in) * " [0, 1]"
        @test constraint_string(REPLMode, c1, in_math_mode = true) == str
        str = "x(par1) + y " * JuMP._math_symbol(IJuliaMode, :leq) * " 2.0, " *
              JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
              JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]"
        @test constraint_string(IJuliaMode, c1, in_math_mode = true) == str
        # test c3 with name
        str = "c3 : x(par1) " * JuMP._math_symbol(REPLMode, :eq) * " 5.0, " *
              JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
              JuMP._math_symbol(REPLMode, :in) * " [0, 0.5]"
        @test constraint_string(REPLMode, c3) == str
        str =  "c3 : \$ x(par1) " * JuMP._math_symbol(IJuliaMode, :eq) * " 5.0, " *
               JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
               JuMP._math_symbol(IJuliaMode, :in) * " [0, 0.5] \$"
        @test constraint_string(IJuliaMode, c3) == str
        # test c3 without name
        str = "x(par1) " * JuMP._math_symbol(REPLMode, :eq) * " 5.0, " *
              JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
              JuMP._math_symbol(REPLMode, :in) * " [0, 0.5]"
        @test constraint_string(REPLMode, ac3) == str
        str =  "\$ x(par1) " * JuMP._math_symbol(IJuliaMode, :eq) * " 5.0, " *
               JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
               JuMP._math_symbol(IJuliaMode, :in) * " [0, 0.5] \$"
        @test constraint_string(IJuliaMode, ac3) == str
    end
    # test constraints_string
    @testset "JuMP.constraints_string" begin
        # test REPLMode
        strings = Vector{String}(undef, 6)
        strings[1] = "c1 : x(par1) + y " * JuMP._math_symbol(REPLMode, :leq) *
                     " 2.0, " * JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
                     JuMP._math_symbol(REPLMode, :in) * " [0, 1]"
        strings[2] = "x(par1) + y " * JuMP._math_symbol(REPLMode, :leq) *
                     " 2.0, " * JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
                     JuMP._math_symbol(REPLMode, :in) * " [0, 1]"
        strings[3] = "c2 : y² " * JuMP._math_symbol(REPLMode, :eq) * " 3.0"
        strings[4] = "y² " * JuMP._math_symbol(REPLMode, :eq) * " 3.0"
        strings[5] = "c3 : x(par1) " * JuMP._math_symbol(REPLMode, :eq) * " 5.0, " *
                     JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
                     JuMP._math_symbol(REPLMode, :in) * " [0, 0.5]"
        strings[6] = "x(par1) " * JuMP._math_symbol(REPLMode, :eq) * " 5.0, " *
                     JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
                     JuMP._math_symbol(REPLMode, :in) * " [0, 0.5]"
        @test constraints_string(REPLMode, m) == strings
        # test IJuliaMode
        strings = Vector{String}(undef, 6)
        strings[1] = "x(par1) + y " * JuMP._math_symbol(IJuliaMode, :leq) *
                     " 2.0, " * JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
                     JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]"
        strings[2] = "x(par1) + y " * JuMP._math_symbol(IJuliaMode, :leq) *
                     " 2.0, " * JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
                     JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]"
        strings[3] = "y^2 " * JuMP._math_symbol(IJuliaMode, :eq) * " 3.0"
        strings[4] = "y^2 " * JuMP._math_symbol(IJuliaMode, :eq) * " 3.0"
        strings[5] = "x(par1) " * JuMP._math_symbol(IJuliaMode, :eq) * " 5.0, " *
                     JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
                     JuMP._math_symbol(IJuliaMode, :in) * " [0, 0.5]"
        strings[6] = "x(par1) " * JuMP._math_symbol(IJuliaMode, :eq) * " 5.0, " *
                     JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
                     JuMP._math_symbol(IJuliaMode, :in) * " [0, 0.5]"
        @test constraints_string(IJuliaMode, m) == strings
    end
    # test objective_function_string
    @testset "JuMP.objective_function_string" begin
        @test objective_function_string(REPLMode, m) == "y + 2"
        @test objective_function_string(IJuliaMode, m) == "y + 2"
    end
end

# Helper function to test IO methods work correctly
function show_test(mode, obj, exp_str; repl=:both)
    if mode == REPLMode
        repl != :show  && @test sprint(print, obj) == exp_str
        repl != :print && @test sprint(show,  obj) == exp_str
    else
        @test sprint(show, "text/latex", obj) == exp_str
    end
end

# Make another helper function for other io methods
io_test(f::Function, exp_str::String, args...) = begin
    io = IOBuffer()
    f(io, args...)
    @test String(take!(io)) == exp_str
end

# test infinite model show
@testset "Show InfiniteModel" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par1 <= 1)
    @infinite_parameter(m, pars[1:2] in MvNormal([1, 1], 1))
    @infinite_variable(m, x(par1))
    @infinite_variable(m, z(pars))
    @hold_variable(m, y)
    @objective(m, Min, 2 + y)
    @constraint(m, c1, x + y -2 <= 0)
    bounds = ParameterBounds(Dict(par1 => IntervalSet(0.1, 1)))
    @BDconstraint(m, c3(par1 in [0, 0.5]), x <= 5)
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    # test Base.show (IntervalSet in REPL)
    @testset "Base.show (REPL IntervalSet)" begin
        show_test(REPLMode, IntervalSet(0, 1), "[0, 1]")
    end
    # test Base.show (IntervalSet in IJulia)
    @testset "Base.show (IJulia IntervalSet)" begin
        show_test(IJuliaMode, IntervalSet(0, 1), "[0, 1]")
    end
    # test Base.show (DistributionSet in REPL)
    @testset "Base.show (REPL DistributionSet)" begin
        show_test(REPLMode, UniDistributionSet(Uniform()), string(Uniform()))
    end
    # test Base.show (DistributionSet in IJulia)
    @testset "Base.show (IJulia DistributionSet)" begin
        show_test(IJuliaMode, UniDistributionSet(Uniform()), string(Uniform()))
    end
    # test Base.show (CollectionSet in REPL)
    @testset "Base.show (REPL CollectionSet)" begin
        show_test(REPLMode, CollectionSet([IntervalSet(0, 0)]),
                  "CollectionSet with 1 set:\n [0, 0]")
    end
    # test Base.show (CollectionSet in IJulia)
    @testset "Base.show (IJulia CollectionSet)" begin
        show_test(IJuliaMode, CollectionSet([IntervalSet(0, 0), IntervalSet(0, 2)]),
                  "CollectionSet with 2 sets:\n [0, 0]\n [0, 2]")
    end
    # test Base.show (ParameterBounds in REPL)
    @testset "Base.show (REPL ParameterBounds)" begin
        str = "Subdomain bounds (1): par1 " * JuMP._math_symbol(REPLMode, :in) *
              " [0.1, 1]"
        show_test(REPLMode, bounds, str)
    end
    # test Base.show (ParameterBounds in IJulia)
    @testset "Base.show (IJulia ParameterBounds)" begin
        str = "Subdomain bounds (1): par1 " * JuMP._math_symbol(IJuliaMode, :in) *
              " [0.1, 1]"
        show_test(IJuliaMode, bounds, str)
    end
    # test Base.show (GeneralVariableRef in IJulia)
    @testset "Base.show (IJulia GeneralVariableRef)" begin
        show_test(IJuliaMode, y, "\$\$ y \$\$")
    end
    # test Base.show (GeneralVariableRef in REPL)
    @testset "Base.show (REPL GeneralVariableRef)" begin
        show_test(REPLMode, y, "y")
    end
    # test Base.show (constraint in REPL)
    @testset "Base.show (REPL Constraint)" begin
        # test normal
        str = "c1 : x(par1) + y " * JuMP._math_symbol(REPLMode, :leq) *
              " 2.0, " * JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
              JuMP._math_symbol(REPLMode, :in) * " [0, 1]"
        show_test(REPLMode, c1, str)
        # test bounded
        str = "c3 : x(par1) " * JuMP._math_symbol(REPLMode, :leq) * " 5.0, " *
              JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
              JuMP._math_symbol(REPLMode, :in) * " [0, 0.5]"
        show_test(REPLMode, c3, str)
    end
    # test Base.show (constraint in IJulia)
    @testset "Base.show (IJulia Constraint)" begin
        # test normal
        str = "c1 : \$ x(par1) + y " * JuMP._math_symbol(IJuliaMode, :leq) *
              " 2.0, " * JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
              JuMP._math_symbol(IJuliaMode, :in) * " [0, 1] \$"
        show_test(IJuliaMode, c1, str)
        # test bounded
        str =  "c3 : \$ x(par1) " * JuMP._math_symbol(IJuliaMode, :leq) * " 5.0, " *
               JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
               JuMP._math_symbol(IJuliaMode, :in) * " [0, 0.5] \$"
        show_test(IJuliaMode, c3, str)
    end
    # test show_backend_summary
    @testset "JuMP.show_backend_summary" begin
        # test without optimizer
        str = "Optimizer model backend information: \nModel mode: AUTOMATIC\n" *
              "CachingOptimizer state: NO_OPTIMIZER\nSolver name: No optimizer" *
              " attached."
        io_test(show_backend_summary, str, m)
        # test with optimizer
        set_optimizer(optimizer_model(m), mockoptimizer)
        str = "Optimizer model backend information: \nModel mode: AUTOMATIC\n" *
              "CachingOptimizer state: EMPTY_OPTIMIZER\nSolver name: Mock"
        io_test(show_backend_summary, str, m)
    end
    # test show_objective_function_summary
    @testset "JuMP.show_objective_function_summary" begin
        str = "Objective function type: GenericAffExpr{Float64,GeneralVariableRef}\n"
        io_test(show_objective_function_summary, str, m)
    end
    # test show_constraints_summary
    @testset "JuMP.show_constraints_summary" begin
        # test the main function
        str = "`GenericAffExpr{Float64,GeneralVariableRef}`-in-`MathOptInter" *
              "face.LessThan{Float64}`: 2 constraints\n"
        io_test(show_constraints_summary, str, m)
    end
    # test show_objective_function_summary
    @testset "Base.show (InfiniteModel)" begin
        # test minimization
        str = "An InfiniteOpt Model\nMinimization problem with:\nFinite " *
              "Parameters: 0\nInfinite Parameters: 3\nVariables: 3" *
              "\nDerivatives: 0\nMeasures: 0" *
              "\nObjective function type: GenericAffExpr{Float64,General" *
              "VariableRef}\n`GenericAffExpr{Float64,GeneralVariableRef}`-in-" *
              "`MathOptInterface.LessThan{Float64}`: 2 constraints" *
              "\nNames registered in the model: c1, c3, par1, " *
              "pars, x, y, z\nOptimizer model backend information: \nModel " *
              "mode: AUTOMATIC\nCachingOptimizer state: EMPTY_OPTIMIZER\n" *
              "Solver name: Mock"
        show_test(REPLMode, m, str, repl=:show)
        # test maximization
        set_objective_sense(m, MOI.MAX_SENSE)
        str = "An InfiniteOpt Model\nMaximization problem with:\nFinite " *
              "Parameters: 0\nInfinite Parameters: 3\nVariables: 3" *
              "\nDerivatives: 0\nMeasures: 0" *
              "\nObjective function type: GenericAffExpr{Float64,General" *
              "VariableRef}\n`GenericAffExpr{Float64,GeneralVariableRef}`-in-" *
              "`MathOptInterface.LessThan{Float64}`: 2 constraints" *
              "\nNames registered in the model: c1, c3, par1, " *
              "pars, x, y, z\nOptimizer model backend information: \nModel " *
              "mode: AUTOMATIC\nCachingOptimizer state: EMPTY_OPTIMIZER\n" *
              "Solver name: Mock"
        show_test(REPLMode, m, str, repl=:show)
        # test feasibility
        set_objective_sense(m, MOI.FEASIBILITY_SENSE)
        str = "An InfiniteOpt Model\nFeasibility problem with:\nFinite " *
              "Parameters: 0\nInfinite Parameters: 3\nVariables: 3" *
              "\nDerivatives: 0\nMeasures: 0" *
              "\n`GenericAffExpr{Float64,GeneralVariableRef}`-in-`MathOpt" *
              "Interface.LessThan{Float64}`: 2 constraints" *
              "\nNames registered in the model: c1, c3, par1, " *
              "pars, x, y, z\nOptimizer model backend information: \nModel " *
              "mode: AUTOMATIC\nCachingOptimizer state: EMPTY_OPTIMIZER\n" *
              "Solver name: Mock"
        show_test(REPLMode, m, str, repl=:show)
    end
end
