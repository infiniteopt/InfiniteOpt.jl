# Test string creation
@testset "String Creators" begin
    # initialize model and attributes
    m = InfiniteModel()
    @infinite_parameter(m, par1 in [0, 1])
    @infinite_parameter(m, pars[1:2] ~ MvNormal([1, 1], [1 0; 0 1]))
    @infinite_parameter(m, pars2[1:2] in [0, 2])
    @infinite_parameter(m, pars3[1:2] in [0, 1], independent = true)
    @variable(m, x, Infinite(par1))
    @variable(m, z, Infinite(pars))
    @variable(m, inf, Infinite(pars, par1, pars3))
    @variable(m, y)
    d1 = @deriv(x, par1)
    d2 = @deriv(inf, pars[1], par1)
    d3 = @deriv(z, pars[2])
    set_name(d3, "d3")
    @objective(m, Min, 2 + y)
    @constraint(m, c1, x + y - 2 <= 0)
    ac1 = @constraint(m, x + y - 2 <= 0)
    @constraint(m, c2, y^2 - 3 == 0)
    ac2 = @constraint(m, y^2 - 3 == 0)
    @constraint(m, c3, x == 5, DomainRestrictions(par1 => [0, 0.5]))
    ac3 = @constraint(m, x == 5, DomainRestrictions(par1 => [0, 0.5]))
    # test _math_symbol (REPL)
    @testset "_math_symbol (REPL)" begin
        if Sys.iswindows()
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :intersect) == "and"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :partial) == "d"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :integral) == "integral"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :expect) == "E"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :leq) == "<="
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :geq) == ">="
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :eq) == "=="
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) == "for all"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :in) == "in"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :dots) == ".."
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :union) == "or"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :infty) == "Inf"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :Vert) == "||"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :sub2) == "_2"
        else
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :intersect) == "∩"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :partial) == "∂"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :integral) == "∫"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :expect) == "𝔼"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :leq) == "≤"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :geq) == "≥"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :eq) == "="
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) == "∀"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :in) == "∈"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :dots) == "…"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :union) == "∪"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :infty) == "∞"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :Vert) == "‖"
            @test InfiniteOpt._math_symbol(MIME("text/plain"), :sub2) == "₂"
        end
        @test InfiniteOpt._math_symbol(MIME("text/plain"), :times) == "*"
        @test InfiniteOpt._math_symbol(MIME("text/plain"), :prop) == "~"
        @test InfiniteOpt._math_symbol(MIME("text/plain"), :sq) == "²"
        @test InfiniteOpt._math_symbol(MIME("text/plain"), :ind_open) == "["
        @test InfiniteOpt._math_symbol(MIME("text/plain"), :ind_close) == "]"
        @test InfiniteOpt._math_symbol(MIME("text/plain"), :open_set) == "{"
        @test InfiniteOpt._math_symbol(MIME("text/plain"), :close_set) == "}"
        @test InfiniteOpt._math_symbol(MIME("text/plain"), :open_rng) == "["
        @test InfiniteOpt._math_symbol(MIME("text/plain"), :close_rng) == "]"
        @test InfiniteOpt._math_symbol(MIME("text/plain"), :integer) == "integer"
        @test InfiniteOpt._math_symbol(MIME("text/plain"), :succeq0) == " is semidefinite"
        @test_throws ErrorException InfiniteOpt._math_symbol(MIME("text/plain"), :bad)
    end
    # test _math_symbol (IJulia)
    @testset "_math_symbol (IJulia)" begin
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :intersect) == "\\cap"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :prop) == "\\sim"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :partial) == "\\partial"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :open_rng) == "\\left["
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :close_rng) == "\\right]"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :integral) == "\\int"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :expect) == "\\mathbb{E}"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :leq) == "\\leq"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :geq) == "\\geq"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :eq) == "="
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :times) == "\\times "
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :sq) == "^2"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :ind_open) == "_{"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :ind_close) == "}"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :for_all) == "\\quad\\forall"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :in) == "\\in"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :open_set) == "\\{"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :close_set) == "\\}"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :dots) == "\\dots"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :union) == "\\cup"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :infty) == "\\infty"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :integer) == "\\in \\mathbb{Z}"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :succeq0) == "\\succeq 0"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :Vert) == "\\Vert"
        @test InfiniteOpt._math_symbol(MIME("text/latex"), :sub2) == "_2"
        @test_throws ErrorException InfiniteOpt._math_symbol(MIME("text/latex"), :bad)
    end
    # test _plural
    @testset "_plural" begin
        @test InfiniteOpt._plural(1) == ""
        @test InfiniteOpt._plural(2) == "s"
    end
    # test _string_round
    @testset "_string_round" begin
        @test InfiniteOpt._string_round(1) == "1"
        @test InfiniteOpt._string_round(1.0) == "1"
        @test InfiniteOpt._string_round(0.0) == "0"
        @test InfiniteOpt._string_round(0.5) == "0.5"
    end
    # test domain_string (IntervalDomain)
    @testset "domain_string (IntervalDomain)" begin
        # test simple case
        domain = IntervalDomain(0, 1)
        @test InfiniteOpt.domain_string(MIME("text/plain"), domain) == "[0, 1]"
        @test InfiniteOpt.domain_string(MIME("text/latex"), domain) == "[0, 1]"
        # test rounding case
        domain = IntervalDomain(-0, 1)
        @test InfiniteOpt.domain_string(MIME("text/plain"), domain) == "[0, 1]"
        @test InfiniteOpt.domain_string(MIME("text/latex"), domain) == "[0, 1]"
        # test decimal case
        domain = IntervalDomain(0.1, 1.3)
        @test InfiniteOpt.domain_string(MIME("text/plain"), domain) == "[0.1, 1.3]"
        @test InfiniteOpt.domain_string(MIME("text/latex"), domain) == "[0.1, 1.3]"
    end
    # test domain_string (DistributionDomain)
    @testset "domain_string (DistributionDomain)" begin
        # test univariate domain
        domain = UniDistributionDomain(Uniform())
        @test InfiniteOpt.domain_string(MIME("text/plain"), domain) == "Uniform{Float64}(a=0.0, b=1.0)"
        @test InfiniteOpt.domain_string(MIME("text/latex"), domain) == "Uniform{Float64}(a=0.0, b=1.0)"
        # test mulivariate domain
        domain = MultiDistributionDomain(MvNormal([1], ones(1, 1)))
        str = "FullNormal(\ndim: 1" # just test first part
        @test InfiniteOpt.domain_string(MIME("text/plain"), domain)[1:length(str)] == str
        @test InfiniteOpt.domain_string(MIME("text/latex"), domain)[1:length(str)]  == str
    end
    # test domain_string (CollectionDomain)
    @testset "domain_string (CollectionDomain)" begin
        domain = CollectionDomain([IntervalDomain(0, 1), IntervalDomain(0, 0.1)])
        str = "CollectionDomain with 2 domains:\n [0, 1]\n [0, 0.1]"
        @test InfiniteOpt.domain_string(MIME("text/plain"), domain) == str
        @test InfiniteOpt.domain_string(MIME("text/latex"), domain) == str
    end
    # test domain_string (Fallback)
    @testset "domain_string (Fallback)" begin
        domain = BadDomain()
        @test InfiniteOpt.domain_string(MIME("text/plain"), domain) == "BadDomain()"
        @test InfiniteOpt.domain_string(MIME("text/latex"), domain) == "BadDomain()"
    end
    # test in_domain_string (IntervalDomain)
    @testset "in_domain_string (Interval)" begin
        # test simple case
        domain = IntervalDomain(0, 1)
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]"
        @test in_domain_string(MIME("text/plain"), domain) == str
        str = InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1]"
        @test in_domain_string(MIME("text/latex"), domain) == str
        # test rounding case
        domain = IntervalDomain(-0, 1)
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]"
        @test in_domain_string(MIME("text/plain"), domain) == str
        str = InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1]"
        @test in_domain_string(MIME("text/latex"), domain) == str
        # test decimal case
        domain = IntervalDomain(0.1, 1.3)
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0.1, 1.3]"
        @test in_domain_string(MIME("text/plain"), domain) == str
        str = InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0.1, 1.3]"
        @test in_domain_string(MIME("text/latex"), domain) == str
        # test finite case
        domain = IntervalDomain(0.1, 0.1)
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :eq) * " 0.1"
        @test in_domain_string(MIME("text/plain"), domain) == str
        str = InfiniteOpt._math_symbol(MIME("text/latex"), :eq) * " 0.1"
        @test in_domain_string(MIME("text/latex"), domain) == str
    end
    # test in_domain_string (Distribution)
    @testset "in_domain_string (Distribution)" begin
        # test univariate domain
        domain = UniDistributionDomain(Uniform())
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :prop) * " Uniform"
        @test in_domain_string(MIME("text/plain"), domain) == str
        str = InfiniteOpt._math_symbol(MIME("text/latex"), :prop) * " Uniform"
        @test in_domain_string(MIME("text/latex"), domain) == str
        # test mulivariate domain
        domain = MultiDistributionDomain(MvNormal([1], ones(1, 1)))
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :prop) * " MvNormal(dim: (1))"
        str2 = InfiniteOpt._math_symbol(MIME("text/plain"), :prop) * " FullNormal(dim: (1))"
        @test in_domain_string(MIME("text/plain"), domain) in [str, str2]
        str = InfiniteOpt._math_symbol(MIME("text/latex"), :prop) * " MvNormal(dim: (1))"
        str2 = InfiniteOpt._math_symbol(MIME("text/latex"), :prop) * " FullNormal(dim: (1))"
        @test in_domain_string(MIME("text/latex"), domain) in [str, str2]
        # test matrix domain
        domain = MultiDistributionDomain(MatrixBeta(2, 2, 2))
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :prop) * " MatrixBeta(dims: (2, 2))"
        @test in_domain_string(MIME("text/plain"), domain) == str
        str = InfiniteOpt._math_symbol(MIME("text/latex"), :prop) * " MatrixBeta(dims: (2, 2))"
        @test in_domain_string(MIME("text/latex"), domain) == str
    end
    # test in_domain_string (Fallback)
    @testset "in_domain_string (Fallback)" begin
        domain = BadDomain()
        in1 = InfiniteOpt._math_symbol(MIME("text/plain"), :in)
        in2 = InfiniteOpt._math_symbol(MIME("text/latex"), :in)
        @test in_domain_string(MIME("text/plain"), domain) == in1 * " BadDomain()"
        @test in_domain_string(MIME("text/latex"), domain) == in2 * " BadDomain()"
    end
    # test in_domain_string (IntervalDomain with Restrictions)
    @testset "in_domain_string (IntervalDomain with Restrictions)" begin
        # test in restrictions
        rs = DomainRestrictions(par1 => 0)
        domain = IntervalDomain(0, 1)
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :eq) * " 0"
        @test in_domain_string(MIME("text/plain"), par1, domain, rs) == str
        str = InfiniteOpt._math_symbol(MIME("text/latex"), :eq) * " 0"
        @test in_domain_string(MIME("text/latex"), par1, domain, rs) == str
        # test not in restrictions
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]"
        @test in_domain_string(MIME("text/plain"), pars[1], domain, rs) == str
        str = InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1]"
        @test in_domain_string(MIME("text/latex"), pars[1], domain, rs) == str
    end
    # test in_domain_string (InfiniteScalarDomain with Restrictions)
    @testset "in_domain_string (InfiniteScalarDomain with Restrictions)" begin
        # test in restrictions
        rs = DomainRestrictions(par1 => 0)
        domain = UniDistributionDomain(Uniform())
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :eq) * " 0"
        @test in_domain_string(MIME("text/plain"), par1, domain, rs) == str
        str = InfiniteOpt._math_symbol(MIME("text/latex"), :eq) * " 0"
        @test in_domain_string(MIME("text/latex"), par1, domain, rs) == str
        # test in restrictions and not equality
        rs = DomainRestrictions(par1 => [0, 1])
        domain = UniDistributionDomain(Uniform())
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :prop) * " Uniform " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :intersect) * " [0, 1]"
        @test in_domain_string(MIME("text/plain"), par1, domain, rs) == str
        str = InfiniteOpt._math_symbol(MIME("text/latex"), :prop) * " Uniform " *
              InfiniteOpt._math_symbol(MIME("text/latex"), :intersect) * " [0, 1]"
        @test in_domain_string(MIME("text/latex"), par1, domain, rs) == str
        # test not in restrictions
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :prop) * " Uniform"
        @test in_domain_string(MIME("text/plain"), pars[1], domain, rs) == str
        str = InfiniteOpt._math_symbol(MIME("text/latex"), :prop) * " Uniform"
        @test in_domain_string(MIME("text/latex"), pars[1], domain, rs) == str
    end
    # test measure_data_string with 1-D DiscreteMeasureData/FunctionalDiscreteMeasureData
    @testset "measure_data_string (1-D)" begin
        # test with bounds
        data = FunctionalDiscreteMeasureData(par1, ones, 0, All, NoGenerativeSupports(),
                                             default_weight, 0, 1, false)
        str = "par1 " * InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]"
        @test InfiniteOpt.measure_data_string(MIME("text/plain"), data) == str
        str = "par1 " * InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1]"
        @test InfiniteOpt.measure_data_string(MIME("text/latex"), data) == str
        # test without bounds
        data = FunctionalDiscreteMeasureData(par1, ones, 0, All, NoGenerativeSupports(),
                                             default_weight, NaN, NaN, false)
        @test InfiniteOpt.measure_data_string(MIME("text/plain"), data) == "par1"
        @test InfiniteOpt.measure_data_string(MIME("text/latex"), data) == "par1"
    end
    # test measure_data_string with Multi-D DiscreteMeasureData/FunctionalDiscreteMeasureData
    @testset "measure_data_string (Multi-D)" begin
        # test with homogenous bounds
        data = FunctionalDiscreteMeasureData(pars2, ones, 0, All, default_weight, [0, 0], [1, 1], false)
        str = "pars2 " * InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]^2"
        @test InfiniteOpt.measure_data_string(MIME("text/plain"), data) == str
        str = "pars2 " * InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1]^2"
        @test InfiniteOpt.measure_data_string(MIME("text/latex"), data) == str
        # test heterogeneous bounds
        data = FunctionalDiscreteMeasureData(pars2, ones, 0, All, default_weight, [0, 0], [0.5, 1], false)
        str = "pars2[1] " * InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 0.5], " *
              "pars2[2] " * InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]"
        @test InfiniteOpt.measure_data_string(MIME("text/plain"), data) == str
        str = "pars2_{1} " * InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 0.5], " *
              "pars2_{2} " * InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1]"
        @test InfiniteOpt.measure_data_string(MIME("text/latex"), data) == str
        # test no bounds with homogenous names
        data = FunctionalDiscreteMeasureData(pars2, ones, 0, All, default_weight, [NaN, NaN], [NaN, NaN], false)
        @test InfiniteOpt.measure_data_string(MIME("text/plain"), data) == "pars2"
        @test InfiniteOpt.measure_data_string(MIME("text/latex"), data) == "pars2"
        # test no bounds with heterogeneous names
        data = FunctionalDiscreteMeasureData([par1, pars[1]], ones, 0, All, default_weight, [NaN, NaN], [NaN, NaN], false)
        @test InfiniteOpt.measure_data_string(MIME("text/plain"), data) == "[par1, pars[1]]"
        @test InfiniteOpt.measure_data_string(MIME("text/latex"), data) == "[par1, pars[1]]"
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
        @test InfiniteOpt.measure_data_string(MIME("text/plain"), data) == "par1"
        @test InfiniteOpt.measure_data_string(MIME("text/latex"), data) == "par1"
        # test with homogenous names
        data = TestData(pars2, 0, 1)
        @test InfiniteOpt.measure_data_string(MIME("text/plain"), data) == "pars2"
        @test InfiniteOpt.measure_data_string(MIME("text/latex"), data) == "pars2"
        # test with heterogeneous names
        data = TestData([par1, pars[1]], 0, 1)
        @test InfiniteOpt.measure_data_string(MIME("text/plain"), data) == "[par1, pars[1]]"
        @test InfiniteOpt.measure_data_string(MIME("text/latex"), data) == "[par1, pars[1]]"
    end
    # test variabel_string (MeasureRef)
    @testset "variable_string (MeasureRef)" begin
        # test non measure toolbox measures
        data = FunctionalDiscreteMeasureData(pars2, ones, 0, All, default_weight, [0, 0], [1, 1], false)
        meas = dispatch_variable_ref(measure(y, data, name = "test"))
        str = "test{pars2 " * InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]^2}[y]"
        @test InfiniteOpt.variable_string(MIME("text/plain"), meas) == str
        str = "\\text{test}_{pars2 " * InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1]^2}\\left[y\\right]"
        @test InfiniteOpt.variable_string(MIME("text/latex"), meas) == str
        # test measure toolbox special cases for integrals and expectations
        @infinite_parameter(m, t in [0, 1])
        meas = dispatch_variable_ref(expect(y, t))
        str = "\\mathbb{E}_{t \\in [0, 1]}\\left[y\\right]"
        @test InfiniteOpt.variable_string(MIME("text/latex"), meas) == str
        str = InfiniteOpt._math_symbol(MIME("text/plain"), :expect) * "{t " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]}[y]"
        @test InfiniteOpt.variable_string(MIME("text/plain"), meas) == str
        meas = dispatch_variable_ref(integral(y, t))
        str = "\\int_{t " * InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1]}ydt"
        @test InfiniteOpt.variable_string(MIME("text/latex"), meas) == str
        int = InfiniteOpt._math_symbol(MIME("text/plain"), :integral)
        str = int * "{t " * InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]}[y]"
        @test InfiniteOpt.variable_string(MIME("text/plain"), meas) == str
        meas = dispatch_variable_ref(integral(y, par1))
        str = "\\int_{par1 " * InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1]}yd(par1)"
        @test InfiniteOpt.variable_string(MIME("text/latex"), meas) == str
    end
    # test _get_base_name
    @testset "_get_base_name (REPL)" begin
        @test InfiniteOpt._get_base_name(MIME("text/plain"), pars[1]) == "pars[1]"
        bad_ref = FiniteVariableRef(m, FiniteVariableIndex(-1))
        @test InfiniteOpt._get_base_name(MIME("text/plain"), bad_ref) == "noname"
    end
    # test _get_base_name
    @testset "_get_base_name (IJulia)" begin
        @test InfiniteOpt._get_base_name(MIME("text/latex"), pars[1]) == "pars_{1}"
        bad_ref = FiniteVariableRef(m, FiniteVariableIndex(-1))
        @test InfiniteOpt._get_base_name(MIME("text/latex"), bad_ref) == "noname"
    end
    # test variable_string (InfiniteVariableRef)
    @testset "variable_string (InfiniteVariableRef)" begin
        # test bad ref
        bad_ref = InfiniteVariableRef(m, InfiniteVariableIndex(-1))
        @test InfiniteOpt.variable_string(MIME("text/plain"), bad_ref) == "noname"
        @test InfiniteOpt.variable_string(MIME("text/latex"), bad_ref) == "noname"
        # test complex normal one
        dvref = dispatch_variable_ref(inf)
        @test InfiniteOpt.variable_string(MIME("text/plain"), dvref) == "inf(pars, par1, pars3)"
        @test InfiniteOpt.variable_string(MIME("text/latex"), dvref) == "inf(pars, par1, pars3)"
        # test mixed types
        @variable(m, inf2, Infinite([par1, pars3[2]]))
        dvref = dispatch_variable_ref(inf2)
        @test InfiniteOpt.variable_string(MIME("text/plain"), dvref) == "inf2([par1, pars3[2]])"
        @test InfiniteOpt.variable_string(MIME("text/latex"), dvref) == "inf2([par1, pars3[2]])"
    end
    # test variable_string (ParameterFunctionRef)
    @testset "variable_string (ParameterFunctionRef)" begin
        # test bad ref
        bad_ref = ParameterFunctionRef(m, ParameterFunctionIndex(-1))
        @test InfiniteOpt.variable_string(MIME("text/plain"), bad_ref) == "noname"
        @test InfiniteOpt.variable_string(MIME("text/latex"), bad_ref) == "noname"
        # test normal one
        dvref = dispatch_variable_ref(parameter_function(sin, par1))
        @test InfiniteOpt.variable_string(MIME("text/plain"), dvref) == "sin(par1)"
        @test InfiniteOpt.variable_string(MIME("text/latex"), dvref) == "sin(par1)"
    end
    # test variable_string (DerivativeRef)
    @testset "variable_string (DerivativeRef)" begin
        # test bad ref
        bad_ref = DerivativeRef(m, DerivativeIndex(-1))
        @test InfiniteOpt.variable_string(MIME("text/plain"), bad_ref) == "noname"
        @test InfiniteOpt.variable_string(MIME("text/latex"), bad_ref) == "noname"
        # test normal one
        dref = dispatch_variable_ref(d1)
        d_re = InfiniteOpt._math_symbol(MIME("text/plain"), :partial)
        @test InfiniteOpt.variable_string(MIME("text/plain"), dref) == "$d_re/$(d_re)par1[x(par1)]"
        @test InfiniteOpt.variable_string(MIME("text/latex"), dref) == "\\frac{\\partial}{\\partial par1}\\left[x(par1)\\right]"
        # test nested one
        dref = dispatch_variable_ref(d2)
        expected = "$d_re/$(d_re)par1[$d_re/$(d_re)pars[1][inf(pars, par1, pars3)]]"
        @test InfiniteOpt.variable_string(MIME("text/plain"), dref) == expected
        expected = "\\frac{\\partial}{\\partial par1}\\left[\\frac{\\partial}{\\partial pars_{1}}\\left[inf(pars, par1, pars3)\\right]\\right]"
        @test InfiniteOpt.variable_string(MIME("text/latex"), dref) == expected
        # test has explicit name
        dref = dispatch_variable_ref(d3)
        @test InfiniteOpt.variable_string(MIME("text/plain"), dref) == "d3(pars)"
        @test InfiniteOpt.variable_string(MIME("text/latex"), dref) == "d3(pars)"
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
        @test InfiniteOpt.variable_string(MIME("text/plain"), bad_ref) == "noname"
        @test InfiniteOpt.variable_string(MIME("text/latex"), bad_ref) == "noname"
        # test short one
        dvref = dispatch_variable_ref(@variable(m, variable_type = Point(x, 0)))
        @test InfiniteOpt.variable_string(MIME("text/plain"), dvref) == "x(0)"
        @test InfiniteOpt.variable_string(MIME("text/latex"), dvref) == "x(0)"
        # test complex one
        dvref = dispatch_variable_ref(@variable(m, variable_type = Point(inf, [0, 0], 0, [0, 0])))
        @test InfiniteOpt.variable_string(MIME("text/plain"), dvref) == "inf([0, 0], 0, [0, 0])"
        @test InfiniteOpt.variable_string(MIME("text/latex"), dvref) == "inf([0, 0], 0, [0, 0])"
        # test with alias name
        dvref = dispatch_variable_ref(@variable(m, z0, Point(z, [0, 0])))
        @test InfiniteOpt.variable_string(MIME("text/plain"), dvref) == "z0"
        @test InfiniteOpt.variable_string(MIME("text/latex"), dvref) == "z0"
        # test with derivative
        dvref = dispatch_variable_ref(@variable(m, variable_type = Point(d1, 0)))
        d_re = InfiniteOpt._math_symbol(MIME("text/plain"), :partial)
        @test InfiniteOpt.variable_string(MIME("text/plain"), dvref) == "$(d_re)/$(d_re)par1[x(par1)](0)"
        # test named derivative
        dvref = dispatch_variable_ref(@variable(m, variable_type = Point(d3, [0, 0])))
        @test InfiniteOpt.variable_string(MIME("text/plain"), dvref) == "d3([0, 0])"
    end
    # test variable_string (SemiInfiniteVariableRef)
    @testset "variable_string (SemiInfiniteVariableRef)" begin
        # test bad ref
        bad_ref = SemiInfiniteVariableRef(m, SemiInfiniteVariableIndex(-1))
        @test InfiniteOpt.variable_string(MIME("text/plain"), bad_ref) == "noname"
        @test InfiniteOpt.variable_string(MIME("text/latex"), bad_ref) == "noname"
        # test short one
        rv = @variable(m, variable_type = SemiInfinite(inf, [0, 0], par1, [0, pars3[2]]))
        dvref = dispatch_variable_ref(rv)
        @test InfiniteOpt.variable_string(MIME("text/plain"), dvref) == "inf([0, 0], par1, [0, pars3[2]])"
        @test InfiniteOpt.variable_string(MIME("text/latex"), dvref) == "inf([0, 0], par1, [0, pars3[2]])"
        # test named derivative
        rv = @variable(m, variable_type = SemiInfinite(d3, [0, pars[2]]))
        dvref = dispatch_variable_ref(rv)
        @test InfiniteOpt.variable_string(MIME("text/plain"), dvref) == "d3([0, pars[2]])"
        # unnamed derivative
        eval_supps = Dict{Int, Float64}(1 => 0)
        var = build_variable(error, d1, eval_supps, check = false)
        rv = @variable(m, variable_type = SemiInfinite(d1, 0))
        dvref = dispatch_variable_ref(rv)
        d_re = InfiniteOpt._math_symbol(MIME("text/plain"), :partial)
        @test InfiniteOpt.variable_string(MIME("text/plain"), dvref) == "$(d_re)/$(d_re)par1[x(par1)](0)"
    end
    # test variable_string (Fallback)
    @testset "variable_string (Fallback)" begin
        @test InfiniteOpt.variable_string(MIME("text/plain"), y) == "y"
        @test InfiniteOpt.variable_string(MIME("text/latex"), y) == "y"
    end
    # test the function_string extensions
    @testset "JuMP.function_string" begin
        @test JuMP.function_string(MIME("text/plain"), dispatch_variable_ref(y)) == "y"
        @test JuMP.function_string(MIME("text/latex"), dispatch_variable_ref(inf)) == "inf(pars, par1, pars3)"
        @test JuMP.function_string(MIME("text/plain"), y) == "y"
        @test JuMP.function_string(MIME("text/latex"), inf) == "inf(pars, par1, pars3)"
        d_re = InfiniteOpt._math_symbol(MIME("text/plain"), :partial)
        @test JuMP.function_string(MIME("text/plain"), d1) == "$d_re/$(d_re)par1[x(par1)]"
    end
    # test restrict_string
    @testset "restrict_string" begin
        # test with single restriction
        rs = DomainRestrictions(par1 => [0.5, 0.7])
        str = "par1 " * InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0.5, 0.7]"
        @test InfiniteOpt.restrict_string(MIME("text/plain"), rs) == str
        str = "par1 " *  InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0.5, 0.7]"
        @test InfiniteOpt.restrict_string(MIME("text/latex"), rs) == str
    end
    # test constraint_string (Finite constraint)
    @testset "JuMP.constraint_string (Finite)" begin
        # test named
        str = "c2 : y² " * in_set_string(MIME("text/plain"), MOI.EqualTo(3.0))
        @test constraint_string(MIME("text/plain"), c2) == str
        str =  "c2 : \$ y^2 " * in_set_string(MIME("text/latex"), MOI.EqualTo(3.0)) * " \$"
        @test constraint_string(MIME("text/latex"), c2) == str
        # test unnamed
        str = "y² " * in_set_string(MIME("text/plain"), MOI.EqualTo(3.0))
        @test constraint_string(MIME("text/plain"), ac2) == str
        str =  "\$ y^2 " * in_set_string(MIME("text/latex"), MOI.EqualTo(3.0)) * " \$"
        @test constraint_string(MIME("text/latex"), ac2) == str
        # test named in math mode
        str = "c2 : y² " * in_set_string(MIME("text/plain"), MOI.EqualTo(3.0))
        @test constraint_string(MIME("text/plain"), c2, in_math_mode = true) == str
        str =  "y^2 " * in_set_string(MIME("text/latex"), MOI.EqualTo(3.0))
        @test constraint_string(MIME("text/latex"), c2, in_math_mode = true) == str
    end
    # test _param_domain_string (IndependentParameter)
    @testset "_param_domain_string (IndependentParameter)" begin
        rs = DomainRestrictions(par1 => 0)
        idx = index(par1)
        str = "par1 " * InfiniteOpt._math_symbol(MIME("text/plain"), :eq) * " 0"
        @test InfiniteOpt._param_domain_string(MIME("text/plain"), m, idx, rs) == str
        str = "par1 " * InfiniteOpt._math_symbol(MIME("text/latex"), :eq) * " 0"
        @test InfiniteOpt._param_domain_string(MIME("text/latex"), m, idx, rs) == str
    end
    # test _param_domain_string (DependentParameters)
    @testset "_param_domain_string (DependentParameters)" begin
        # Collection set
        rs = DomainRestrictions(pars2[1] => [0, 1])
        idx = index(pars2[1]).object_index
        str = "pars2[1] " * InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1], " *
              "pars2[2] " * InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 2]"
        @test InfiniteOpt._param_domain_string(MIME("text/plain"), m, idx, rs) == str
        str = "pars2_{1} " * InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1], " *
              "pars2_{2} " * InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 2]"
        @test InfiniteOpt._param_domain_string(MIME("text/latex"), m, idx, rs) == str
        # other set with equalities
        rs = DomainRestrictions(pars[1] => 0, pars[2] => 1)
        idx = index(pars[1]).object_index
        str = InfiniteOpt.restrict_string(MIME("text/plain"), rs)
        str2 = string(split(str, ", ")[2], ", ", split(str, ", ")[1])
        @test InfiniteOpt._param_domain_string(MIME("text/plain"), m, idx, rs) in [str, str2]
        str = InfiniteOpt.restrict_string(MIME("text/latex"), rs)
        str2 = string(split(str, ", ")[2], ", ", split(str, ", ")[1])
        @test InfiniteOpt._param_domain_string(MIME("text/latex"), m, idx, rs) in [str, str2]
        # other set without equalities and including in the restrictions
        rs = DomainRestrictions(pars[1] => [0, 1])
        str = "pars " * InfiniteOpt._math_symbol(MIME("text/plain"), :prop) *
              " MvNormal(dim: (2)) " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :intersect) *
              " (pars[1] " * InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1])"
        str2 = "pars " * InfiniteOpt._math_symbol(MIME("text/plain"), :prop) *
              " FullNormal(dim: (2)) " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :intersect) *
              " (pars[1] " * InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1])"
        @test InfiniteOpt._param_domain_string(MIME("text/plain"), m, idx, rs) in [str, str2]
        str = "pars " * InfiniteOpt._math_symbol(MIME("text/latex"), :prop) *
              " MvNormal(dim: (2)) " *
              InfiniteOpt._math_symbol(MIME("text/latex"), :intersect) *
              " (pars_{1} " * InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1])"
        str2 = "pars " * InfiniteOpt._math_symbol(MIME("text/latex"), :prop) *
              " FullNormal(dim: (2)) " *
              InfiniteOpt._math_symbol(MIME("text/latex"), :intersect) *
              " (pars_{1} " * InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1])"
        @test InfiniteOpt._param_domain_string(MIME("text/latex"), m, idx, rs) in [str, str2]
        # other set without equalities and not included in restrictions
        rs = DomainRestrictions(par1 => [0, 1])
        str = "pars " * InfiniteOpt._math_symbol(MIME("text/plain"), :prop) *
              " MvNormal(dim: (2))"
        str2 = "pars " * InfiniteOpt._math_symbol(MIME("text/plain"), :prop) *
              " FullNormal(dim: (2))"
        @test InfiniteOpt._param_domain_string(MIME("text/plain"), m, idx, rs) in [str, str2]
        str = "pars " * InfiniteOpt._math_symbol(MIME("text/latex"), :prop) *
              " MvNormal(dim: (2))"
        str2 = "pars " * InfiniteOpt._math_symbol(MIME("text/latex"), :prop) *
              " FullNormal(dim: (2))"
        @test InfiniteOpt._param_domain_string(MIME("text/latex"), m, idx, rs) in [str, str2]
    end
    # test constraint_string (infinite constraint)
    @testset "JuMP.constraint_string (Infinite)" begin
        # test c1 with name
        str = "c1 : x(par1) + y " * in_set_string(MIME("text/plain"), MOI.LessThan(2.0)) * ", " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) * " par1 " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]"
        @test constraint_string(MIME("text/plain"), c1) == str
        str = "c1 : \$ x(par1) + y " * in_set_string(MIME("text/latex"), MOI.LessThan(2.0)) * ", " *
              InfiniteOpt._math_symbol(MIME("text/latex"), :for_all) * " par1 " *
              InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1] \$"
        @test constraint_string(MIME("text/latex"), c1) == str
        # test c1 without name
        str = "x(par1) + y " * in_set_string(MIME("text/plain"), MOI.LessThan(2.0)) * ", " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) * " par1 " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]"
        @test constraint_string(MIME("text/plain"), ac1) == str
        str = "\$ x(par1) + y " * in_set_string(MIME("text/latex"), MOI.LessThan(2.0)) * ", " *
              InfiniteOpt._math_symbol(MIME("text/latex"), :for_all) * " par1 " *
              InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1] \$"
        @test constraint_string(MIME("text/latex"), ac1) == str
        # test c1 with name and in_math_mode
        str = "c1 : x(par1) + y " * in_set_string(MIME("text/plain"), MOI.LessThan(2.0)) * ", " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) * " par1 " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]"
        @test constraint_string(MIME("text/plain"), c1, in_math_mode = true) == str
        str = "x(par1) + y " * in_set_string(MIME("text/latex"), MOI.LessThan(2.0)) * ", " *
              InfiniteOpt._math_symbol(MIME("text/latex"), :for_all) * " par1 " *
              InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1]"
        @test constraint_string(MIME("text/latex"), c1, in_math_mode = true) == str
        # test c3 with name
        str = "c3 : x(par1) " * in_set_string(MIME("text/plain"), MOI.EqualTo(5.0)) * ", " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) * " par1 " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 0.5]"
        @test constraint_string(MIME("text/plain"), c3) == str
        str =  "c3 : \$ x(par1) " * in_set_string(MIME("text/latex"), MOI.EqualTo(5.0)) * ", " *
               InfiniteOpt._math_symbol(MIME("text/latex"), :for_all) * " par1 " *
               InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 0.5] \$"
        @test constraint_string(MIME("text/latex"), c3) == str
        # test c3 without name
        str = "x(par1) " * in_set_string(MIME("text/plain"), MOI.EqualTo(5.0)) * ", " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) * " par1 " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 0.5]"
        @test constraint_string(MIME("text/plain"), ac3) == str
        str =  "\$ x(par1) " * in_set_string(MIME("text/latex"), MOI.EqualTo(5.0)) * ", " *
               InfiniteOpt._math_symbol(MIME("text/latex"), :for_all) * " par1 " *
               InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 0.5] \$"
        @test constraint_string(MIME("text/latex"), ac3) == str
    end
    # test constraints_string
    @testset "JuMP.constraints_string" begin
        # test MIME("text/plain")
        strings = Vector{String}(undef, 6)
        strings[1] = "c1 : x(par1) + y " * in_set_string(MIME("text/plain"), MOI.LessThan(2.0)) * 
                     ", " * InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) * " par1 " *
                     InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]"
        strings[2] = "x(par1) + y " * in_set_string(MIME("text/plain"), MOI.LessThan(2.0)) *
                     ", " * InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) * " par1 " *
                     InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]"
        strings[3] = "c2 : y² " * in_set_string(MIME("text/plain"), MOI.EqualTo(3.0))
        strings[4] = "y² " * in_set_string(MIME("text/plain"), MOI.EqualTo(3.0))
        strings[5] = "c3 : x(par1) " * in_set_string(MIME("text/plain"), MOI.EqualTo(5.0)) * ", " *
                     InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) * " par1 " *
                     InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 0.5]"
        strings[6] = "x(par1) " * in_set_string(MIME("text/plain"), MOI.EqualTo(5.0)) * ", " *
                     InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) * " par1 " *
                     InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 0.5]"
        @test constraints_string(MIME("text/plain"), m) == strings
        # test MIME("text/latex")
        strings = Vector{String}(undef, 6)
        strings[1] = "x(par1) + y " * in_set_string(MIME("text/latex"), MOI.LessThan(2.0)) *
                     ", " * InfiniteOpt._math_symbol(MIME("text/latex"), :for_all) * " par1 " *
                     InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1]"
        strings[2] = "x(par1) + y " * in_set_string(MIME("text/latex"), MOI.LessThan(2.0)) *
                     ", " * InfiniteOpt._math_symbol(MIME("text/latex"), :for_all) * " par1 " *
                     InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1]"
        strings[3] = "y^2 " * in_set_string(MIME("text/latex"), MOI.EqualTo(3.0))
        strings[4] = "y^2 " * in_set_string(MIME("text/latex"), MOI.EqualTo(3.0))
        strings[5] = "x(par1) " * in_set_string(MIME("text/latex"), MOI.EqualTo(5.0)) * ", " *
                     InfiniteOpt._math_symbol(MIME("text/latex"), :for_all) * " par1 " *
                     InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 0.5]"
        strings[6] = "x(par1) " * in_set_string(MIME("text/latex"), MOI.EqualTo(5.0)) * ", " *
                     InfiniteOpt._math_symbol(MIME("text/latex"), :for_all) * " par1 " *
                     InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 0.5]"
        @test constraints_string(MIME("text/latex"), m) == strings
    end
    # test objective_function_string
    @testset "JuMP.objective_function_string" begin
        @test objective_function_string(MIME("text/plain"), m) == "y + 2"
        @test objective_function_string(MIME("text/latex"), m) == "y + 2"
    end
end

# test infinite model show
@testset "Show InfiniteModel" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, par1 in [0, 1])
    @infinite_parameter(m, pars[1:2] ~ MvNormal([1, 1], [1 0; 0 1]))
    @variable(m, x, Infinite(par1))
    @variable(m, z, Infinite(pars))
    @variable(m, y)
    @objective(m, Min, 2 + y)
    @constraint(m, c1, x + y -2 <= 0)
    rs = DomainRestrictions(par1 => [0.1, 1])
    @constraint(m, c3, x <= 5, DomainRestrictions(par1 => [0, 0.5]))
    mockoptimizer = () -> MOIU.MockOptimizer(MOIU.UniversalFallback(MOIU.Model{Float64}()),
                                             eval_objective_value=false)
    # test Base.show (IntervalDomain in REPL)
    @testset "Base.show (REPL IntervalDomain)" begin
        show_test(MIME("text/plain"), IntervalDomain(0, 1), "[0, 1]")
    end
    # test Base.show (IntervalDomain in IJulia)
    @testset "Base.show (IJulia IntervalDomain)" begin
        show_test(MIME("text/latex"), IntervalDomain(0, 1), "[0, 1]")
    end
    # test Base.show (DistributionDomain in REPL)
    @testset "Base.show (REPL DistributionDomain)" begin
        show_test(MIME("text/plain"), UniDistributionDomain(Uniform()), string(Uniform()))
    end
    # test Base.show (DistributionDomain in IJulia)
    @testset "Base.show (IJulia DistributionDomain)" begin
        show_test(MIME("text/latex"), UniDistributionDomain(Uniform()), string(Uniform()))
    end
    # test Base.show (CollectionDomain in REPL)
    @testset "Base.show (REPL CollectionDomain)" begin
        show_test(MIME("text/plain"), CollectionDomain([IntervalDomain(0, 0)]),
                  "CollectionDomain with 1 domain:\n [0, 0]")
    end
    # test Base.show (CollectionDomain in IJulia)
    @testset "Base.show (IJulia CollectionDomain)" begin
        show_test(MIME("text/latex"), CollectionDomain([IntervalDomain(0, 0), IntervalDomain(0, 2)]),
                  "CollectionDomain with 2 domains:\n [0, 0]\n [0, 2]")
    end
    # test Base.show (DomainRestrictions in REPL)
    @testset "Base.show (REPL DomainRestrictions)" begin
        str = "Subdomain restrictions (1): par1 " * InfiniteOpt._math_symbol(MIME("text/plain"), :in) *
              " [0.1, 1]"
        show_test(MIME("text/plain"), rs, str)
    end
    # test Base.show (DomainRestrictions in IJulia)
    @testset "Base.show (IJulia DomainRestrictions)" begin
        str = "Subdomain restrictions (1): par1 " * InfiniteOpt._math_symbol(MIME("text/latex"), :in) *
              " [0.1, 1]"
        show_test(MIME("text/latex"), rs, str)
    end
    # test Base.show (GeneralVariableRef in IJulia)
    @testset "Base.show (IJulia GeneralVariableRef)" begin
        show_test(MIME("text/latex"), y, "\$\$ y \$\$")
    end
    # test Base.show (GeneralVariableRef in REPL)
    @testset "Base.show (REPL GeneralVariableRef)" begin
        show_test(MIME("text/plain"), y, "y")
    end
    # test Base.show (constraint in REPL)
    @testset "Base.show (REPL Constraint)" begin
        # test normal
        str = "c1 : x(par1) + y " * in_set_string(MIME("text/plain"), MOI.LessThan(2.0)) *
              ", " * InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) * " par1 " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 1]"
        show_test(MIME("text/plain"), c1, str)
        # test restricted
        str = "c3 : x(par1) " * in_set_string(MIME("text/plain"), MOI.LessThan(5.0)) * ", " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :for_all) * " par1 " *
              InfiniteOpt._math_symbol(MIME("text/plain"), :in) * " [0, 0.5]"
        show_test(MIME("text/plain"), c3, str)
    end
    # test Base.show (constraint in IJulia)
    @testset "Base.show (IJulia Constraint)" begin
        # test normal
        str = "c1 : \$ x(par1) + y " * in_set_string(MIME("text/latex"), MOI.LessThan(2.0)) *
              ", " * InfiniteOpt._math_symbol(MIME("text/latex"), :for_all) * " par1 " *
              InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 1] \$"
        show_test(MIME("text/latex"), c1, str)
        # test restricted
        str =  "c3 : \$ x(par1) " * in_set_string(MIME("text/latex"), MOI.LessThan(5.0)) * ", " *
               InfiniteOpt._math_symbol(MIME("text/latex"), :for_all) * " par1 " *
               InfiniteOpt._math_symbol(MIME("text/latex"), :in) * " [0, 0.5] \$"
        show_test(MIME("text/latex"), c3, str)
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
        str2 = "Objective function type: GenericAffExpr{Float64, GeneralVariableRef}\n"
        io_test(show_objective_function_summary, [str, str2], m)
    end
    # test show_constraints_summary
    @testset "JuMP.show_constraints_summary" begin
        # test the main function
        str = "`GenericAffExpr{Float64,GeneralVariableRef}`-in-`MathOptInter" *
              "face.LessThan{Float64}`: 2 constraints\n"
        str2 = "`GenericAffExpr{Float64, GeneralVariableRef}`-in-`MathOptInter" *
              "face.LessThan{Float64}`: 2 constraints\n"
        io_test(show_constraints_summary, [str, str2], m)
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
        str2 = "An InfiniteOpt Model\nMinimization problem with:\nFinite " *
              "Parameters: 0\nInfinite Parameters: 3\nVariables: 3" *
              "\nDerivatives: 0\nMeasures: 0" *
              "\nObjective function type: GenericAffExpr{Float64, General" *
              "VariableRef}\n`GenericAffExpr{Float64, GeneralVariableRef}`-in-" *
              "`MathOptInterface.LessThan{Float64}`: 2 constraints" *
              "\nNames registered in the model: c1, c3, par1, " *
              "pars, x, y, z\nOptimizer model backend information: \nModel " *
              "mode: AUTOMATIC\nCachingOptimizer state: EMPTY_OPTIMIZER\n" *
              "Solver name: Mock"
        show_test(MIME("text/plain"), m, [str, str2], repl=:show)
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
        str2 = "An InfiniteOpt Model\nMaximization problem with:\nFinite " *
              "Parameters: 0\nInfinite Parameters: 3\nVariables: 3" *
              "\nDerivatives: 0\nMeasures: 0" *
              "\nObjective function type: GenericAffExpr{Float64, General" *
              "VariableRef}\n`GenericAffExpr{Float64, GeneralVariableRef}`-in-" *
              "`MathOptInterface.LessThan{Float64}`: 2 constraints" *
              "\nNames registered in the model: c1, c3, par1, " *
              "pars, x, y, z\nOptimizer model backend information: \nModel " *
              "mode: AUTOMATIC\nCachingOptimizer state: EMPTY_OPTIMIZER\n" *
              "Solver name: Mock"
        show_test(MIME("text/plain"), m, [str, str2], repl=:show)
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
        str2 = "An InfiniteOpt Model\nFeasibility problem with:\nFinite " *
              "Parameters: 0\nInfinite Parameters: 3\nVariables: 3" *
              "\nDerivatives: 0\nMeasures: 0" *
              "\n`GenericAffExpr{Float64, GeneralVariableRef}`-in-`MathOpt" *
              "Interface.LessThan{Float64}`: 2 constraints" *
              "\nNames registered in the model: c1, c3, par1, " *
              "pars, x, y, z\nOptimizer model backend information: \nModel " *
              "mode: AUTOMATIC\nCachingOptimizer state: EMPTY_OPTIMIZER\n" *
              "Solver name: Mock"
        show_test(MIME("text/plain"), m, [str, str2], repl=:show)
    end
end
