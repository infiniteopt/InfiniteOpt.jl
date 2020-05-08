using JuMP: REPLMode, IJuliaMode

# Test string creation
@testset "String Creators" begin
    # initialize model and attributes
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par1 <= 1)
    @infinite_parameter(m, pars[1:2] in MvNormal([1, 1], 1))
    @infinite_parameter(m, pars2[1:2] in [0, 2])
    @infinite_variable(m, x(par1))
    @infinite_variable(m, z(pars))
    @hold_variable(m, y)
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
        else
            @test InfiniteOpt._infopt_math_symbol(REPLMode, :intersect) == "∩"
        end
        @test InfiniteOpt._infopt_math_symbol(REPLMode, :times) == "*"
        @test InfiniteOpt._infopt_math_symbol(REPLMode, :prop) == "~"
    end
    # test _infopt_math_symbol (IJulia)
    @testset "_infopt_math_symbol (IJulia)" begin
        @test InfiniteOpt._infopt_math_symbol(IJuliaMode, :intersect) == "\\cap"
        @test InfiniteOpt._infopt_math_symbol(IJuliaMode, :eq) == "="
        @test InfiniteOpt._infopt_math_symbol(IJuliaMode, :prop) == "\\sim"
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
        @test set_string(REPLMode, set) == "[0, 1]"
        @test set_string(IJuliaMode, set) == "[0, 1]"
        # test rounding case
        set = IntervalSet(-0, 1)
        @test set_string(REPLMode, set) == "[0, 1]"
        @test set_string(IJuliaMode, set) == "[0, 1]"
        # test decimal case
        set = IntervalSet(0.1, 1.3)
        @test set_string(REPLMode, set) == "[0.1, 1.3]"
        @test set_string(IJuliaMode, set) == "[0.1, 1.3]"
    end
    # test set_string (DistributionSet)
    @testset "set_string (DistributionSet)" begin
        # test univariate set
        set = UniDistributionSet(Uniform())
        @test set_string(REPLMode, set) == "Uniform{Float64}(a=0.0, b=1.0)"
        @test set_string(IJuliaMode, set) == "Uniform{Float64}(a=0.0, b=1.0)"
        # test mulivariate set
        set = MultiDistributionSet(MvNormal([1], 1))
        str = "IsoNormal(\ndim: 1\nμ: [1.0]\nΣ: [1.0]\n)\n"
        @test set_string(REPLMode, set) == str
        @test set_string(IJuliaMode, set) == str
    end
    # test set_string (CollectionSet)
    @testset "set_string (CollectionSet)" begin
        set = CollectionSet([IntervalSet(0, 1), IntervalSet(0, 0.1)])
        str = "CollectionSet with 2 sets:\n [0, 1]\n [0, 0.1]"
        @test set_string(REPLMode, set) == str
        @test set_string(IJuliaMode, set) == str
    end
    # test set_string (Fallback)
    @testset "set_string (Fallback)" begin
        set = BadSet()
        @test set_string(REPLMode, set) == "BadSet()"
        @test set_string(IJuliaMode, set) == "BadSet()"
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
        str = bound_string(REPLMode, bounds)
        @test InfiniteOpt._param_domain_string(REPLMode, m, idx, bounds) == str
        str = bound_string(IJuliaMode, bounds)
        @test InfiniteOpt._param_domain_string(IJuliaMode, m, idx, bounds) == str
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
              "Parameters: 0\nInfinite Parameters: 3\nVariables: " *
              "3\nObjective function type: GenericAffExpr{Float64,General" *
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
              "Parameters: 0\nInfinite Parameters: 3\nVariables: " *
              "3\nObjective function type: GenericAffExpr{Float64,General" *
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
              "\n`GenericAffExpr{Float64,GeneralVariableRef}`-in-`MathOpt" *
              "Interface.LessThan{Float64}`: 2 constraints" *
              "\nNames registered in the model: c1, c3, par1, " *
              "pars, x, y, z\nOptimizer model backend information: \nModel " *
              "mode: AUTOMATIC\nCachingOptimizer state: EMPTY_OPTIMIZER\n" *
              "Solver name: Mock"
        show_test(REPLMode, m, str, repl=:show)
    end
end
