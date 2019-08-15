using JuMP: REPLMode, IJuliaMode

# Test string creation
@testset "String Creators" begin
    # initialize model and attributes
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par1 <= 1)
    @infinite_parameter(m, pars[1:2] in MvNormal([1, 1], 1))
    @infinite_variable(m, x(par1))
    @infinite_variable(m, z(pars))
    @global_variable(m, y)
    @objective(m, Min, 2 + y)
    @constraint(m, c1, x + y -2 <= 0)
    @constraint(m, c2, y^2 - 3 == 0)
    @constraint(m, c3, x == 5,
                parameter_bounds = Dict(par1 => IntervalSet(0, 0.5)))
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
    end
    # test in_set_string (Distribution)
    @testset "JuMP.in_set_string (Distribution)" begin
        # test univariate set
        set = DistributionSet(Uniform())
        str = JuMP._math_symbol(REPLMode, :in) * " Uniform(a=0.0, b=1.0)"
        @test in_set_string(REPLMode, set) == str
        str = JuMP._math_symbol(IJuliaMode, :in) * " Uniform(a=0.0, b=1.0)"
        @test in_set_string(IJuliaMode, set) == str
        # test mulivariate set
        set = DistributionSet(MvNormal([1], 1))
        str = JuMP._math_symbol(REPLMode, :in) * " IsoNormal(dim: 1, μ: " *
                 "[1.0], Σ: [1.0])"
        @test in_set_string(REPLMode, set) == str
        str = JuMP._math_symbol(IJuliaMode, :in) * " IsoNormal(dim: 1, μ: " *
                 "[1.0], Σ: [1.0])"
        @test in_set_string(IJuliaMode, set) == str
    end
    # test bound_string
    @testset "bound_string" begin
        # test with single bound
        bounds = Dict(par1 => IntervalSet(0.5, 0.7))
        str = JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
                 JuMP._math_symbol(REPLMode, :in) * " [0.5, 0.7]"
        @test InfiniteOpt.bound_string(REPLMode, bounds) == str
        str = JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
                 JuMP._math_symbol(IJuliaMode, :in) * " [0.5, 0.7]"
        @test InfiniteOpt.bound_string(IJuliaMode, bounds) == str
    end
    # test constraint_string (bounded constraints)
    @testset "JuMP.constraint_string (Bounded)" begin
        str = "x(par1) " * JuMP._math_symbol(REPLMode, :eq) * " 5.0, " *
              JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
              JuMP._math_symbol(REPLMode, :in) * " [0, 0.5]"
        @test constraint_string(REPLMode, constraint_object(c3)) == str
        str =  "x(par1) " * JuMP._math_symbol(IJuliaMode, :eq) * " 5.0, " *
               JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
               JuMP._math_symbol(IJuliaMode, :in) * " [0, 0.5]"
        @test constraint_string(IJuliaMode, constraint_object(c3)) == str
    end
    # test constraint_string of GeneralConstraintRefs
    @testset "JuMP.constraint_string (Reference)" begin
        # test normal constraint
        str = "c1 : x(par1) + y " * JuMP._math_symbol(REPLMode, :leq) *
              " 2.0"
        @test constraint_string(REPLMode, c1) == str
        str = "c1 : \$ x(par1) + y " * JuMP._math_symbol(IJuliaMode, :leq) *
              " 2.0 \$"
        @test constraint_string(IJuliaMode, c1) == str
        # test bounded constraint
        str = "c3 : x(par1) " * JuMP._math_symbol(REPLMode, :eq) * " 5.0, " *
              JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
              JuMP._math_symbol(REPLMode, :in) * " [0, 0.5]"
        @test constraint_string(REPLMode, c3) == str
        str =  "c3 : \$ x(par1) " * JuMP._math_symbol(IJuliaMode, :eq) * " 5.0, " *
               JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
               JuMP._math_symbol(IJuliaMode, :in) * " [0, 0.5] \$"
        @test constraint_string(IJuliaMode, c3) == str
    end
    # test constraints_string
    @testset "JuMP.constraints_string" begin
        # test REPLMode
        strings = Vector{String}(undef, 5)
        strings[1] = "x(par1) + y " * JuMP._math_symbol(REPLMode, :leq) * " 2.0"
        strings[2] = "y² " * JuMP._math_symbol(REPLMode, :eq) * " 3.0"
        strings[3] = "x(par1) " * JuMP._math_symbol(REPLMode, :eq) * " 5.0, " *
                     JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
                     JuMP._math_symbol(REPLMode, :in) * " [0, 0.5]"
        strings[4] = "par1 " * JuMP._math_symbol(REPLMode, :in) * " [0, 1]"
        strings[5] = "pars " * JuMP._math_symbol(REPLMode, :in) *
                     " IsoNormal(dim: 2, μ: [1.0, 1.0], Σ: [1.0 0.0; 0.0 1.0])"
        @test constraints_string(REPLMode, m) == strings
        # test IJuliaMode
        strings = Vector{String}(undef, 5)
        strings[1] = "x(par1) + y " * JuMP._math_symbol(IJuliaMode, :leq) * " 2.0"
        strings[2] = "y^2 " * JuMP._math_symbol(IJuliaMode, :eq) * " 3.0"
        strings[3] = "x(par1) " * JuMP._math_symbol(IJuliaMode, :eq) * " 5.0, " *
                     JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
                     JuMP._math_symbol(IJuliaMode, :in) * " [0, 0.5]"
        strings[4] = "par1 " * JuMP._math_symbol(IJuliaMode, :in) * " [0, 1]"
        strings[5] = "pars " * JuMP._math_symbol(IJuliaMode, :in) *
                     " IsoNormal(dim: 2, μ: [1.0, 1.0], Σ: [1.0 0.0; 0.0 1.0])"
        @test constraints_string(IJuliaMode, m) == strings
    end
    # test objective_function_string
    @testset "JuMP.objective_function_string" begin
        @test objective_function_string(REPLMode, m) == "y + 2"
        @test objective_function_string(IJuliaMode, m) == "y + 2"
    end
end

# Helper function to test IO methods work correctly
function io_test(mode, obj, exp_str; repl=:both)
    if mode == REPLMode
        repl != :show  && @test sprint(print, obj) == exp_str
        repl != :print && @test sprint(show,  obj) == exp_str
    else
        @test sprint(show, "text/latex", obj) == string("\$\$ ",exp_str," \$\$")
    end
end

# test infinite model show
@testset "Show InfiniteModel" begin
    # initialize models
    m = InfiniteModel()
    @infinite_parameter(m, 0 <= par1 <= 1)
    @infinite_parameter(m, pars[1:2] in MvNormal([1, 1], 1))
    @infinite_variable(m, x(par1))
    @infinite_variable(m, z(pars))
    @global_variable(m, y)
    @objective(m, Min, 2 + y)
    @constraint(m, c1, x + y -2 <= 0)
    @constraint(m, c2, y^2 - 3 == 0)
    @constraint(m, c3, x == 5,
                parameter_bounds = Dict(par1 => IntervalSet(0, 0.5)))
    # test Base.show (constraint in REPL)
    @testset "Base.show (REPL Constraint)" begin
        # test normal
        str = "c1 : x(par1) + y " * JuMP._math_symbol(REPLMode, :leq) *
              " 2.0"
        io_test(REPLMode, c1, str)
        # test bounded
        str = "c3 : x(par1) " * JuMP._math_symbol(REPLMode, :eq) * " 5.0, " *
              JuMP._math_symbol(REPLMode, :for_all) * " par1 " *
              JuMP._math_symbol(REPLMode, :in) * " [0, 0.5]"
        io_test(REPLMode, c3, str)
    end
    # test Base.show (constraint in IJulia)
    @testset "Base.show (IJulia Constraint)" begin
        # test normal
        # str = "c1 : \$ x(par1) + y " * JuMP._math_symbol(IJuliaMode, :leq) *
        #       " 2.0 \$"
        # io_test(IJuliaMode, c1, str)
        # # test bounded
        # str =  "c3 : \$\$ x(par1) " * JuMP._math_symbol(IJuliaMode, :eq) * " 5.0, " *
        #        JuMP._math_symbol(IJuliaMode, :for_all) * " par1 " *
        #        JuMP._math_symbol(IJuliaMode, :in) * " [0, 0.5] \$\$"
        # io_test(IJuliaMode, c3, str)
    end
    # test show_backend_summary
    @testset "JuMP.show_backend_summary" begin

    end
end
