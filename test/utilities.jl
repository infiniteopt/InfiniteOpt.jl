using JuMP: REPLMode, IJuliaMode

macro test_expression(expr)
    esc(quote
            @test JuMP.isequal_canonical(@expression(m, $expr), $expr)
    end)
end

macro test_expression_with_string(expr, str)
    esc(quote
            @test string(@inferred $expr) == $str
            @test_expression $expr
    end)
end

# Test that the macro call `m` throws an error exception during pre-compilation
macro test_macro_throws(errortype, m)
    # See https://discourse.julialang.org/t/test-throws-with-macros-after-pr-23533/5878
    :(@test_throws $errortype try @eval $m catch err; throw(err.error) end)
end

# Define test data structures
struct BadSet <: AbstractInfiniteSet end
struct BadScalarSet <: InfiniteScalarSet end
struct BadArraySet <: InfiniteArraySet end
struct TestBridge{C} <: MOI.Bridges.AbstractBridge where {C} end
struct BadData <: AbstractMeasureData end
struct Bad end
struct NotASetType end
struct TestIndex <: ObjectIndex
    value::Int
end
struct TestVariableRef <: DispatchVariableRef
    model::InfiniteModel
    index::TestIndex
end
InfiniteOpt.dispatch_variable_ref(m::InfiniteModel, i::TestIndex) = TestVariableRef(m, i)

# Define test functions
function new_fn end

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

# Make method for sorting matrices
sortcols(A) = sortslices(A, dims=2, lt=(x,y)->isless(x[:],y[:]))
