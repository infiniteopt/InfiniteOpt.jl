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

# Make macro for container construction testing 
macro gen_container(expr, vals)
    idxvars, inds = InfiniteOpt._build_ref_sets(error, expr)
    code = JuMPC.container_code(idxvars, inds, esc(vals), :Auto)
    return code
end

# Define test data structures
struct BadDomain <: AbstractInfiniteDomain end
struct BadScalarDomain <: InfiniteScalarDomain end
struct BadArrayDomain <: InfiniteArrayDomain end
struct TestBridge{C} <: MOI.Bridges.AbstractBridge where {C} end
struct TestMethod <: NonGenerativeDerivativeMethod end
struct TestGenMethod <: GenerativeDerivativeMethod end
struct BadFiniteTech <: FDTechnique end
struct TestGenInfo <: AbstractGenerativeInfo end
struct BadData <: AbstractMeasureData end
struct Bad end
struct NotADomainType end
struct TestJuMPTag <: AbstractJuMPTag end
struct TestBackend <: AbstractTransformationBackend end
struct TestIndex <: ObjectIndex
    value::Int
end
struct TestVariableRef <: DispatchVariableRef
    model::InfiniteModel
    index::TestIndex
end
InfiniteOpt.dispatch_variable_ref(m::InfiniteModel, i::TestIndex) = TestVariableRef(m, i)

struct TestIndex2 <: ObjectIndex
    value::Int
end
struct TestVariableRef2 <: DispatchVariableRef
    model::InfiniteModel
    index::TestIndex2
end
InfiniteOpt.dispatch_variable_ref(m::InfiniteModel, i::TestIndex2) = TestVariableRef2(m, i)
JuMP.name(::TestVariableRef2) = "test"

InfiniteOpt.support_label(::TestGenMethod) = InternalLabel

# Define test functions
function new_fn end

# Define some fallbacks for TestBackend
JuMP.get_attribute(::TestBackend, ::MOI.TerminationStatus) = MOI.ALMOST_LOCALLY_SOLVED
JuMP.get_attribute(::TestBackend, ::MOI.PrimalStatus) = MOI.NEARLY_FEASIBLE_POINT
JuMP.get_attribute(::TestBackend, ::MOI.DualStatus) = MOI.NEARLY_FEASIBLE_POINT

# Helper function to test IO methods work correctly
function show_test(mode, obj, exp_str::String; repl=:both)
    if mode == MIME("text/plain")
        repl != :show  && @test sprint(print, obj) == exp_str
        repl != :print && @test sprint(show,  obj) == exp_str
    else
        @test sprint(show, "text/latex", obj) == exp_str
    end
end

# Helper function for IO methods with different possibilities
function show_test(mode, obj, exp_str::Vector{String}; repl=:both)
    if mode == MIME("text/plain")
        repl != :show  && @test sprint(print, obj) in exp_str
        repl != :print && @test sprint(show,  obj) in exp_str
    else
        @test sprint(show, "text/latex", obj) in exp_str
    end
end

# Make another helper function for other io methods
io_test(f::Function, exp_str::String, args...) = begin
    io = IOBuffer()
    f(io, args...)
    @test String(take!(io)) == exp_str
end

# Make another helper function for other io methods with multiple options
io_test(f::Function, exp_str::Vector{String}, args...) = begin
    io = IOBuffer()
    f(io, args...)
    @test String(take!(io)) in exp_str
end

# Test the output of a function that prints to stdout
function stdout_test(f::Function, exp_str, args...)
    original_stdout = stdout
    (read_pipe, write_pipe) = redirect_stdout()
    @test f(args...) isa Nothing
    redirect_stdout(original_stdout)
    close(write_pipe)
    @test read(read_pipe, String) == exp_str
end

# Make method for sorting matrices
sortcols(A) = sortslices(A, dims=2, lt=(x,y)->isless(x[:],y[:]))

# Extend comparison of SparseAxisArrays
function Base.:(==)(a::JuMPC.SparseAxisArray, b::JuMPC.SparseAxisArray)::Bool
    return a.data == b.data
end

## Make dumby measure data for the purpose of testing finite variables
struct TestData{P, A} <: AbstractMeasureData
    pref::P
    lb::A
    ub::A
end
InfiniteOpt.parameter_refs(d::TestData) = d.pref
function InfiniteOpt.add_supports_to_parameters(d::TestData)::Nothing
    return
end
InfiniteOpt.support_label(d::TestData) = UniqueMeasure{:a}


function Base.isequal(nlp1::GenericNonlinearExpr, nlp2::GenericNonlinearExpr)
    return nlp1.head == nlp2.head && isequal(nlp1.args, nlp2.args)
end


