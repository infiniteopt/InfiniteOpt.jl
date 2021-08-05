################################################################################
#                                 DATATYPES
################################################################################
function Base.copy(node::LCRST.Node)
    new_node = LCRST.Node(node.data) # we don't make copies of the data (should we?)
    for child in node 
        LCRST.addchild(new_node, copy(child))
    end
    return new_node
end

function LCRST.addchild(parent::LCRST.Node{T}, newc::LCRST.Node{T}) where T
    # copy the new node if it is not a root
    # otherwise, we are just merging 2 graphs together
    if !LCRST.isroot(newc)
        newc = copy(newc)
    end
    # add it on to the tree
    newc.parent = parent
    prevc = parent.child
    if prevc == parent
        parent.child = newc
    else
        prevc = LCRST.lastsibling(prevc)
        prevc.sibling = newc
    end
    return newc
end

# This is ambiguous but faster than the concrete alternatives tested
# Even better than using Node{Any}...
"""

"""
struct NodeData
    value
end

"""

"""
struct NLPExpr <: JuMP.AbstractJuMPScalar
    expr::LCRST.Node{NodeData}
end

# Extend basic functions
Base.broadcastable(nlp::NLPExpr) = Ref(nlp)
Base.copy(nlp::NLPExpr)::NLPExpr = NLPExpr(copy(nlp.expr))
Base.zero(::Type{NLPExpr}) = NLPExpr(LCRST.Node(NodeData(0.0)))
Base.one(::Type{NLPExpr}) = NLPExpr(LCRST.Node(NodeData(1.0)))

# Convenient expression alias
const AbstractInfOptExpr = Union{
    NLPExpr, 
    JuMP.GenericQuadExpr{Float64, GeneralVariableRef}, 
    JuMP.GenericAffExpr{Float64, GeneralVariableRef},
    GeneralVariableRef
}

## Make convenient dispatch methods for raw child input
# NLPExpr
function _process_child_input(nlp::NLPExpr)
    return nlp.expr
end

# An InfiniteOpt expression (not general nonlinear)
function _process_child_input(v::AbstractInfOptExpr)
    return NodeData(v)
end

# Function symbol
function _process_child_input(f::Symbol)
    return NodeData(f)
end

# Real number 
function _process_child_input(c::Real)
    return NodeData(Float64(c))
end

# Fallback
function _process_child_input(v)
    error("Unrecognized algebraic expression input `$v`.")
end

# Generic graph builder
function _call_graph(func::Symbol, args...)::LCRST.Node{NodeData}
    root = LCRST.Node(NodeData(func))
    for a in args 
        LCRST.addchild(root, _process_child_input(a))
    end
    return root
end

################################################################################
#                           MULTIPLICATION OPERATORS
################################################################################
function Base.:*(
    quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef},
    expr::AbstractInfOptExpr
    )::NLPExpr
    return NLPExpr(_call_graph(:*, quad, expr))
end

function Base.:*(
    expr::AbstractInfOptExpr,
    quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef}
    )::NLPExpr
    return NLPExpr(_call_graph(:*, expr, quad))
end

function Base.:*(
    quad1::JuMP.GenericQuadExpr{Float64, GeneralVariableRef},
    quad2::JuMP.GenericQuadExpr{Float64, GeneralVariableRef}
    )::NLPExpr
    return NLPExpr(_call_graph(:*, quad1, quad2))
end

function Base.:*(
    nlp::NLPExpr, 
    expr::Union{AbstractInfOptExpr, Real}
    )::NLPExpr
    return NLPExpr(_call_graph(:*, nlp, expr))
end

function Base.:*(
    expr::Union{AbstractInfOptExpr, Real}, 
    nlp::NLPExpr
    )::NLPExpr
    return NLPExpr(_call_graph(:*, expr, nlp))
end

function Base.:*(nlp1::NLPExpr, nlp2::NLPExpr)::NLPExpr
    return NLPExpr(_call_graph(:*, nlp1, nlp2))
end

function Base.:*(
    expr1::AbstractInfOptExpr,
    expr2::AbstractInfOptExpr,
    expr3::AbstractInfOptExpr,
    exprs::Vararg{AbstractInfOptExpr}
    )::NLPExpr
    return NLPExpr(_call_graph(:*, expr1, expr2, expr3, exprs...))
end

function Base.:*(nlp::NLPExpr)::NLPExpr
    return nlp
end

################################################################################
#                              DIVISION OPERATORS
################################################################################
function Base.:/(
    expr1::Union{AbstractInfOptExpr, Real}, 
    expr2::AbstractInfOptExpr
    )::NLPExpr
    return NLPExpr(_call_graph(:/, expr1, expr2))
end

function Base.:/(nlp::NLPExpr, c::Real)::NLPExpr
    return NLPExpr(_call_graph(:/, nlp, c))
end

################################################################################
#                                POWER OPERATORS
################################################################################
function Base.:^(expr::AbstractInfOptExpr, c::Integer)
    if iszero(c)
        return 1.0
    elseif isone(c)
        return expr 
    elseif c == 2 
        return expr * expr
    else 
        return NLPExpr(_call_graph(:^, expr, c))
    end
end

function Base.:^(expr::AbstractInfOptExpr, c::Real)
    if iszero(c)
        return 1.0
    elseif isone(c)
        return expr 
    elseif c == 2 
        return expr * expr
    else 
        return NLPExpr(_call_graph(:^, expr, c))
    end
end

function Base.:^(
    expr1::Union{AbstractInfOptExpr, Real}, 
    expr2::AbstractInfOptExpr
    )::NLPExpr
    return NLPExpr(_call_graph(:^, expr1, expr2))
end

################################################################################
#                             SUBTRACTION OPERATORS
################################################################################
function Base.:-(nlp::NLPExpr, expr::Union{AbstractInfOptExpr, Real})::NLPExpr
    return NLPExpr(_call_graph(:-, nlp, expr))
end

function Base.:-(expr::Union{AbstractInfOptExpr, Real}, nlp::NLPExpr)::NLPExpr
    return NLPExpr(_call_graph(:-, expr, nlp))
end

function Base.:-(nlp1::NLPExpr, nlp2::NLPExpr)::NLPExpr
    return NLPExpr(_call_graph(:-, nlp1, nlp2))
end

function Base.:-(nlp::NLPExpr)::NLPExpr
    return NLPExpr(_call_graph(:-, nlp))
end

################################################################################
#                              ADDITION OPERATORS
################################################################################
function Base.:+(nlp::NLPExpr, expr::Union{AbstractInfOptExpr, Real})::NLPExpr
    return NLPExpr(_call_graph(:+, nlp, expr))
end

function Base.:+(expr::Union{AbstractInfOptExpr, Real}, nlp::NLPExpr)::NLPExpr
    return NLPExpr(_call_graph(:+, expr, nlp))
end

function Base.:+(nlp1::NLPExpr, nlp2::NLPExpr)::NLPExpr
    return NLPExpr(_call_graph(:+, nlp1, nlp2))
end

function Base.:+(nlp::NLPExpr)::NLPExpr
    return NLPExpr(_call_graph(:+, nlp))
end

################################################################################
#                             NATIVE NLP FUNCTIONS
################################################################################
"""

"""
const NativeNLPFunctions = Dict{Symbol, Function}()

"""

"""
const Base1ArgFuncList = (
    :sqrt => sqrt,
    :cbrt => cbrt,
    :abs => abs,
    :abs2 => abs2,
    :inv => inv,
    :log => log,
    :log10 => log10,
    :log2 => log2,
    :log1p => log1p,
    :exp => exp,
    :exp2 => exp2,
    :expm1 => expm1,
    :sin => sin,
    :cos => cos,
    :tan => tan, 
    :sec => sec, 
    :csc => csc,
    :cot => cot,
    :sind => sind,
    :cosd => cosd,
    :tand => tand, 
    :secd => secd,
    :cscd => cscd,
    :cotd => cotd,
    :asin => asin,
    :acos => acos, 
    :atan => atan, 
    :asec => asec, 
    :acsc => acsc, 
    :acot => acot,
    :asind => asind,
    :acosd => acosd,
    :atand => atand,
    :asecd => asecd,
    :acscd => acscd,
    :acotd => acotd, 
    :sinh => sinh,
    :cosh => cosh, 
    :tanh => tanh, 
    :sech => sech,
    :csch => csch,
    :coth => coth,
    :asinh => asinh, 
    :acosh => acosh,
    :atanh => atanh,
    :asech => asech,
    :acsch => acsch,
    :acoth => acoth,
    :deg2rad => deg2rad,
    :rad2deg => rad2deg
)

# Setup the base 1 argument functions
for (name, func) in Base1ArgFuncList
    # add it to the main storage dict
    NativeNLPFunctions[name] = func
    # make an expression constructor
    @eval begin 
        function Base.$name(v::AbstractInfOptExpr)::NLPExpr
            return NLPExpr(_call_graph($(quot(name)), v))
        end
    end
end

# Setup the Base functions with 2 arguments
for (name, func) in (:min => min, :max => max)
    # add it to the main storage dict
    NativeNLPFunctions[name] = func
    # make an expression constructor
    @eval begin 
        function Base.$name(
            v1::Union{AbstractInfOptExpr, Real}, 
            v2::Union{AbstractInfOptExpr, Real}
            )::NLPExpr
            return NLPExpr(_call_graph($(quot(name)), v1, v2))
        end
    end
end

# TODO figure out how to do the logical functions properly

# Setup the ifelse function
# NativeNLPFunctions[:ifelse] = Base.ifelse
# function ifelse(
#     cond::NLPExpr,
#     v1::Union{AbstractInfOptExpr, Real}, 
#     v2::Union{AbstractInfOptExpr, Real}
#     )::NLPExpr
#     return NLPExpr(_call_graph(:ifelse, cond, v1, v2))
# end

# Setup the Base comparison functions
# for (name, func) in (:< => Base.:(<), :(==) => Base.:(==))
#     # add it to the main storage dict
#     NativeNLPFunctions[name] = func
#     # make an expression constructor
#     @eval begin 
#         function Base.$name(
#             v1::AbstractInfOptExpr, 
#             v2::Union{AbstractInfOptExpr, Real}
#             )::NLPExpr
#             return NLPExpr(_call_graph($(quot(name)), v1, v2))
#         end
#     end
# end

"""

"""
const Special1ArgFuncList = (
    :erf => SpecialFunctions.erf,
    :erfinv => SpecialFunctions.erfinv,
    :erfc => SpecialFunctions.erfc,
    :erfcinv => SpecialFunctions.erfcinv,
    :erfi => SpecialFunctions.erfi,
    :gamma => SpecialFunctions.gamma,
    :lgamma => SpecialFunctions.lgamma,
    :digamma => SpecialFunctions.digamma,
    :invdigamma => SpecialFunctions.invdigamma,
    :trigamma => SpecialFunctions.trigamma,
    :airyai => SpecialFunctions.airyai,
    :airybi => SpecialFunctions.airybi,
    :airyaiprime => SpecialFunctions.airyaiprime,
    :airybiprime => SpecialFunctions.airybiprime,
    :besselj0 => SpecialFunctions.besselj0,
    :besselj1 => SpecialFunctions.besselj1,
    :bessely0 => SpecialFunctions.bessely0,
    :bessely1 => SpecialFunctions.bessely1,
    :erfcx => SpecialFunctions.erfcx,
    :dawson => SpecialFunctions.dawson
)

# Setup the SpecialFunctions 1 argument functions
for (name, func) in Special1ArgFuncList
    # add it to the main storage dict
    NativeNLPFunctions[name] = func
    # make an expression constructor
    @eval begin 
        function SpecialFunctions.$name(v::AbstractInfOptExpr)::NLPExpr
            return NLPExpr(_call_graph($(quot(name)), v))
        end
    end
end

################################################################################
#                               USER FUNCTIONS
################################################################################


################################################################################
#                               LINEAR ALGEBRA
################################################################################
LinearAlgebra.dot(lhs::AbstractInfOptExpr, rhs::AbstractInfOptExpr) = lhs * rhs
LinearAlgebra.dot(lhs::AbstractInfOptExpr, rhs::Real) = lhs * rhs
LinearAlgebra.dot(lhs::Real, rhs::AbstractInfOptExpr) = lhs * rhs

# TODO figure out what needs to be added to make this work (involves MutableArithmetics)

################################################################################
#                                  PRINTING
################################################################################
function Base.show(io::IO, data::NodeData)
    return print(io, string(data.value))
end

const _Precedence = (; :^ => 4, Symbol("+u") => 3, Symbol("-u") => 3, :* => 2, 
                     :/ => 2, :+ => 1, :- => 1)

## Make functions to determine the precedence of a leaf
# AffExpr
function _leaf_precedence(aff::JuMP.GenericAffExpr)::Int
    has_const = !iszero(JuMP.constant(aff))
    itr = JuMP.linear_terms(aff)
    num_terms = length(itr) 
    if iszero(num_terms)
        # we have only a constant
        return 10 # will always have precedence
    elseif has_const || num_terms > 1
        # we have an expr with multiple terms
        return 1
    elseif isone(first(itr)[1])
        # we have a single variable
        return 10 # will always have precedence
    elseif first(itr)[1] == -1
        # we have a single unary negative variable
        return 3
    else
        # we have a single variable multiplied by some coefficient
        return 2
    end
end

# QuadExpr
function _leaf_precedence(quad::JuMP.GenericQuadExpr)::Int
    has_aff = !iszero(quad.aff)
    itr = JuMP.quad_terms(quad)
    num_terms = length(itr) 
    if iszero(num_terms)
        # we have an affine expression
        return _leaf_precedence(quad.aff)
    elseif has_aff || num_terms > 1
        # we have a general quadratic expression
        return 1
    else
        # we only have a single quadratic term
        return 2
    end
end

# Other
function _leaf_precedence(v)::Int
    return 10
end

# Recursively build an expression string, starting with a root node
function _expr_string(
    node::LCRST.Node{NodeData}, 
    str::String = "";
    simple::Bool = false,
    prev_prec = 0
    )::String
    # prepocess the raw value
    raw_value = node.data.value
    is_op = !simple && raw_value in keys(_Precedence)
    data_str = string(raw_value)
    # make a string according to the node structure
    if LCRST.isleaf(node) && _leaf_precedence(raw_value) > prev_prec
        # we have a leaf that doesn't require parentheses
        return str * data_str
    elseif LCRST.isleaf(node)
        # we have a leaf that requires parentheses
        return str * string("(", data_str, ")")
    elseif is_op && !LCRST.islastsibling(node.child)
        # we have a binary operator
        curr_prec = _Precedence[raw_value]
        has_prec = curr_prec > prev_prec
        if !has_prec
            str *= "("
        end
        op_str = string(" ", data_str, " ")
        for child in node
            str = string(_expr_string(child, str, prev_prec = curr_prec), op_str)
        end
        str = str[1:end-length(op_str)]
        return has_prec ? str : str * ")"
    elseif is_op
        # we have a unary operator
        curr_prec = _Precedence[Symbol(raw_value, :u)]
        has_prec = curr_prec > prev_prec
        if !has_prec 
            str *= "("
        end
        str *= string(data_str, _expr_string(node.child, str, prev_prec = curr_prec))
        return has_prec ? str : str * ")"
    else 
        # we have a function with 1 or 2 inputs
        str *= string(data_str, "(")
        for child in node
            str = _expr_string(child, str, simple = simple)
            str *= ", "
        end
        return str[1:end-2] * ")"
    end
end

# Extend JuMP.function_string for nonlinear expressions
function JuMP.function_string(mode, nlp::NLPExpr)::String
    return _expr_string(nlp.expr)
end
