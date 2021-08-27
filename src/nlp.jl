################################################################################
#                                 DATATYPES
################################################################################
# Extend addchild to take the root of another graph as input
function _LCRST.addchild(parent::_LCRST.Node{T}, newc::_LCRST.Node{T}) where T
    # copy the new node if it is not a root
    # otherwise, we are just merging 2 graphs together
    if !_LCRST.isroot(newc)
        newc = copy(newc)
    end
    # add it on to the tree
    newc.parent = parent
    prevc = parent.child
    if prevc == parent
        parent.child = newc
    else
        prevc = _LCRST.lastsibling(prevc)
        prevc.sibling = newc
    end
    return newc
end

# Extend addchild with convenient nothing dispatch for empty previous child
function _LCRST.addchild(
    parent::_LCRST.Node{T}, 
    oldc::Nothing, 
    newc::_LCRST.Node{T}
    ) where T
    return _LCRST.addchild(parent, newc)
end

# Extend addchild to efficiently add multiple children if the previous is known
function _LCRST.addchild(
    parent::_LCRST.Node{T}, 
    prevc::_LCRST.Node{T}, 
    data::T
    ) where T
    # add it on to the tree
    newc = _LCRST.Node(data, parent)
    prevc.sibling = newc
    return newc
end

# Extend addchild to efficiently add multiple children if the previous is known
function _LCRST.addchild(
    parent::_LCRST.Node{T}, 
    prevc::_LCRST.Node{T}, 
    newc::_LCRST.Node{T}
    ) where T
    # copy the new node if it is not a root
    # otherwise, we are just merging 2 graphs together
    is_first = prevc === newc
    if !_LCRST.isroot(newc)
        newc = copy(newc)
    end
    # add it on to the tree
    newc.parent = parent
    if is_first
        parent.child = newc
    else
        prevc.sibling = newc
    end
    return newc
end

# Map a LCRST tree based by operating each node with a function
function _map_tree(map_func::Function, node::_LCRST.Node)
    new_node = map_func(node)
    prev = nothing
    for child in node
        prev = _LCRST.addchild(new_node, prev, _map_tree(map_func, child))
    end
    return new_node
end

# Extend copying for graph nodes
function Base.copy(node::_LCRST.Node)
    return _map_tree(n -> _LCRST.Node(n.data), node)
end

# This is ambiguous but faster than the concrete alternatives tested so far
# Even better than using Node{Any}...
"""

"""
struct NodeData
    value
end

# Getter function for the node value (so it is easy to change later on if needed)
function _node_value(data::NodeData)
    return data.value
end

"""

"""
struct NLPExpr <: JuMP.AbstractJuMPScalar
    expr::_LCRST.Node{NodeData}
end

# Extend basic functions
Base.broadcastable(nlp::NLPExpr) = Ref(nlp)
Base.copy(nlp::NLPExpr)::NLPExpr = NLPExpr(copy(nlp.expr))
Base.zero(::Type{NLPExpr}) = NLPExpr(_LCRST.Node(NodeData(0.0)))
Base.one(::Type{NLPExpr}) = NLPExpr(_LCRST.Node(NodeData(1.0)))

# Convenient expression alias
const AbstractInfOptExpr = Union{
    NLPExpr, 
    JuMP.GenericQuadExpr{Float64, GeneralVariableRef}, 
    JuMP.GenericAffExpr{Float64, GeneralVariableRef},
    GeneralVariableRef
}

# TODO create conversion function to AST 

# TODO make remove zeros function 

# TODO more intelligently operate with constants

################################################################################
#                          EXPRESSION CREATION HELPERS
################################################################################
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
function _call_graph(func::Symbol, arg1, args...)
    root = _LCRST.Node(NodeData(func))
    prevc = _LCRST.addchild(root, _process_child_input(arg1))
    for a in args 
        prevc = _LCRST.addchild(root, prevc, _process_child_input(a))
    end
    return root
end

################################################################################
#                               MUTABLE ARITHMETICS
################################################################################
_MA.mutability(::Type{NLPExpr}) = _MA.IsMutable()

function _MA.mutable_operate!(::typeof(_MA.add_mul), nlp::NLPExpr, args...) 
    return nlp + *(args...)
end

# TODO continue

################################################################################
#                               SUMS AND PRODUCTS
################################################################################
## Define helper functions for sum reductions 
# Container of NLPExprs
function _reduce_by_first(::typeof(sum), first_itr::NLPExpr, itr)
    root = _LCRST.Node(NodeData(:+))
    prevc = _LCRST.addchild(root, first_itr.expr, first_itr.expr)
    for ex in itr
        prevc = _LCRST.addchild(root, prevc, ex.expr)
    end
    return NLPExpr(root)
end

# Container of InfiniteOpt exprs
function _reduce_by_first(
    ::typeof(sum), 
    first_itr::JuMP.AbstractJuMPScalar, 
    itr
    )
    result = first_itr
    for i in itr 
        result = _MA.operate!(_MA.add_mul, result, i)
    end
    return result
end

# Fallback
function _reduce_by_first(::typeof(sum), first_itr, itr; kw...)
    return first_itr + sum(identity, itr; kw...)
end

# Hyjack Base.sum for better efficiency on iterators --> this is type piracy...
function Base.sum(itr::Base.Generator; kw...)
    isempty(itr) && return sum(identity, itr; kw...)
    itr1, new_itr = Iterators.peel(itr)
    return _reduce_by_first(sum, itr1, new_itr; kw...)
end

# Extend Base.sum for container of NLPExprs
function Base.sum(arr::AbstractArray{<:NLPExpr}; init = 0.0)
    isempty(arr) && return init
    itr1, new_itr = Iterators.peel(arr)
    return _reduce_by_first(sum, itr1, new_itr; kw...)
end

# Extend Base.sum for container of InfiniteOpt exprs
function Base.sum(arr::AbstractArray{<:AbstractInfOptExpr}; init = 0.0)
    isempty(arr) && return init
    result = _MA.Zero()
    for i in arr 
        result = _MA.operate!(_MA.add_mul, result, i)
    end
    return result
end

## Define helper functions for reducing products 
# Container of InfiniteOpt exprs
function _reduce_by_first(::typeof(prod), first_itr::AbstractInfOptExpr, itr)
    root = _LCRST.Node(NodeData(:*))
    prevc = _LCRST.addchild(root, first_itr.expr, first_itr.expr)
    for ex in itr
        prevc = _LCRST.addchild(root, prevc, ex.expr)
    end
    return NLPExpr(root)
end

# Fallback
function _reduce_by_first(::typeof(prod), first_itr, itr; kw...)
    return first_itr * prod(identity, itr; kw...)
end

# Hyjack Base.prod for better efficiency on iterators --> this is type piracy...
function Base.prod(itr::Base.Generator; kw...)
    isempty(itr) && return prod(identity, itr; kw...)
    itr1, new_itr = Iterators.peel(itr)
    return _reduce_by_first(prod, itr1, new_itr; kw...)
end

# Extend Base.prod for container of InfiniteOpt exprs
function Base.prod(arr::AbstractArray{<:AbstractInfOptExpr}; init = 0.0)
    isempty(arr) && return init
    itr1, new_itr = Iterators.peel(arr)
    return _reduce_by_first(prod, itr1, new_itr; kw...)
end

################################################################################
#                           MULTIPLICATION OPERATORS
################################################################################
function Base.:*(
    quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef},
    expr::AbstractInfOptExpr
    )
    return NLPExpr(_call_graph(:*, quad, expr))
end

function Base.:*(
    expr::AbstractInfOptExpr,
    quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef}
    )
    return NLPExpr(_call_graph(:*, expr, quad))
end

function Base.:*(
    quad1::JuMP.GenericQuadExpr{Float64, GeneralVariableRef},
    quad2::JuMP.GenericQuadExpr{Float64, GeneralVariableRef}
    )
    return NLPExpr(_call_graph(:*, quad1, quad2))
end

function Base.:*(
    nlp::NLPExpr,
    quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef}
    )
    return NLPExpr(_call_graph(:*, nlp, quad))
end

function Base.:*(
    quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef}, 
    nlp::NLPExpr
    )
    return NLPExpr(_call_graph(:*, quad, nlp))
end

function Base.:*(
    nlp::NLPExpr, 
    expr::Union{AbstractInfOptExpr, Real}
    )
    return NLPExpr(_call_graph(:*, nlp, expr))
end

function Base.:*(
    expr::Union{AbstractInfOptExpr, Real}, 
    nlp::NLPExpr
    )
    return NLPExpr(_call_graph(:*, expr, nlp))
end

function Base.:*(nlp1::NLPExpr, nlp2::NLPExpr)
    return NLPExpr(_call_graph(:*, nlp1, nlp2))
end

function Base.:*(
    expr1::AbstractInfOptExpr,
    expr2::AbstractInfOptExpr,
    expr3::AbstractInfOptExpr,
    exprs::Vararg{AbstractInfOptExpr}
    )
    return NLPExpr(_call_graph(:*, expr1, expr2, expr3, exprs...))
end

function Base.:*(nlp::NLPExpr)
    return nlp
end

################################################################################
#                              DIVISION OPERATORS
################################################################################
function Base.:/(
    expr1::Union{AbstractInfOptExpr, Real}, 
    expr2::AbstractInfOptExpr
    )
    return NLPExpr(_call_graph(:/, expr1, expr2))
end

function Base.:/(nlp::NLPExpr, c::Real)
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
    )
    return NLPExpr(_call_graph(:^, expr1, expr2))
end

################################################################################
#                             SUBTRACTION OPERATORS
################################################################################
function Base.:-(nlp::NLPExpr, expr::Union{AbstractInfOptExpr, Real})
    return NLPExpr(_call_graph(:-, nlp, expr))
end

function Base.:-(expr::Union{AbstractInfOptExpr, Real}, nlp::NLPExpr)
    return NLPExpr(_call_graph(:-, expr, nlp))
end

function Base.:-(nlp1::NLPExpr, nlp2::NLPExpr)
    return NLPExpr(_call_graph(:-, nlp1, nlp2))
end

function Base.:-(nlp::NLPExpr)
    return NLPExpr(_call_graph(:-, nlp))
end

################################################################################
#                              ADDITION OPERATORS
################################################################################
function Base.:+(nlp::NLPExpr, expr::Union{AbstractInfOptExpr, Real})
    return NLPExpr(_call_graph(:+, nlp, expr))
end

function Base.:+(expr::Union{AbstractInfOptExpr, Real}, nlp::NLPExpr)
    return NLPExpr(_call_graph(:+, expr, nlp))
end

function Base.:+(nlp1::NLPExpr, nlp2::NLPExpr)
    return NLPExpr(_call_graph(:+, nlp1, nlp2))
end

function Base.:+(nlp::NLPExpr)
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
        function Base.$name(v::AbstractInfOptExpr)
            return NLPExpr(_call_graph($(quot(name)), v))
        end
    end
end

# TODO add norm functions

# Setup the Base functions with 2 arguments
for (name, func) in (:min => min, :max => max)
    # add it to the main storage dict
    NativeNLPFunctions[name] = func
    # make an expression constructor
    @eval begin 
        function Base.$name(
            v1::Union{AbstractInfOptExpr, Real}, 
            v2::Union{AbstractInfOptExpr, Real}
            )
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
        function SpecialFunctions.$name(v::AbstractInfOptExpr)
            return NLPExpr(_call_graph($(quot(name)), v))
        end
    end
end

################################################################################
#                               USER FUNCTIONS
################################################################################
# TODO add function registration 

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
# Define better printing for NodeData
function Base.show(io::IO, data::NodeData)
    return print(io, string(_node_value(data)))
end

# Map operators to their respective precedence (largest is highest priority)
const _Precedence = (; :^ => 4, Symbol("+u") => 3, Symbol("-u") => 3, :* => 2, 
                     :/ => 2, :+ => 1, :- => 1)

## Make functions to determine the precedence of a leaf
# AffExpr
function _leaf_precedence(aff::JuMP.GenericAffExpr)
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
function _leaf_precedence(quad::JuMP.GenericQuadExpr)
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
function _leaf_precedence(v)
    return 10
end

# Recursively build an expression string, starting with a root node
function _expr_string(
    node::_LCRST.Node{NodeData}, 
    str::String = "";
    simple::Bool = false,
    prev_prec = 0,
    prev_comm = false
    )
    # prepocess the raw value
    raw_value = _node_value(node.data)
    is_op = !simple && raw_value in keys(_Precedence)
    data_str = _string_round(raw_value)
    # make a string according to the node structure
    if _LCRST.isleaf(node) && _leaf_precedence(raw_value) > prev_prec
        # we have a leaf that doesn't require parentheses
        return str * data_str
    elseif _LCRST.isleaf(node)
        # we have a leaf that requires parentheses
        return str * string("(", data_str, ")")
    elseif is_op && !_LCRST.islastsibling(node.child)
        # we have a binary operator
        curr_prec = _Precedence[raw_value]
        has_prec = curr_prec > prev_prec || (prev_comm && curr_prec == prev_prec)
        if !has_prec
            str *= "("
        end
        op_str = data_str == "^" ? data_str : string(" ", data_str, " ")
        is_comm = raw_value == :* || raw_value == :+
        for child in node
            str = string(_expr_string(child, str, prev_prec = curr_prec, 
                                      prev_comm = is_comm), op_str)
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
        str *= string(data_str, _expr_string(node.child, str, 
                                             prev_prec = curr_prec))
        return has_prec ? str : str * ")"
    else 
        # we have a function
        str *= string(data_str, "(")
        for child in node
            str = _expr_string(child, str, simple = simple)
            str *= ", "
        end
        return str[1:end-2] * ")"
    end
end

# Extend JuMP.function_string for nonlinear expressions
function JuMP.function_string(mode, nlp::NLPExpr)
    return _expr_string(nlp.expr)
end
