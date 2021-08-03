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
#                            NATIVE NLP FUNCTIONS
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

const _Operators = (:*, :-, :+, :^, :/)

function _expr_string(
    node::LCRST.Node{NodeData}, 
    str::String = "";
    simple::Bool = false,
    prev_op::Bool = false
    )::String
    # determine if the node is an operator 
    is_op = !simple && node.data.value in _Operators
    # make the node data into a string
    data_str = string(node.data.value)
    is_jump_expr = node.data.value isa Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr}
    # make a string according to the node structure
    if LCRST.isleaf(node) && prev_op && is_jump_expr
        # we have a leaf that is an affine/quadratic expr so we are done here
        return str * string("(", data_str, ")")
    elseif LCRST.isleaf(node)
        # we have a leaf that is a variable or number
        return str * data_str
    elseif is_op && !LCRST.islastsibling(node.child)
        # we have a binary operator
        str *= "("
        op_str = string(" ", data_str, " ")
        for child in node
            str = string(_expr_string(child, str, prev_op = true), op_str)
        end
        return str[1:end-length(op_str)] * ")"
    elseif is_op
        # we have a unary operator
        return string("(", data_str,
                      _expr_string(node.child, str, prev_op = true), ")")
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
