################################################################################
#                                 DATATYPES
################################################################################
# This is ambiguous but faster than the concrete alternatives tested
# Even better than using BinaryNode{Any}...
struct NodeData
    value
end

struct NLPExpr <: JuMP.AbstractJuMPScalar
    expr::Collections.BinaryNode{NodeData}
end

Base.broadcastable(nlp::NLPExpr) = Ref(nlp)

Base.copy(nlp::NLPExpr)::NLPExpr = NLPExpr(copy(nlp.expr))

const AbstractInfOptExpr = Union{
        NLPExpr, 
        JuMP.GenericQuadExpr{Float64, GeneralVariableRef}, 
        JuMP.GenericAffExpr{Float64, GeneralVariableRef},
        GeneralVariableRef
    }

function _process_child_input(nlp::NLPExpr)
    # copy if the expression root has a parent 
    # meaning it was already embedded in building an expression
    # this helps only copy if necessary
    if isdefined(nlp.expr, :parent)
        return copy(nlp.expr)
    else
        return nlp.expr
    end
end

function _process_child_input(v)
    return NodeData(v)
end

function _series_graph(op::Symbol, v1, v2, vars...)::Collections.BinaryNode{NodeData}
    root = Collections.BinaryNode(NodeData(op))
    Collections.rightchild(_process_child_input(vars[end]), root)
    left = Collections.leftchild(NodeData(op), root)
    for v in Iterators.rest(Iterators.reverse(vars), 2)
        Collections.rightchild(_process_child_input(v), left)
        left = Collections.leftchild(NodeData(op), left)
    end
    Collections.rightchild(_process_child_input(v2), left)
    Collections.leftchild(_process_child_input(v1), left)
    return root
end

function _series_graph(op::Symbol, v1, v2, v3)::Collections.BinaryNode{NodeData}
    root = Collections.BinaryNode(NodeData(op))
    Collections.rightchild(_process_child_input(v3), root)
    left = Collections.leftchild(NodeData(op), root)
    Collections.rightchild(_process_child_input(v2), left)
    Collections.leftchild(_process_child_input(v1), left)
    return root
end

function _series_graph(op::Symbol, v1, v2)::Collections.BinaryNode{NodeData}
    root = Collections.BinaryNode(NodeData(op))
    Collections.leftchild(_process_child_input(v1), root)
    Collections.rightchild(_process_child_input(v2), root)
    return root
end

################################################################################
#                           MULTIPLICATION OPERATORS
################################################################################
function Base.:*(
    quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef},
    expr::AbstractInfOptExpr
    )::NLPExpr
    return NLPExpr(_series_graph(:*, quad, expr))
end

function Base.:*(
    expr::AbstractInfOptExpr,
    quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef}
    )::NLPExpr
    return NLPExpr(_series_graph(:*, expr, quad))
end

function Base.:*(
    quad1::JuMP.GenericQuadExpr{Float64, GeneralVariableRef},
    quad2::JuMP.GenericQuadExpr{Float64, GeneralVariableRef}
    )::NLPExpr
    return NLPExpr(_series_graph(:*, quad1, quad2))
end

function Base.:*(
    nlp::NLPExpr, 
    expr::Union{AbstractInfOptExpr, Real}
    )::NLPExpr
    return NLPExpr(_series_graph(:*, nlp, expr))
end

function Base.:*(
    expr::Union{AbstractInfOptExpr, Real}, 
    nlp::NLPExpr
    )::NLPExpr
    return NLPExpr(_series_graph(:*, expr, nlp))
end

function Base.:*(nlp1::NLPExpr, nlp2::NLPExpr)::NLPExpr
    return NLPExpr(_series_graph(:*, nlp1, nlp2))
end

function Base.:*(
    expr1::AbstractInfOptExpr,
    expr2::AbstractInfOptExpr,
    expr3::AbstractInfOptExpr,
    exprs::Vararg{AbstractInfOptExpr}
    )::NLPExpr
    return NLPExpr(_series_graph(:*, expr1, expr2, expr3, exprs...))
end

################################################################################
#                              DIVISION OPERATORS
################################################################################
function Base.:/(
    expr1::Union{AbstractInfOptExpr, Real}, 
    expr2::AbstractInfOptExpr
    )::NLPExpr
    return NLPExpr(_series_graph(:/, expr1, expr2))
end

function Base.:/(nlp::NLPExpr, c::Real)::NLPExpr
    return NLPExpr(_series_graph(:/, nlp, c))
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
        return NLPExpr(_series_graph(:^, expr, c))
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
        return NLPExpr(_series_graph(:^, expr, c))
    end
end

function Base.:^(
    expr1::Union{AbstractInfOptExpr, Real}, 
    expr2::AbstractInfOptExpr
    )::NLPExpr
    return NLPExpr(_series_graph(:^, expr1, expr2))
end

################################################################################
#                             SUBTRACTION OPERATORS
################################################################################
function Base.:-(nlp::NLPExpr, expr::Union{AbstractInfOptExpr, Real})::NLPExpr
    return NLPExpr(_series_graph(:-, nlp, expr))
end

function Base.:-(expr::Union{AbstractInfOptExpr, Real}, nlp::NLPExpr)::NLPExpr
    return NLPExpr(_series_graph(:-, expr, nlp))
end

function Base.:-(nlp1::NLPExpr, nlp2::NLPExpr)::NLPExpr
    return NLPExpr(_series_graph(:-, nlp1, nlp2))
end

################################################################################
#                              ADDITION OPERATORS
################################################################################
function Base.:+(nlp::NLPExpr, expr::Union{AbstractInfOptExpr, Real})::NLPExpr
    return NLPExpr(_series_graph(:+, nlp, expr))
end

function Base.:+(expr::Union{AbstractInfOptExpr, Real}, nlp::NLPExpr)::NLPExpr
    return NLPExpr(_series_graph(:+, expr, nlp))
end

function Base.:+(nlp1::NLPExpr, nlp2::NLPExpr)::NLPExpr
    return NLPExpr(_series_graph(:+, nlp1, nlp2))
end

################################################################################
#                                  PRINTING
################################################################################
function Base.show(io::IO, data::NodeData)
    return print(io, string(data.value))
end

const _Operators = (:*, :-, :+, :^, :/)

function _expr_string(
    root::Collections.BinaryNode{NodeData}, 
    str::String = "";
    simple::Bool = false,
    prev_op::Bool = false
    )
    # determine if the node is an operator 
    is_op = !simple && root.data.value in _Operators
    # make the node data into a string
    data_str = string(root.data.value)
    is_jump_expr = root.data.value isa Union{JuMP.GenericAffExpr, JuMP.GenericQuadExpr}
    # make a string according to the node structure
    if !Collections.has_children(root) && prev_op && is_jump_expr
        # we have a leaf that is an affine/quadratic expr so we are done here
        return str * string("(", data_str, ")")
    elseif !Collections.has_children(root)
        # we have a leaf that is a variable or number
        return str * data_str
    elseif is_op && isdefined(root, :right)
        # we have a binary operator
        return string("(", _expr_string(root.left, str, prev_op = true), " ", 
                      data_str, " ", 
                      _expr_string(root.right, str, prev_op = true), ")")
    elseif is_op
        # we have a unary operator
        return string("(", data_str,
                      _expr_string(root.left, str, prev_op = true), ")")
    else 
        # we have a function with 1 or 2 inputs
        str *= string(data_str, "(")
        str = _expr_string(root.left, str, simple = simple)
        if isdefined(root, :right)
            str *= ", "
            str = _expr_string(root.right, str, simple = simple)
        end
        return str * ")"
    end
end

# Extend JuMP.function_string for nonlinear expressions
function JuMP.function_string(mode, nlp::NLPExpr)::String
    return _expr_string(nlp.expr)
end
