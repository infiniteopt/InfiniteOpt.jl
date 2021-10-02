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
    # check if the prev is actually a child of the parent 
    @assert prevc.parent === parent "Previous child doesn't belong to parent."
    # copy the new node if it is not a root
    # otherwise, we are just merging 2 graphs together
    if !_LCRST.isroot(newc)
        newc = copy(newc)
    end
    # add it on to the tree
    newc.parent = parent
    prevc.sibling = newc
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

# Replace a node with its only child if it only has 1 child
function _merge_parent_and_child(node::_LCRST.Node)
    if _LCRST.islastsibling(node.child)
        child = node.child
        node.data = child.data
        for n in child
            n.parent = node
        end
        node.child = child.child
        child.child = child
        child.parent = child
    end
    return node
end

# This is ambiguous but faster than the concrete alternatives tested so far
# Even better than using Node{Any}...
"""
    NoteData

A `DataType` for storing values in an expression tree that is used in a 
[`NLPExpr`](@ref). Acceptable value types include:
- `Real`: Constants
- `GeneralVariableRef`: Optimization variables
- `JuMP.GenericAffExpr{Float64, GeneralVariableRef}`: Affine expressions
- `JuMP.GenericQuadExpr{Float64, GeneralVariableRef}`: Quadratic expressions
- `Symbol`: Registered NLP function name.

**Fields**
- `value`: The stored value.
"""
struct NodeData
    value
end

# Getter function for the node value (so it is easy to change later on if needed)
function _node_value(data::NodeData)
    return data.value
end

# Recursively determine if node is effectively zero
function _is_zero(node::_LCRST.Node{NodeData})
    raw = _node_value(node.data)
    if isequal(raw, 0)
        return true
    elseif _LCRST.isleaf(node)
        return false
    elseif raw in (:+, :-) && all(_is_zero(n) for n in node)
        return true
    elseif raw == :* && any(_is_zero(n) for n in node)
        return true
    elseif raw in (:/, :^) && _is_zero(node.child)
        return true
    elseif all(_is_zero(n) for n in node) && iszero(get(_NativeNLPFunctions, raw, (i...) -> true)((0.0 for n in node)...))
        return true
    else
        return false
    end
end

# Prone any nodes that are effectively zero
function _drop_zeros!(node::_LCRST.Node{NodeData})
    if _LCRST.isleaf(node)
        return node
    elseif _is_zero(node)
        node.data = NodeData(0.0)
        _LCRST.makeleaf!(node)
        return node
    end
    raw = _node_value(node.data)
    if raw == :+
        for n in node
            if _is_zero(n)
                _LCRST.prunebranch!(n)
            end
        end
        _merge_parent_and_child(node)
    elseif raw == :-
        if _is_zero(node.child)
            _LCRST.prunebranch!(node.child)
        end
    end
    for n in node 
        _drop_zeros!(n)
    end
    return node
end

# Extend Base.isequal for our node types
function Base.isequal(n1::_LCRST.Node{NodeData}, n2::_LCRST.Node{NodeData})
    isequal(_node_value(n1.data), _node_value(n2.data)) || return false
    count(i -> true, n1) != count(i -> true, n2) && return false
    for (c1, c2) in zip(n1, n2)
        if !isequal(c1, c2)
            return false
        end
    end
    return true
end

"""
    NLPExpr <: JuMP.AbstractJuMPScalar

A `DataType` for storing scalar nonlinear expressions. It stores the expression 
algebraically via an expression tree where each node contains [`NodeData`](@ref) 
that can store one of the following:
- a registered function name (stored as a `Symbol`)
- a constant
- a variable 
- an affine expression
- a quadratic expression.
Specifically, it employs a left-child right-sibling tree 
(from `LeftChildRightSiblingTrees.jl`) to represent the expression tree.

**Fields**
- `tree_root::LeftChildRightSiblingTrees.Node{NodeData}`: The root node of the 
  expression tree.
"""
struct NLPExpr <: JuMP.AbstractJuMPScalar
    tree_root::_LCRST.Node{NodeData}

    # Constructor
    function NLPExpr(tree_root::_LCRST.Node{NodeData})
        @warn "General nonlinear expression support is experimental. Please, " *
              "notify us on GitHub if you run into unexpected behavior." maxlog = 1
        return new(tree_root)
    end
end

# Extend basic functions
Base.broadcastable(nlp::NLPExpr) = Ref(nlp)
Base.copy(nlp::NLPExpr) = NLPExpr(copy(nlp.tree_root))
Base.zero(::Type{NLPExpr}) = NLPExpr(_LCRST.Node(NodeData(0.0)))
Base.one(::Type{NLPExpr}) = NLPExpr(_LCRST.Node(NodeData(1.0)))
function Base.isequal(nlp1::NLPExpr, nlp2::NLPExpr) 
    return isequal(nlp1.tree_root, nlp2.tree_root)
end

"""
    JuMP.drop_zeros!(nlp::NLPExpr)::NLPExpr

Removes the zeros (possibly introduced by deletion) from an nonlinear expression. 
Note this only uses a few simple heuristics and will not remove more complex 
relationships like `cos(π/2)`. 

**Example**
```julia-repl
julia> expr = x^2.3 * max(0, zero(NLPExpr)) - exp(1/x + 0)
x^2.3 * max(0, 0) - exp(1 / x + 0)

julia> drop_zeros!(expr)
-exp(1 / x)
```
"""
function JuMP.drop_zeros!(nlp::NLPExpr)
    _drop_zeros!(nlp.tree_root) # uses a basic simplification scheme
    return nlp
end

# Extend JuMP.isequal_canonical (uses some heuristics but is not perfect)
function JuMP.isequal_canonical(nlp1::NLPExpr, nlp2::NLPExpr)
    n1 = _drop_zeros!(copy(nlp1.tree_root))
    n2 = _drop_zeros!(copy(nlp2.tree_root))
    return isequal(n1, n2)
end

# Print the tree structure of the expression tree
function print_expression_tree(io::IO, nlp::NLPExpr) 
    return AbstractTrees.print_tree(io, nlp.tree_root)
end
print_expression_tree(io::IO, expr) = println(io, expr)

"""
    print_expression_tree(nlp::NLPExpr)

Print a tree representation of the nonlinear expression `nlp`.

**Example**
```julia-repl
julia> expr = (x * sin(x)^3) / 2
(x * sin(x)^3) / 2

julia> print_expression_tree(expr)
/
├─ *
│  ├─ x
│  └─ ^
│     ├─ sin
│     │  └─ x
│     └─ 3
└─ 2
```
"""
print_expression_tree(nlp::NLPExpr) = print_expression_tree(stdout::IO, nlp) 

# Convenient expression alias
const AbstractInfOptExpr = Union{
    NLPExpr, 
    JuMP.GenericQuadExpr{Float64, GeneralVariableRef}, 
    JuMP.GenericAffExpr{Float64, GeneralVariableRef},
    GeneralVariableRef
}

## Dispatch function for ast mapping 
# Constant 
function _ast_process_node(map_func::Function, c)
    return c
end

# Variable
function _ast_process_node(map_func::Function, v::GeneralVariableRef)
    return map_func(v)
end

# AffExpr
function _ast_process_node(map_func::Function, aff::JuMP.GenericAffExpr)
    ex = Expr(:call, :+)
    if !iszero(aff.constant)
        push!(ex.args, aff.constant)
    end
    for (v, c) in aff.terms
        if isone(c)
            push!(ex.args, map_func(v))
        else
            push!(ex.args, Expr(:call, :*, c, map_func(v)))
        end
    end
    return ex
end

# QuadExpr
function _ast_process_node(map_func::Function, quad::JuMP.GenericQuadExpr)
    ex = Expr(:call, :+)
    append!(ex.args, _ast_process_node(map_func, quad.aff).args[2:end])
    for (xy, c) in quad.terms
        if isone(c)
            push!(ex.args, Expr(:call, :*, map_func(xy.a), map_func(xy.b)))
        else
            push!(ex.args, Expr(:call, :*, c, map_func(xy.a), map_func(xy.b)))
        end
    end
    return ex
end

# Map an expression tree to a Julia AST tree that is compatible with JuMP
function _tree_map_to_ast(map_func::Function, node::_LCRST.Node)    
    if _LCRST.isleaf(node)
        return _ast_process_node(map_func, _node_value(node.data))
    else
        ex = Expr(:call, _node_value(node.data)) # will be function symbol name 
        append!(ex.args, (_tree_map_to_ast(map_func, n) for n in node))
        return ex
    end
end

"""
    map_nlp_to_ast(map_func::Function, nlp::NLPExpr)::Expr

Map the nonlinear expression `nlp` to a Julia AST expression where each variable 
is mapped via `map_func` and is directly interpolated into the AST expression. 
This is intended as an internal method that can be helpful for developers that 
wish to map a `NLPExpr` to a Julia AST expression that is compatible with 
`JuMP.add_NL_expression`. 
"""
function map_nlp_to_ast(map_func::Function, nlp::NLPExpr)
    return _tree_map_to_ast(map_func, nlp.tree_root)
end

################################################################################
#                          EXPRESSION CREATION HELPERS
################################################################################
## Make convenient dispatch methods for raw child input
# NLPExpr
function _process_child_input(nlp::NLPExpr)
    return nlp.tree_root
end

# An InfiniteOpt expression (not general nonlinear)
function _process_child_input(v::AbstractInfOptExpr)
    return NodeData(v)
end

# Function symbol
function _process_child_input(f::Symbol)
    return NodeData(f)
end

# A constant
function _process_child_input(c::Union{Real, Bool})
    return NodeData(c)
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
# Define NLPExpr as a mutable type for MA
_MA.mutability(::Type{NLPExpr}) = _MA.IsMutable()

# Extend MA.promote_operation for bettered efficiency (TODO fix ambiguousnous)
for type in (:Real, :GeneralVariableRef,
             :(JuMP.GenericAffExpr{Float64, GeneralVariableRef}), 
             :(JuMP.GenericQuadExpr{Float64, GeneralVariableRef}))
    @eval begin
        function _MA.promote_operation(
            ::Union{typeof(+),typeof(-),typeof(*),typeof(/),typeof(^)},
            ::Type{<:$type},
            ::Type{NLPExpr}
            )
            return NLPExpr
        end
        function _MA.promote_operation(
            ::Union{typeof(+),typeof(-),typeof(*),typeof(/),typeof(^)},
            ::Type{NLPExpr},
            ::Type{<:$type}
            )
            return NLPExpr
        end
    end
end
function _MA.promote_operation(
    ::Union{typeof(+),typeof(-),typeof(*),typeof(/),typeof(^)},
    ::Type{NLPExpr},
    ::Type{NLPExpr}
    )
    return NLPExpr
end
for type in (:GeneralVariableRef, 
             :(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @eval begin
        function _MA.promote_operation(
            ::Union{typeof(*),typeof(/),typeof(^)},
            ::Type{<:$type},
            ::Type{JuMP.GenericQuadExpr{Float64, GeneralVariableRef}}
            )
            return NLPExpr
        end
        function _MA.promote_operation(
            ::Union{typeof(*),typeof(/),typeof(^)},
            ::Type{JuMP.GenericQuadExpr{Float64, GeneralVariableRef}},
            ::Type{<:$type}
            )
            return NLPExpr
        end
    end
end
function _MA.promote_operation(
    ::Union{typeof(*),typeof(/),typeof(^)},
    ::Type{<:JuMP.GenericQuadExpr{Float64, GeneralVariableRef}},
    ::Type{<:JuMP.GenericQuadExpr{Float64, GeneralVariableRef}}
    )
    return NLPExpr
end
function _MA.promote_operation(
    ::Union{typeof(/),typeof(^)},
    ::Type{Union{<:Real, <:AbstractInfOptExpr}},
    ::Type{<:AbstractInfOptExpr}
    )
    return NLPExpr
end

# Extend MA.scaling in case an NLPExpr needs to be converted to a number
function _MA.scaling(nlp::NLPExpr)
    c = _node_value(expr.data)
    if !(c isa Real) 
        throw(InexactError("Cannot convert `$nlp` to `$Float64`."))
    end
    return _MA.scaling(c)
end

# Extend MA.mutable_Copy to avoid unnecessary copying
function _MA.mutable_copy(nlp::NLPExpr) 
    return nlp # we don't need to copy since we build from the leaves up
end

# Extend MA.mutable_operate! as required 
function _MA.mutable_operate!(
    op::Union{typeof(zero), typeof(one)}, 
    ::NLPExpr
    ) 
    return op(NLPExpr) # not actually mutable for safety and efficiency
end
function _MA.mutable_operate!(
    op::Union{typeof(+), typeof(-), typeof(*), typeof(/), typeof(/)}, 
    nlp::NLPExpr,
    v
    ) 
    return op(nlp, v)
end
function _MA.mutable_operate!(
    op::Union{typeof(+), typeof(-), typeof(*), typeof(/), typeof(/)}, 
    v,
    nlp::NLPExpr
    ) 
    return op(v, nlp)
end
function _MA.mutable_operate!(
    op::Union{typeof(+), typeof(-), typeof(*), typeof(/), typeof(/)}, 
    nlp1::NLPExpr,
    nlp2::NLPExpr
    ) 
    return op(nlp1, nlp2)
end
function _MA.mutable_operate!(op::_MA.AddSubMul, nlp::NLPExpr, args...) 
    return _MA.add_sub_op(op)(nlp, *(args...))
end

# TODO maybe extend _MA.add_mul/_MA_.sub_mul as well

################################################################################
#                               SUMS AND PRODUCTS
################################################################################
## Define helper functions for sum reductions 
# Container of NLPExprs
function _reduce_by_first(::typeof(sum), first_itr::NLPExpr, itr)
    root = _LCRST.Node(NodeData(:+))
    prevc = _LCRST.addchild(root, first_itr.tree_root)
    for ex in itr
        prevc = _LCRST.addchild(root, prevc, _process_child_input(ex))
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
    return isempty(itr) ? first_itr : first_itr + sum(identity, itr; kw...)
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
    return _reduce_by_first(sum, itr1, new_itr)
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
    prevc = _LCRST.addchild(root, first_itr.tree_root)
    for ex in itr
        prevc = _LCRST.addchild(root, prevc, _process_child_input(ex))
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
    return _reduce_by_first(prod, itr1, new_itr)
end

################################################################################
#                           MULTIPLICATION OPERATORS
################################################################################
# TODO more intelligently operate with constants

# QuadExpr * expr
function Base.:*(
    quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef},
    expr::AbstractInfOptExpr
    )
    return NLPExpr(_call_graph(:*, quad, expr))
end

# expr * QuadExpr
function Base.:*(
    expr::AbstractInfOptExpr,
    quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef}
    )
    return NLPExpr(_call_graph(:*, expr, quad))
end

# QuadExpr * QuadExpr
function Base.:*(
    quad1::JuMP.GenericQuadExpr{Float64, GeneralVariableRef},
    quad2::JuMP.GenericQuadExpr{Float64, GeneralVariableRef}
    )
    return NLPExpr(_call_graph(:*, quad1, quad2))
end

# NLPExpr * QuadExpr
function Base.:*(
    nlp::NLPExpr,
    quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef}
    )
    return NLPExpr(_call_graph(:*, nlp, quad))
end

# QuadExpr * NLPExpr
function Base.:*(
    quad::JuMP.GenericQuadExpr{Float64, GeneralVariableRef}, 
    nlp::NLPExpr
    )
    return NLPExpr(_call_graph(:*, quad, nlp))
end

# NLPExpr * expr/constant
function Base.:*(nlp::NLPExpr, expr::Union{AbstractInfOptExpr, Real})
    return NLPExpr(_call_graph(:*, nlp, expr))
end

# expr/constant * NLPExpr
function Base.:*(expr::Union{AbstractInfOptExpr, Real}, nlp::NLPExpr)
    return NLPExpr(_call_graph(:*, expr, nlp))
end

# NLPExpr * NLPExpr
function Base.:*(nlp1::NLPExpr, nlp2::NLPExpr)
    return NLPExpr(_call_graph(:*, nlp1, nlp2))
end

# expr * expr * expr ...
function Base.:*(
    expr1::AbstractInfOptExpr,
    expr2::AbstractInfOptExpr,
    expr3::AbstractInfOptExpr,
    exprs::Vararg{AbstractInfOptExpr}
    )
    return NLPExpr(_call_graph(:*, expr1, expr2, expr3, exprs...))
end

# *NLPExpr
function Base.:*(nlp::NLPExpr)
    return nlp
end

################################################################################
#                              DIVISION OPERATORS
################################################################################
# expr/constant / expr
function Base.:/(
    expr1::Union{AbstractInfOptExpr, Real}, 
    expr2::AbstractInfOptExpr
    )
    return NLPExpr(_call_graph(:/, expr1, expr2))
end

# NLPExpr / constant
function Base.:/(nlp::NLPExpr, c::Real)
    if iszero(c)
        error("Cannot divide by zero.")
    elseif isone(c)
        return nlp
    else
        return NLPExpr(_call_graph(:/, nlp, c))
    end
end

################################################################################
#                                POWER OPERATORS
################################################################################
# expr ^ Integer
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

# expr ^ Real
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

# NLPExpr ^ Integer
function Base.:^(expr::NLPExpr, c::Integer)
    if iszero(c)
        return 1.0
    elseif isone(c)
        return expr 
    else 
        return NLPExpr(_call_graph(:^, expr, c))
    end
end

# NLPExpr ^ Real
function Base.:^(expr::NLPExpr, c::Real)
    if iszero(c)
        return 1.0
    elseif isone(c)
        return expr 
    else 
        return NLPExpr(_call_graph(:^, expr, c))
    end
end

# expr/constant ^ expr
function Base.:^(
    expr1::Union{AbstractInfOptExpr, Real}, 
    expr2::AbstractInfOptExpr
    )
    return NLPExpr(_call_graph(:^, expr1, expr2))
end

################################################################################
#                             SUBTRACTION OPERATORS
################################################################################
# TODO more intelligently operate with constants

# NLPExpr - expr/constant
function Base.:-(nlp::NLPExpr, expr::Union{AbstractInfOptExpr, Real})
    return NLPExpr(_call_graph(:-, nlp, expr))
end

# expr/constant - NLPExpr
function Base.:-(expr::Union{AbstractInfOptExpr, Real}, nlp::NLPExpr)
    return NLPExpr(_call_graph(:-, expr, nlp))
end

# NLPExpr - NLPExpr
function Base.:-(nlp1::NLPExpr, nlp2::NLPExpr)
    return NLPExpr(_call_graph(:-, nlp1, nlp2))
end

# -NLPExpr
function Base.:-(nlp::NLPExpr)
    return NLPExpr(_call_graph(:-, nlp))
end

# Var - Var (to avoid using v == v)
function Base.:-(lhs::V, rhs::V) where {V<:GeneralVariableRef}
    if isequal(lhs, rhs)
        return zero(JuMP.GenericAffExpr{Float64,V})
    else
        return JuMP.GenericAffExpr(0.0, 
                DataStructures.OrderedDict(lhs => 1.0, rhs => -1.0))
    end
end

################################################################################
#                              ADDITION OPERATORS
################################################################################
# TODO more intelligently operate with constants

# NLPExpr + expr/constant
function Base.:+(nlp::NLPExpr, expr::Union{AbstractInfOptExpr, Real})
    return NLPExpr(_call_graph(:+, nlp, expr))
end

# expr/constant + NLPExpr
function Base.:+(expr::Union{AbstractInfOptExpr, Real}, nlp::NLPExpr)
    return NLPExpr(_call_graph(:+, expr, nlp))
end

# NLPExpr + NLPExpr
function Base.:+(nlp1::NLPExpr, nlp2::NLPExpr)
    return NLPExpr(_call_graph(:+, nlp1, nlp2))
end

# +NLPExpr
function Base.:+(nlp::NLPExpr)
    return NLPExpr(_call_graph(:+, nlp))
end

################################################################################
#                             NATIVE NLP FUNCTIONS
################################################################################
# Store all of the native registered functions
const _NativeNLPFunctions = Dict{Symbol, Function}()

# List of 1 argument base functions to register
const _Base1ArgFuncList = (
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
for (name, func) in _Base1ArgFuncList
    # add it to the main storage dict
    _NativeNLPFunctions[name] = func
    # make an expression constructor
    @eval begin 
        function Base.$name(v::AbstractInfOptExpr)
            return NLPExpr(_call_graph($(quot(name)), v))
        end
    end
end

# Setup the Base functions with 2 arguments
for (name, func) in (:min => min, :max => max)
    # add it to the main storage dict
    _NativeNLPFunctions[name] = func
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

# Setup the ifelse function
_NativeNLPFunctions[:ifelse] = Core.ifelse

"""
    InfiniteOpt.ifelse(cond::NLPExpr, v1::Union{AbstractInfOptExpr, Real}, 
                       v2::Union{AbstractInfOptExpr, Real})::NLPExpr

A symbolic version of `Core.ifelse` that can be used to establish symbolic 
expressions with logic conditions. Note that is must be written 
`InfiniteOpt.ifelse` since it conflicts with `Core.ifelse`.

**Example**
```julia 
julia> InfiniteOpt.ifelse(x >= y, 0, y^3)
ifelse(x >= y, 0, y^3)
```
"""
function ifelse(
    cond::NLPExpr,
    v1::Union{AbstractInfOptExpr, Real}, 
    v2::Union{AbstractInfOptExpr, Real}
    )::NLPExpr
    return NLPExpr(_call_graph(:ifelse, cond, v1, v2))
end

# Setup the Base comparison functions
for (name, func) in (:< => Base.:(<), :(==) => Base.:(==), :> => Base.:(>),
                     :<= => Base.:(<=), :>= => Base.:(<=))
    # add it to the main storage dict
    _NativeNLPFunctions[name] = func
    # make an expression constructor
    @eval begin 
        function Base.$name(v::AbstractInfOptExpr, c::Real)
            return NLPExpr(_call_graph($(quot(name)), v, c))
        end
        function Base.$name(c::Real, v::AbstractInfOptExpr)
            return NLPExpr(_call_graph($(quot(name)), c, v))
        end
        function Base.$name(v1::AbstractInfOptExpr, v2::AbstractInfOptExpr)
            return NLPExpr(_call_graph($(quot(name)), v1, v2))
        end
        if $(quot(name)) in (:<, :>)
            function Base.$name(v1::GeneralVariableRef, v2::GeneralVariableRef)
                if isequal(v1, v2)
                    return false
                else
                    return NLPExpr(_call_graph($(quot(name)), v1, v2))
                end
            end
        else
            function Base.$name(v1::GeneralVariableRef, v2::GeneralVariableRef)
                if isequal(v1, v2)
                    return true
                else
                    return NLPExpr(_call_graph($(quot(name)), v1, v2))
                end
            end
        end
    end
end

# Setup the Base logical functions (we cannot extend && and || directly)
_NativeNLPFunctions[:&&] = Base.:&
_NativeNLPFunctions[:||] = Base.:|

# Logical And
function Base.:&(v::AbstractInfOptExpr, c::Bool)
    if c
        return v
    else
        return false
    end
end
function Base.:&(c::Bool, v::AbstractInfOptExpr)
    if c
        return v
    else
        return false
    end
end
function Base.:&(v1::AbstractInfOptExpr, v2::AbstractInfOptExpr)
    return NLPExpr(_call_graph(:&&, v1, v2))
end

# Logical Or
function Base.:|(v::AbstractInfOptExpr, c::Bool)
    if c
        return true
    else
        return v
    end
end
function Base.:|(c::Bool, v::AbstractInfOptExpr)
    if c
        return true
    else
        return v
    end
end
function Base.:|(v1::AbstractInfOptExpr, v2::AbstractInfOptExpr)
    return NLPExpr(_call_graph(:||, v1, v2))
end

const _Special1ArgFuncList = (
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
for (name, func) in _Special1ArgFuncList
    # add it to the main storage dict
    _NativeNLPFunctions[name] = func
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
"""
    name_to_function(model::InfiniteModel, name::Symbol)::Function 

Return the registered function that corresponds to `n`. Returns `nothing` if 
`n` does not correspond to a registered function. This helps retrieve of the 
functions of function names stored in `NLPExpr`s.
"""
function name_to_function(model::InfiniteModel, name::Symbol) 
    # TODO update to look for functions registered to the model as well
    return get(_NativeNLPFunctions, name, nothing)
end

"""
    all_registered_functions(model::InfiniteModel)::Vector{Function}

Retrieve all the functions that are currently registered to 
"""
function all_registered_functions(model::InfiniteModel) 
    # TODO update to look for functions registered to the model as well
    return collect(values(_NativeNLPFunctions))
end

# TODO add function registration 

################################################################################
#                               LINEAR ALGEBRA
################################################################################
# Extend LinearAlgebra.dot for increased efficiency
LinearAlgebra.dot(lhs::AbstractInfOptExpr, rhs::AbstractInfOptExpr) = lhs * rhs
LinearAlgebra.dot(lhs::AbstractInfOptExpr, rhs::Real) = lhs * rhs
LinearAlgebra.dot(lhs::Real, rhs::AbstractInfOptExpr) = lhs * rhs

# Implement promote_rule to help build better containers
function Base.promote_rule(::Type{NLPExpr}, ::Type{<:Real})
    return NLPExpr
end
function Base.promote_rule(::Type{NLPExpr}, ::Type{GeneralVariableRef})
    return NLPExpr
end
function Base.promote_rule(::Type{NLPExpr}, ::Type{JuMP.GenericAffExpr})
    return NLPExpr
end
function Base.promote_rule(::Type{NLPExpr}, ::Type{JuMP.GenericQuadExpr})
    return NLPExpr
end

################################################################################
#                                  PRINTING
################################################################################
# Define better printing for NodeData
function Base.show(io::IO, data::NodeData)
    return print(io, string(_node_value(data)))
end

# Map operators to their respective precedence (largest is highest priority)
const _Precedence = (; :^ => 6, Symbol("+u") => 5, Symbol("-u") => 5, :* => 4, 
                     :/ => 4, :+ => 3, :- => 3, :(==) => 2, :<= => 2, :>= => 2, 
                     :> => 2, :< => 2, :&& => 1, :|| => 1)

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
        return 3
    elseif isone(first(itr)[1])
        # we have a single variable
        return 10 # will always have precedence
    elseif first(itr)[1] == -1
        # we have a single unary negative variable
        return 5
    else
        # we have a single variable multiplied by some coefficient
        return 4
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
        return 3
    else
        # we only have a single quadratic term
        return 4
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
    is_op = !simple && raw_value isa Symbol && haskey(_Precedence, raw_value)
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
        str = str[1:prevind(str, end, length(op_str))]
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
        return str[1:prevind(str, end, 2)] * ")"
    end
end

# Extend JuMP.function_string for nonlinear expressions
function JuMP.function_string(mode, nlp::NLPExpr)
    return _expr_string(nlp.tree_root)
end
