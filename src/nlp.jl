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
    NodeData

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
    elseif all(_is_zero(n) for n in node) && iszero(Base.get(_NativeNLPFunctions, (raw, length(collect(node))), (i...) -> true)((0.0 for n in node)...))
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
        elseif _is_zero(node.child.sibling)
            _LCRST.prunebranch!(node.child.sibling)
            _merge_parent_and_child(node)
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
    for (v, c) in aff.terms
        if isone(c)
            push!(ex.args, map_func(v))
        else
            push!(ex.args, Expr(:call, :*, c, map_func(v)))
        end
    end
    if !iszero(aff.constant)
        push!(ex.args, aff.constant)
    end
    return ex
end

# QuadExpr
function _ast_process_node(map_func::Function, quad::JuMP.GenericQuadExpr)
    ex = Expr(:call, :+)
    for (xy, c) in quad.terms
        if isone(c)
            push!(ex.args, Expr(:call, :*, map_func(xy.a), map_func(xy.b)))
        else
            push!(ex.args, Expr(:call, :*, c, map_func(xy.a), map_func(xy.b)))
        end
    end
    append!(ex.args, _ast_process_node(map_func, quad.aff).args[2:end])
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
#                               SUMS AND PRODUCTS
################################################################################
## Define helper functions for sum reductions 
# Container of NLPExprs
function _reduce_by_first(::typeof(sum), first_itr::NLPExpr, itr, orig_itr; kws...)
    for kw in kws
        error("Unexpected keyword argument `$kw`.")
    end
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
    itr,
    orig_itr;
    kws...
    )
    for kw in kws
        error("Unexpected keyword argument `$kw`.")
    end
    result = first_itr
    for i in itr 
        result = _MA.operate!!(_MA.add_mul, result, i)
    end
    return result
end

# Fallback
function _reduce_by_first(::typeof(sum), first_itr, itr, orig_itr; kws...)
    return sum(identity, orig_itr; kws...)
end

# Hyjack Base.sum for better efficiency on iterators --> this is type piracy...
function Base.sum(itr::Base.Generator; kws...)
    isempty(itr) && return sum(identity, itr; kws...)
    itr1, new_itr = Iterators.peel(itr)
    return _reduce_by_first(sum, itr1, new_itr, itr; kws...)
end

# Extend Base.sum for container of NLPExprs
function Base.sum(arr::AbstractArray{<:NLPExpr}; init = zero(NLPExpr))
    isempty(arr) && return init
    itr1, new_itr = Iterators.peel(arr)
    return _reduce_by_first(sum, itr1, new_itr, arr)
end

# Extend Base.sum for container of InfiniteOpt exprs
function Base.sum(
    arr::AbstractArray{<:AbstractInfOptExpr}; 
    init = zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
    )
    isempty(arr) && return init
    result = _MA.Zero()
    for i in arr 
        result = _MA.operate!!(_MA.add_mul, result, i)
    end
    return result
end

## Define helper functions for reducing products 
# Container of InfiniteOpt exprs
function _reduce_by_first(::typeof(prod), first_itr::AbstractInfOptExpr, itr, orig_itr; kws...)
    for kw in kws
        error("Unexpected keyword argument `$kw`.")
    end
    root = _LCRST.Node(NodeData(:*))
    prevc = _LCRST.addchild(root, _process_child_input(first_itr))
    for ex in itr
        prevc = _LCRST.addchild(root, prevc, _process_child_input(ex))
    end
    return NLPExpr(root)
end

# Fallback
function _reduce_by_first(::typeof(prod), first_itr, itr, orig_itr; kws...)
    return prod(identity, orig_itr; kws...)
end

# Hyjack Base.prod for better efficiency on iterators --> this is type piracy...
function Base.prod(itr::Base.Generator; kws...)
    isempty(itr) && return prod(identity, itr; kws...)
    itr1, new_itr = Iterators.peel(itr)
    return _reduce_by_first(prod, itr1, new_itr, itr; kws...)
end

# Extend Base.prod for container of InfiniteOpt exprs
function Base.prod(arr::AbstractArray{<:AbstractInfOptExpr}; init = one(NLPExpr))
    isempty(arr) && return init
    itr1, new_itr = Iterators.peel(arr)
    return _reduce_by_first(prod, itr1, new_itr, arr)
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
        return one(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
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
        return one(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
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
        return one(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
    elseif isone(c)
        return expr 
    else 
        return NLPExpr(_call_graph(:^, expr, c))
    end
end

# NLPExpr ^ Real
function Base.:^(expr::NLPExpr, c::Real)
    if iszero(c)
        return one(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
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
    return nlp
end

################################################################################
#                               MUTABLE ARITHMETICS
################################################################################
# Define NLPExpr as a mutable type for MA
_MA.mutability(::Type{NLPExpr}) = _MA.IsMutable()

# Extend MA.promote_operation for bettered efficiency
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
for type in (:GeneralVariableRef,
             :(JuMP.GenericAffExpr{Float64, GeneralVariableRef}), 
             :(JuMP.GenericQuadExpr{Float64, GeneralVariableRef}))
    @eval begin
        function _MA.promote_operation(
            ::Union{typeof(/),typeof(^)},
            ::Type{<:Real},
            ::Type{<:$type}
            )
            return NLPExpr
        end
    end
end
for type in (:GeneralVariableRef,
             :(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @eval begin
        function _MA.promote_operation(
            ::Union{typeof(/),typeof(^)},
            ::Type{GeneralVariableRef},
            ::Type{<:$type}
            )
            return NLPExpr
        end
    end
end
for type in (:GeneralVariableRef,
             :(JuMP.GenericAffExpr{Float64, GeneralVariableRef}))
    @eval begin
        function _MA.promote_operation(
            ::Union{typeof(/),typeof(^)},
            ::Type{JuMP.GenericAffExpr{Float64, GeneralVariableRef}},
            ::Type{<:$type}
            )
            return NLPExpr
        end
    end
end

# Extend MA.scaling in case an NLPExpr needs to be converted to a number
function _MA.scaling(nlp::NLPExpr)
    c = _node_value(nlp.tree_root.data)
    if !(c isa Real) 
        error("Cannot convert `$nlp` to `$Float64`.")
    end
    return _MA.scaling(c)
end

# Extend MA.mutable_copy to avoid unnecessary copying
function _MA.mutable_copy(nlp::NLPExpr) 
    return nlp # we don't need to copy since we build from the leaves up
end

# Extend MA.operate! as required 
function _MA.operate!(
    op::Union{typeof(zero), typeof(one)}, 
    ::NLPExpr
    ) 
    return op(NLPExpr) # not actually mutable for safety and efficiency
end
function _MA.operate!(
    op::Union{typeof(+), typeof(-), typeof(*), typeof(/), typeof(^)}, 
    nlp::NLPExpr,
    v
    ) 
    return op(nlp, v)
end
function _MA.operate!(
    op::Union{typeof(+), typeof(-), typeof(*), typeof(/), typeof(^)}, 
    v,
    nlp::NLPExpr
    ) 
    return op(v, nlp)
end
function _MA.operate!(
    op::typeof(+), 
    v::Union{JuMP.GenericAffExpr{Float64, GeneralVariableRef}, 
             JuMP.GenericQuadExpr{Float64, GeneralVariableRef}},
    nlp::NLPExpr
    ) 
    return op(v, nlp)
end
function _MA.operate!(
    op::typeof(-), 
    v::Union{JuMP.GenericAffExpr{Float64, GeneralVariableRef}, 
             JuMP.GenericQuadExpr{Float64, GeneralVariableRef}},
    nlp::NLPExpr
    ) 
    return op(v, nlp)
end
function _MA.operate!(
    op::Union{typeof(+), typeof(-), typeof(*), typeof(/), typeof(^)}, 
    nlp1::NLPExpr,
    nlp2::NLPExpr
    ) 
    return op(nlp1, nlp2)
end
function _MA.operate!(op::_MA.AddSubMul, nlp::NLPExpr, args...) 
    return _MA.add_sub_op(op)(nlp, *(args...))
end

# TODO maybe extend _MA.add_mul/_MA_.sub_mul as well

################################################################################
#                             NATIVE NLP FUNCTIONS
################################################################################
# Store all of the native registered functions
const _NativeNLPFunctions = Dict{Tuple{Symbol, Int}, Function}(
    (:-, 2) => -, 
    (:/, 2) => /, 
    (:^, 2) => ^
)

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
    _NativeNLPFunctions[(name, 1)] = func
    # make an expression constructor
    @eval begin 
        function Base.$name(v::AbstractInfOptExpr)
            return NLPExpr(_call_graph($(quot(name)), v))
        end
    end
end

# Setup the ifelse function
_NativeNLPFunctions[(:ifelse, 3)] = Core.ifelse

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
    )
    return NLPExpr(_call_graph(:ifelse, cond, v1, v2))
end
function ifelse(
    cond::Bool,
    v1::Union{AbstractInfOptExpr, Real}, 
    v2::Union{AbstractInfOptExpr, Real}
    )
    return cond ? v1 : v2
end

# Setup the Base comparison functions
for (name, func) in (:< => Base.:(<), :(==) => Base.:(==), :> => Base.:(>),
                     :<= => Base.:(<=), :>= => Base.:(>=))
    # add it to the main storage dict
    _NativeNLPFunctions[(name, 2)] = func
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
_NativeNLPFunctions[(:&&, 2)] = Base.:&
_NativeNLPFunctions[(:||, 2)] = Base.:|

# Logical And
function Base.:&(v::Union{GeneralVariableRef, NLPExpr}, c::Bool)
    return c ? v : false
end
function Base.:&(c::Bool, v::Union{GeneralVariableRef, NLPExpr})
    return c ? v : false
end
function Base.:&(
    v1::Union{GeneralVariableRef, NLPExpr}, 
    v2::Union{GeneralVariableRef, NLPExpr}
    )
    return NLPExpr(_call_graph(:&&, v1, v2))
end

# Logical Or
function Base.:|(v::Union{GeneralVariableRef, NLPExpr}, c::Bool)
    return c ? true : v
end
function Base.:|(c::Bool, v::Union{GeneralVariableRef, NLPExpr})
    return c ? true : v
end
function Base.:|(
    v1::Union{GeneralVariableRef, NLPExpr}, 
    v2::Union{GeneralVariableRef, NLPExpr})
    return NLPExpr(_call_graph(:||, v1, v2)
    )
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
    _NativeNLPFunctions[(name, 1)] = func
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
    RegisteredFunction{F <: Function, G <: Union{Function, Nothing}, 
                       H <: Union{Function, Nothing}}

A type for storing used defined registered functions and their information that 
is needed by JuMP for build an `NLPEvaluator`. The constructor is of the form:
```julia
    RegisteredFunction(name::Symbol, num_args::Int, func::Function, 
                       [gradient::Function, hessian::Function])
```

**Fields**
- `name::Symbol`: The name of the function that is used in `NLPExpr`s.
- `num_args::Int`: The number of function arguments.
- `func::F`: The function itself.
- `gradient::G`: The gradient function if one is given.
- `hessian::H`: The hessian function if one is given.
"""
struct RegisteredFunction{F <: Function, G, H}
    name::Symbol
    num_args::Int
    func::F
    gradient::G
    hessian::H

    # Constructors
    function RegisteredFunction(
        name::Symbol, 
        num_args::Int, 
        func::F
        ) where {F <: Function}
        return new{F, Nothing, Nothing}(name, num_args, func, nothing, nothing)
    end
    function RegisteredFunction(
        name::Symbol, 
        num_args::Int, 
        func::F,
        gradient::G
        ) where {F <: Function, G <: Function}
        if isone(num_args) && !hasmethod(gradient, Tuple{Real})
            error("Invalid gradient function form, see the docs for details.")
        elseif !isone(num_args) && !hasmethod(gradient, Tuple{AbstractVector{Real}, ntuple(_->Real, num_args)...})
            error("Invalid multi-variate gradient function form, see the docs for details.")
        end
        return new{F, G, Nothing}(name, num_args, func, gradient, nothing)
    end
    function RegisteredFunction(
        name::Symbol, 
        num_args::Int, 
        func::F,
        gradient::G,
        hessian::H
        ) where {F <: Function, G <: Function, H <: Function}
        if isone(num_args) && !hasmethod(gradient, Tuple{Real})
            error("Invalid gradient function form, see the docs for details.")
        elseif isone(num_args) && !hasmethod(hessian, Tuple{Real})
            error("Invalid hessian function form, see the docs for details.")
        end 
        return new{F, G, H}(name, num_args, func, gradient, hessian)
    end
end

# Helper function for @register
function _register(
    _error::Function, 
    call_mod::Module,
    model::InfiniteModel, 
    name::Symbol, 
    num_args::Int, 
    funcs...
    )
    if !all(f -> f isa Function, funcs) 
        _error("Gradient and/or hessian must be functions.")
    elseif haskey(_NativeNLPFunctions, (name, num_args)) || 
           haskey(model.func_lookup, (name, num_args))
        _error("A function with name `$name` and $num_args arguments is already " *
               "registered. Please use a function with a different name.")
    elseif !hasmethod(funcs[1], NTuple{num_args, Real})
        _error("The function `$name` is not defined for arguments of type `Real`.")
    elseif length(unique!([m.module for m in methods(funcs[1])])) > 1 || 
           first(methods(funcs[1])).module !== call_mod
        _error("Cannot register function names that are used by packages. Try " * 
               "wrapping `$(funcs[1])` in a user-defined function.")
    end
    push!(model.registrations, RegisteredFunction(name, num_args, funcs...))
    model.func_lookup[name, num_args] = funcs[1]
    return
end

# Helper function to check the inputs of created functions 
function _check_function_args(model::InfiniteModel, f_name, args...)
    for a in args 
        m = _model_from_expr(a)
        if !isnothing(m) && m !== model
            error("`$f_name` is a registered function in a different model than " *
                  "`$a` belongs to. Try registering `$f_name` to the current " * 
                  "model.")
        end
    end
    return
end

"""
    @register(model::InfiniteModel, func_expr, [gradient::Function], [hessian::Function])

Register a user-defined function in accordance with `func_expr` such that it can 
be used in `NLPExpr`s that are used with `model` without being traced.

**Argument Information**
Here `func_expr` is of the form: `myfunc(a, b)` where `myfunc` is the function 
name and the number of arguments are given symbolically. Note that the choice 
of argument symbols is arbitrary. Each function argument must support anything 
of type `Real` to specified.

Here we can also specify a gradient function `gradient` which for 1 argument 
functions must taken in the 1 argument and return its derivative. For 
multi-argument functions the gradient function must be of the form:
```julia
function gradient(g::AbstractVector{T}, args::T...) where {T <: Real}
    # fill g vector with the gradient of the function
end
```

For 1 argument functions we can also specify a hessian function with takes that 
argument and return the 2nd derivative. Hessians can ge specified for 
multi-argument functions, but `JuMP` backends do not currently support this.

If no gradient and/or hessian is given, the automatic differentation capabilities 
of the backend (e.g., `JuMP`) will be used to determine them. Note that the 
`JuMP` backend does not use Hessian's for user-defined multi-argument functions.

**Notes**
- When possible, tracing is preferred over registering a function (see 
  [Function Tracing](@ref) for more info).
- Only user-defined functions can be specified. If the function is used by a 
  package then it can not be used directly. However, we can readily wrap it in a 
  new function `newfunc(a) = pkgfunc(a)`.
- We can only register functions in the same scope that they are defined in.
- Registered functions can only be used in or below the scope in which they are 
  registered. For instance, if we register some function inside of another 
  function then we can only use it inside that function (not outside of it). 
- A function with a given name and number of arguments can only be registered 
  once in a particular model.

**Examples**
```julia-repl
julia> @variable(model, x)
x

julia> f(a) = a^3;

julia> f(x) # user-function gets traced
x^3

julia> @register(model, f(a)) # register function
f (generic function with 2 methods)

julia> f(x) # function is no longer traced and autodifferentiation will be used
f(x)

julia> f2(a) = a^2; g2(a) = 2 * a; h2(a) = 2;

julia> @register(model, f2(a), g2, h2) # register with explicit gradient and hessian
f2 (generic function with 2 methods)

julia> f2(x)
f2(x)

julia> f3(a, b) = a * b^2;

julia> function g3(v, a, b)
          v[1] = b^2
          v[2] = 2 * a * b
          return
       end;

julia> @register(model, f3(a, b), g3) # register multi-argument function
f3 (generic function with 4 methods)

julia> f3(42, x)
f3(42, x)
```
"""
macro register(model, f, args...)
     # define error message function
     _error(str...) = _macro_error(:register, (f, args...), __source__, str...)

    # parse the arguments and check
    pos_args, extra_kwargs, _, _ = _extract_kwargs(args)
    if !isempty(extra_kwargs)
        _error("Keyword arguments were given, but none are accepted.")
    elseif length(pos_args) > 2
        _error("Too many position arguments given, should be of form " * 
               "`@register(myfunc(a), [gradient], [hessian])` where " * 
               "`gradient` and `hessian` are optional arguments.") 
    end

    # process the function input
    if isexpr(f, :call) && all(a -> a isa Symbol, f.args)
        f_name = f.args[1]
        f_args = f.args[2:end]
        num_args = length(f_args)
    else
        _error("Unexpected function format, should be of form `myfunc(a, b)`.")
    end

    # start creating the register code and register 
    code = Expr(:block)
    push!(code.args, quote 
        $model isa InfiniteModel || $_error("Expected an `InfiniteModel`.")
    end)
    calling_mod = __module__ # get the module the macro is being called from
    push!(code.args, quote 
        InfiniteOpt._register($_error, $calling_mod, $model, $(quot(f_name)), 
                              $num_args, $(f_name), $(args...))
    end)

    # define the function overloads needed to create expressions
    Ts = [Real, AbstractInfOptExpr]
    type_combos = vec(collect(Iterators.product(ntuple(_->Ts, num_args)...)))
    filter!(ts -> !all(T -> T == Real, ts), type_combos) # remove combo with only Reals
    annotype(name, T) = :($name :: $T)
    set_args(xs, vs) = (xs = map(annotype, xs, vs); xs)
    for ts in type_combos
        push!(code.args, quote
            function $(f_name)($(set_args(f_args, ts)...))
                InfiniteOpt._check_function_args($model, $(quot(f_name)), $(f_args...))
                return NLPExpr(InfiniteOpt._call_graph($(quot(f_name)), $(f_args...)))
            end
        end)
    end

    # return the code 
    return esc(code)
end

"""
    name_to_function(model::InfiniteModel, name::Symbol, num_args::Int)::Union{Function, Nothing} 

Return the registered function that corresponds to `name` with `num_args`. 
Returns `nothing` if no such registered function exists. This helps retrieve the 
functions of function names stored in `NLPExpr`s.
"""
function name_to_function(model::InfiniteModel, name::Symbol, num_args::Int)
    if name == :+
        return +
    elseif name == :*
        return *
    else
        return Base.get(_NativeNLPFunctions, (name, num_args), 
                   Base.get(model.func_lookup, (name, num_args), nothing))
    end
end

"""
    all_registered_functions(model::InfiniteModel)::Vector{Function}

Retrieve all the functions that are currently registered to `model`.
"""
function all_registered_functions(model::InfiniteModel) 
    funcs = append!(collect(values(_NativeNLPFunctions)), (+, *))
    return append!(funcs, values(model.func_lookup))
end

"""
    user_registered_functions(model::InfiniteModel)::Vector{RegisteredFunction}

Return all the functions (and their associated information) that the user has 
registered to `model`. Each is stored as a [`RegisteredFunction`](@ref).
"""
function user_registered_functions(model::InfiniteModel)
    return model.registrations
end

## Define helper function to add registered functions to JuMP
# No gradient or hessian
function _add_func_data_to_jump(
    model::JuMP.Model, 
    data::RegisteredFunction{F, Nothing, Nothing}
    ) where {F <: Function}
    JuMP.register(model, data.name, data.num_args, data.func, autodiff = true)
    return
end

# Only gradient information
function _add_func_data_to_jump(
    model::JuMP.Model, 
    data::RegisteredFunction{F, G, Nothing}
    ) where {F <: Function, G <: Function}
    JuMP.register(model, data.name, data.num_args, data.func, data.gradient, 
                  autodiff = isone(data.num_args))
    return
end

# Gradient and hessian information
function _add_func_data_to_jump(model::JuMP.Model, data::RegisteredFunction)
    if data.num_args > 1
        error("JuMP does not support hessians for multi-argument registered " * 
              "functions.")
    end
    JuMP.register(model, data.name, data.num_args, data.func, data.gradient, 
                  data.hessian)
    return
end

"""
    add_registered_to_jump(opt_model::JuMP.Model, inf_model::InfiniteModel)::Nothing

Add the user registered functions in `inf_model` to a `JuMP` model `opt_model`. 
This is intended as an internal method, but it is provided for developers that 
extend `InfiniteOpt` to use other optimizer models.
"""
function add_registered_to_jump(opt_model::JuMP.Model, inf_model::InfiniteModel)
    for data in user_registered_functions(inf_model)
        _add_func_data_to_jump(opt_model, data)
    end
    return
end

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
function Base.promote_rule(::Type{NLPExpr}, ::Type{<:JuMP.GenericAffExpr})
    return NLPExpr
end
function Base.promote_rule(::Type{NLPExpr}, ::Type{<:JuMP.GenericQuadExpr})
    return NLPExpr
end

# TODO make proper MA extensions to enable efficient definition

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
    prev_prec = 0,
    prev_comm = false
    )
    # prepocess the raw value
    raw_value = _node_value(node.data)
    is_op = raw_value isa Symbol && haskey(_Precedence, raw_value)
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
            str = _expr_string(child, str)
            str *= ", "
        end
        return str[1:prevind(str, end, 2)] * ")"
    end
end

# Extend JuMP.function_string for nonlinear expressions
function JuMP.function_string(mode, nlp::NLPExpr)
    return _expr_string(nlp.tree_root)
end
