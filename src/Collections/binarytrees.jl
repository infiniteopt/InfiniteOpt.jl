mutable struct BinaryNode{T}
    data::T
    parent::BinaryNode{T}
    left::BinaryNode{T}
    right::BinaryNode{T}

    # Root constructor
    BinaryNode{T}(data) where T = new{T}(data)
    # Child node constructor
    BinaryNode{T}(data, parent::BinaryNode{T}) where T = new{T}(data, parent)
end
BinaryNode(data) = BinaryNode{typeof(data)}(data)

function leftchild(data::T, parent::BinaryNode{T}) where {T}
    !isdefined(parent, :left) || error("left child is already assigned")
    node = BinaryNode{T}(data, parent)
    parent.left = node
end
function rightchild(data::T, parent::BinaryNode{T}) where {T}
    !isdefined(parent, :right) || error("right child is already assigned")
    node = BinaryNode{T}(data, parent)
    parent.right = node
end

function leftchild(child::BinaryNode{T}, parent::BinaryNode{T}) where {T}
    !isdefined(parent, :left) || error("left child is already assigned")
    !isdefined(child, :parent) || error("child already has parent")
    child.parent = parent
    parent.left = child
end
function rightchild(child::BinaryNode{T}, parent::BinaryNode{T}) where {T}
    !isdefined(parent, :right) || error("right child is already assigned")
    !isdefined(child, :parent) || error("child already has parent")
    child.parent = parent
    parent.right = child
end

# Implement iteration over the immediate children of a node
function Base.iterate(node::BinaryNode)
    isdefined(node, :left) && return (node.left, false)
    isdefined(node, :right) && return (node.right, true)
    return nothing
end
function Base.iterate(node::BinaryNode, state::Bool)
    state && return nothing
    isdefined(node, :right) && return (node.right, true)
    return nothing
end
Base.IteratorSize(::Type{BinaryNode{T}}) where T = Base.SizeUnknown()
Base.eltype(::Type{BinaryNode{T}}) where T = BinaryNode{T}

## Things we need to define to leverage the native iterator over children
## for the purposes of AbstractTrees.
# Set the traits of this kind of tree
Base.eltype(::Type{<:AbstractTrees.TreeIterator{BinaryNode{T}}}) where T = BinaryNode{T}
Base.IteratorEltype(::Type{<:AbstractTrees.TreeIterator{BinaryNode{T}}}) where T = Base.HasEltype()
AbstractTrees.parentlinks(::Type{BinaryNode{T}}) where T = AbstractTrees.StoredParents()
AbstractTrees.siblinglinks(::Type{BinaryNode{T}}) where T = AbstractTrees.StoredSiblings()
# Use the native iteration for the children
AbstractTrees.children(node::BinaryNode) = node
has_children(node::BinaryNode) = isdefined(node, :left) || isdefined(node, :right)

# Base.parent(root::BinaryNode, node::BinaryNode) = isdefined(node, :parent) ? node.parent : nothing

function AbstractTrees.nextsibling(tree::BinaryNode, child::BinaryNode)
    isdefined(child, :parent) || return nothing
    p = child.parent
    if isdefined(p, :right)
        child === p.right && return nothing
        return p.right
    end
    return nothing
end

# We also need `pairs` to return something sensible.
# If you don't like integer keys, you could do, e.g.,
#   Base.pairs(node::BinaryNode) = BinaryNodePairs(node)
# and have its iteration return, e.g., `:left=>node.left` and `:right=>node.right` when defined.
# But the following is easy:
Base.pairs(node::BinaryNode) = enumerate(node)

AbstractTrees.printnode(io::IO, node::BinaryNode) = print(io, node.data)

function Base.copy(node::BinaryNode)
    new_node = BinaryNode(node.data)
    if isdefined(node, :left)
        leftchild(copy(node.left), new_node)
    end
    if isdefined(node, :right)
        rightchild(copy(node.right), new_node)
    end
    return new_node
end
