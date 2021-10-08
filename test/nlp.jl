# Test the extensions to LCRST
@testset "LCRST Extensions" begin 
    # test additions to addchild
    @testset "addchild(parent, child)" begin 
        # setup nodes
        p = Node(1)
        c1 = Node(2)
        c2 = Node(3)
        # test first addition of first child
        @test addchild(p, c1) == c1
        @test c1.parent == p
        @test lastsibling(c1) == c1
        @test c1.child == c1
        @test p.child == c1
        # test first addition of 2nd child
        @test addchild(p, c2) == c2
        @test c1.parent == p
        @test c2.parent == p
        @test c1.sibling == c2
        @test c2.sibling == c2
        @test c2.child == c2
        @test c1.child == c1
        @test p.child == c1
        # test subsequent addition of first child (invokes copy)
        @test addchild(p, c1) !== c1
        @test c1.parent == p
        @test lastsibling(c1) !== c2 && lastsibling(c1) !== c1
        @test c1.child == c1
        @test p.child == c1
        c3 = lastsibling(c1)
        @test c3.parent == p
        @test c2.sibling == c3 
        @test c3.child == c3 
        @test c3 !== c1
        @test c3.data == 2
        # test subsequent addition of child to a node with no children 
        p2 = Node(4)
        @test addchild(p2, c1) !== c1
        c4 = p2.child
        @test c4.parent == p2
        @test c4.data == 2 
        @test c4.sibling == c4 
        @test c4.child == c4
    end
    @testset "addchild(parent, nothing, child)" begin 
        # setup nodes
        p = Node(1)
        c = Node(2)
        # test first addition of first child
        @test addchild(p, nothing, c) == c
        @test c.parent == p
        @test lastsibling(c) == c
        @test c.child == c
        @test p.child == c
        # test subsequent addition of child to a node with no children 
        p2 = Node(4)
        @test addchild(p2, nothing, c) !== c
        c2 = p2.child
        @test c2.parent == p2
        @test c2.data == 2 
        @test c2.sibling == c2
        @test c2.child == c2
    end
    @testset "addchild(parent, prevchild, data)" begin 
        # setup nodes
        p = Node(1)
        c1 = Node(2)
        p.child = c1
        c1.parent = p
        # test addition on first next sibling
        @test addchild(p, c1, 3) isa Node 
        c2 = c1.sibling
        @test c1 !== c2
        @test c2.data == 3
        @test c2.parent == p
        @test c2.sibling == c2
        @test c2.child == c2
        @test c1.sibling == c2
        @test c1.parent == p
        @test c1.child == c1
        # test adding more in for loop
        prev = c2
        for i in 4:5
            prev = addchild(p, prev, i)
        end
        @test [n.data for n in p] == [2, 3, 4, 5]
        c3 = c2.sibling 
        c4 = c3.sibling 
        @test c2 !== c3 && c3 !== c4
        @test c3.parent == p 
        @test c3.sibling == c4
        @test c3.child == c3
        @test c4.parent == p
        @test c4.sibling == c4
    end
    @testset "addchild(parent, prevchild, child)" begin 
        # setup nodes
        p = Node(1)
        c1 = Node(2)
        p.child = c1
        c1.parent = p
        c2 = Node(3)
        # test invalid previous 
        @test_throws AssertionError addchild(p, c2, c1)
        # test first addition of 2nd child
        @test addchild(p, c1, c2) == c2
        @test c1.parent == p
        @test c2.parent == p
        @test c1.sibling == c2
        @test c2.sibling == c2
        @test c2.child == c2
        @test c1.child == c1
        @test p.child == c1
        # test subsequent addition of first child (invokes copy)
        @test addchild(p, c2, c1) !== c1
        @test c1.parent == p
        @test lastsibling(c1) !== c2 && lastsibling(c1) !== c1
        @test c1.child == c1
        @test p.child == c1
        c3 = lastsibling(c1)
        @test c3.parent == p
        @test c2.sibling == c3 
        @test c3.child == c3 
        @test c3 !== c1
        @test c3.data == 2
        # test subsequent double addition of a child to a root
        p = Node(1)
        c1 = Node(2)
        p.child = c1
        c1.parent = p
        @test addchild(p, c1, c1) !== c1
        @test [n.data for n in p] == [2, 2]
        c2 = c1.sibling
        @test c1 !== c2
        @test c1.parent == p
        @test c2.parent == p
        @test c1.sibling == c2
        @test c2.sibling == c2
        @test c2.child == c2
        @test c1.child == c1
        @test p.child == c1
        # test for loop addition of children 
        p = Node(1)
        cs = [Node(2), Node(3), Node(4)]
        prev = nothing 
        for i in 1:3
            prev = addchild(p, prev, cs[i])
        end
        @test p.child == cs[1]
        for i in 1:3
            @test cs[i].data == i + 1
            @test cs[i].parent == p
            @test lastsibling(cs[i]) == cs[3]
        end
    end
    # test _map_tree
    @testset "_map_tree" begin 
        # setup tree
        p = Node(1)
        c1 = addchild(p, 2)
        c2 = addchild(p, 3)
        c3 = addchild(c2, 4)
        # test simple map 
        str_p = InfiniteOpt._map_tree(n -> Node(string(n.data)), p)
        @test isroot(str_p)
        @test str_p.data == "1"
        @test str_p.child.data == "2"
        @test isleaf(str_p.child)
        @test str_p.sibling == str_p
        @test str_p.child.sibling.data == "3"
        @test islastsibling(str_p.child.sibling)
        @test isleaf(str_p.child.sibling.child)
        @test str_p.child.sibling.child.data == "4"
        # test empty root 
        p = Node(1)
        str_p = InfiniteOpt._map_tree(n -> Node(string(n.data)), p)
        @test isroot(str_p)
        @test isleaf(str_p)
        @test str_p.data == "1"
    end
    # test copy(Node)
    @testset "Base.copy" begin 
        # setup tree
        p = Node(1)
        c1 = addchild(p, 2)
        c2 = addchild(p, 3)
        c3 = addchild(c2, 4)
        # test simple copy
        p2 = copy(p)
        @test p2 !== p
        @test isroot(p2)
        @test p2.data == 1
        @test p2.child.data == 2
        @test isleaf(p2.child)
        @test p2.child !== p.child
        @test p2.sibling == p2
        @test p2.child.sibling.data == 3
        @test islastsibling(p2.child.sibling)
        @test isleaf(p2.child.sibling.child)
        @test p2.child.sibling.child.data == 4
        # test empty root 
        p = Node(1)
        p2 = copy(p)
        @test p !== p2
        @test isroot(p2)
        @test isleaf(p2)
        @test p2.data == 1
    end
    # test _merge_parent_and_child
    @testset "_merge_parent_and_child" begin 
        # setup tree
        p = Node(1)
        c1 = addchild(p, 2)
        c2 = addchild(c1, 3)
        c3 = addchild(c2, 4)
        c4 = addchild(c2, 5)
        # test normal 
        @test InfiniteOpt._merge_parent_and_child(c1) == c1 
        @test c1.data == 3
        @test c1.child == c3
        @test c3.parent == c1 
        @test c4.parent == c1
        @test c1.parent == p
        # test at root node 
        @test InfiniteOpt._merge_parent_and_child(p) == p
        @test p.data == 3
        @test p.child == c3
        @test c3.parent == p
        @test c4.parent == p
        @test p.parent == p
        # test nothing happens 
        @test InfiniteOpt._merge_parent_and_child(p) == p
        @test p.data == 3
        @test p.child == c3
        @test c3.parent == p
        @test c4.parent == p
        @test p.parent == p
    end
end

# Test the core NLP data structures and methods 
@testset "Core Data Structures" begin 
    # setup model data 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @variable(m, y, Infinite(t))
    @variable(m, z)
    # test NodeData
    @testset "NodeData" begin 
        @test NodeData(1).value == 1
    end
    # test _node_value
    @testset "_node_value" begin 
        @test InfiniteOpt._node_value(NodeData(1)) == 1
    end
    # test _is_zero
    @testset "_is_zero" begin 
        # test easy case 
        @test InfiniteOpt._is_zero(Node(NodeData(0)))
        # test is leaf that isn't zero
        @test !InfiniteOpt._is_zero(Node(NodeData(1)))
        @test !InfiniteOpt._is_zero(Node(NodeData(y)))
        # test addition block 
        p = Node(NodeData(:+)) 
        addchild(p, NodeData(y)) 
        addchild(p, NodeData(0))
        @test !InfiniteOpt._is_zero(p)
        p.child.data = NodeData(0)
        @test InfiniteOpt._is_zero(p)
        # test multiplication block
        p = Node(NodeData(:*)) 
        addchild(p, NodeData(y)) 
        addchild(p, NodeData(0))
        @test InfiniteOpt._is_zero(p)
        p.child.sibling.data = NodeData(t)
        @test !InfiniteOpt._is_zero(p)
        # test power/divide
        p = Node(NodeData(:^))
        addchild(p, NodeData(y)) 
        addchild(p, NodeData(2))
        @test !InfiniteOpt._is_zero(p)
        p.child.data = NodeData(0)
        @test InfiniteOpt._is_zero(p)
        # test function 
        p = Node(NodeData(:max))
        addchild(p, NodeData(0))
        addchild(p, NodeData(0)) 
        @test InfiniteOpt._is_zero(p)
        p.child.data = NodeData(y)
        @test !InfiniteOpt._is_zero(p)
    end
    # test _drop_zeros!
    @testset "_drop_zeros!" begin 
        # test simple case 
        n = Node(NodeData(0))
        @test InfiniteOpt._drop_zeros!(n) == n
        @test n.data.value == 0
        # test addition with zeros
        p = Node(NodeData(:+))
        addchild(p, NodeData(0))
        addchild(p, NodeData(y))
        @test InfiniteOpt._drop_zeros!(p) == p
        @test isequal(p.data.value, y)
        # test case with subtraction and zero function 
        p = Node(NodeData(:-))
        n1 = addchild(p, NodeData(:sin))
        n2 = addchild(p, NodeData(y))
        n3 = addchild(n1, NodeData(0))
        @test InfiniteOpt._drop_zeros!(p) == p
        @test isequal(p.child.data.value, y)
        # test truncating leaf
        p = Node(NodeData(:^))
        n1 = addchild(p, NodeData(y))
        n2 = addchild(p, NodeData(:*))
        n3 = addchild(n2, NodeData(0))
        n4 = addchild(n2, NodeData(t))
        @test InfiniteOpt._drop_zeros!(p) == p
        @test p.child.sibling.data.value == 0
    end
    # test isequal for Nodes
    @testset "Base.isequal (Nodes)" begin 
        # test simple case 
        n1 = Node(NodeData(1))
        n2 = Node(NodeData(2))
        @test !isequal(n1, n2)
        n2.data = NodeData(1)
        @test isequal(n1, n2)
        # test unequal count case 
        n1 = Node(NodeData(:+))
        addchild(n1, NodeData(1))
        n2 = Node(NodeData(:+))
        addchild(n2, NodeData(1))
        addchild(n2, NodeData(1))
        @test !isequal(n1, n2)
        # test complicated isequal case 
        n1 = Node(NodeData(:+))
        addchild(n1, NodeData(1))
        addchild(n1, NodeData(y))
        n2 = Node(NodeData(:+))
        addchild(n2, NodeData(1))
        addchild(n2, NodeData(y))
        @test isequal(n1, n2)
        # test more complicated not equal case 
        n2 = Node(NodeData(:+))
        addchild(n2, NodeData(1))
        addchild(n2, NodeData(t))
        @test !isequal(n1, n2)
    end
    # test NLPExpr 
    @testset "NLPExpr" begin
        @test NLPExpr(Node(NodeData(0))).tree_root.data.value == 0
    end
    # test Base basics
    @testset "Base Basics (NLPExpr)" begin 
        # setup expr
        p = Node(NodeData(:+))
        addchild(p, NodeData(y))
        addchild(p, NodeData(2))
        nlp = NLPExpr(p)
        # test broadcastable
        @test Base.broadcastable(nlp) isa Base.RefValue
        # test copy
        @test copy(nlp) !== nlp
        @test isequal(copy(nlp), nlp)
        # test zero and one 
        @test zero(NLPExpr).tree_root.data.value == 0 
        @test one(NLPExpr).tree_root.data.value == 1
        # test isequal
        p = Node(NodeData(:+))
        addchild(p, NodeData(y))
        addchild(p, NodeData(3))
        nlp2 = NLPExpr(p)
        @test !isequal(nlp, nlp2)
        p.child.sibling.data = NodeData(2)
        @test isequal(nlp, nlp2)
    end
    # test drop_zeros!
    @testset "drop_zeros!" begin 
        # setup expr
        p = Node(NodeData(:^))
        n1 = addchild(p, NodeData(y))
        n2 = addchild(p, NodeData(:*))
        n3 = addchild(n2, NodeData(0))
        n4 = addchild(n2, NodeData(t))
        nlp = NLPExpr(p)
        @test drop_zeros!(nlp) === nlp
        @test nlp.tree_root.child.sibling.data.value == 0
    end
    # test isequal_canonical
    @testset "isequal_canonical" begin 
        # test that they are equal
        p = Node(NodeData(:^))
        n1 = addchild(p, NodeData(y))
        n2 = addchild(p, NodeData(:*))
        n3 = addchild(n2, NodeData(0))
        n4 = addchild(n2, NodeData(t))
        nlp1 = NLPExpr(p)
        p = Node(NodeData(:^))
        n1 = addchild(p, NodeData(y))
        n2 = addchild(p, NodeData(0))
        nlp2 = NLPExpr(p)
        @test isequal_canonical(nlp1, nlp2)
        # test that they are not equal
        p = Node(NodeData(:^))
        n1 = addchild(p, NodeData(y))
        n2 = addchild(p, NodeData(t))
        nlp3 = NLPExpr(p)
        @test !isequal_canonical(nlp1, nlp3)
        @test !isequal_canonical(nlp2, nlp3)
    end
    # test tree printing
    @testset "Tree Printing" begin
        # setup
        p = Node(NodeData(:^))
        n1 = addchild(p, NodeData(y))
        n2 = addchild(p, NodeData(2))
        nlp = NLPExpr(p)
        # test with NLPExpre
        expected = "^\n├─ y(t)\n└─ 2\n"
        io_test(print_expression_tree, expected, nlp)
        # test with other 
        io_test(print_expression_tree, "y(t)\n", y)
    end
    # test ast mapping 
    @testset "ast mapping" begin 
        @variable(Model(), x)
        # map tree of vars and constants 
        p = Node(NodeData(:^))
        n1 = addchild(p, NodeData(y))
        n2 = addchild(p, NodeData(:*))
        addchild(n2, NodeData(t))
        addchild(n2, NodeData(0.2))
        nlp = NLPExpr(p)
        @test map_nlp_to_ast(v -> x, nlp) == :($x ^ ($x * 0.2))
        # map tree with affine expression
        p = Node(NodeData(:^))
        n1 = addchild(p, NodeData(2y + t - 42))
        n2 = addchild(p, NodeData(:*))
        addchild(n2, NodeData(t))
        addchild(n2, NodeData(0.2))
        nlp = NLPExpr(p)
        @test map_nlp_to_ast(v -> x, nlp) == :((2 * $x + $x + -42) ^ ($x * 0.2))
        # map tree with quadratic expression (no affine)
        p = Node(NodeData(:^))
        n1 = addchild(p, NodeData(2y^2 + t^2))
        n2 = addchild(p, NodeData(:*))
        addchild(n2, NodeData(t))
        addchild(n2, NodeData(0.2))
        nlp = NLPExpr(p)
        @test map_nlp_to_ast(v -> x, nlp) == :((2 * $x * $x + $x * $x) ^ ($x * 0.2))
        # map tree with quadratic expression (w/ affine)
        p = Node(NodeData(:^))
        n1 = addchild(p, NodeData(2y^2 + t^2 - 42))
        n2 = addchild(p, NodeData(:*))
        addchild(n2, NodeData(t))
        addchild(n2, NodeData(0.2))
        nlp = NLPExpr(p)
        @test map_nlp_to_ast(v -> x, nlp) == :((2 * $x * $x + $x * $x + -42) ^ ($x * 0.2))
    end
end

# Test the basic expression generation via operators
@testset "Operator Definition" begin 
    # setup model data 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @variable(m, y, Infinite(t))
    @variable(m, z)
    # test _process_child_input
    @testset "_process_child_input" begin 
        # nlp 
        p = Node(NodeData(:sin))
        addchild(p, NodeData(y))
        nlp = NLPExpr(p)
        @test InfiniteOpt._process_child_input(nlp) == p
        # expressions
        @test InfiniteOpt._process_child_input(y) == NodeData(y)
        @test isequal_canonical(InfiniteOpt._process_child_input(2y).value, 2y)
        @test isequal_canonical(InfiniteOpt._process_child_input(y^2).value, y^2)
        # function symbol 
        @test InfiniteOpt._process_child_input(:^) == NodeData(:^)
        # constant 
        @test InfiniteOpt._process_child_input(true) == NodeData(true)
        @test InfiniteOpt._process_child_input(3) == NodeData(3)
        # fallback 
        @test_throws ErrorException InfiniteOpt._process_child_input("bad")
    end
    # test the basic graph builder
    @testset "_call_graph" begin 
        # test simple 
        p = Node(NodeData(:*))
        addchild(p, NodeData(y))
        addchild(p, NodeData(t))
        addchild(p, NodeData(1))
        @test isequal(InfiniteOpt._call_graph(:*, y, t, 1), p)
        # test with other expression 
        p = Node(NodeData(:sin))
        addchild(p, NodeData(y))
        nlp = NLPExpr(p)
        p = Node(NodeData(:^))
        addchild(p, nlp.tree_root)
        addchild(p, NodeData(y^2))
        @test isequal(InfiniteOpt._call_graph(:^, nlp, y^2), p)
    end
end

# TODO carry on...
