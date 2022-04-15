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
        p = Node(NodeData(:abs))
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
        # test case with substraction of zero 
        p = Node(NodeData(:-))
        n1 = addchild(p, NodeData(:sin))
        n2 = addchild(p, NodeData(0))
        n3 = addchild(n1, NodeData(y))
        @test InfiniteOpt._drop_zeros!(p) == p
        @test p.data.value == :sin 
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
        # test again 
        test_output = @capture_out print_expression_tree(nlp)
        @test test_output == expected
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
    # test the sum 
    @testset "sum" begin 
        # test empty sums
        @test sum(i for i in Int[]) == 0
        @test isequal(sum(NLPExpr[]), zero(NLPExpr))
        @test isequal(sum(GeneralVariableRef[]), 
                      zero(GenericAffExpr{Float64, GeneralVariableRef}))
        # test NLP sums 
        p = Node(NodeData(:sin))
        addchild(p, NodeData(y))
        nlp = NLPExpr(p)
        new = NLPExpr(InfiniteOpt._call_graph(:+, nlp, nlp))
        @test isequal(sum(i for i in [nlp, nlp]), new)
        @test isequal(sum([nlp, nlp]), new)
        @test_throws ErrorException sum(i for i in [nlp, nlp]; bad = 42)
        @test_throws ErrorException sum(i for i in [y, t]; bad = 42)
        # test other expressions 
        @test isequal(sum(i for i in [y, t]), y + t)
        @test isequal(sum([y, t]), y + t)
        # normal sums 
        @test sum(i for i in 1:3) == 6
    end
    # test the product
    @testset "prod" begin 
        # test empty sums
        @test prod(i for i in Int[]) == 1
        @test isequal(prod(NLPExpr[]), one(NLPExpr))
        # test NLP sums 
        p = Node(NodeData(:sin))
        addchild(p, NodeData(y))
        nlp = NLPExpr(p)
        new = NLPExpr(InfiniteOpt._call_graph(:*, nlp, nlp))
        @test isequal(prod(i for i in [nlp, nlp]), new)
        @test isequal(prod([nlp, nlp]), new)
        @test_throws ErrorException prod(i for i in [nlp, nlp]; bad = 42)
        @test_throws ErrorException prod(i for i in [y, t, z]; bad = 42)
        # test other expressions 
        new = NLPExpr(InfiniteOpt._call_graph(:*, y, t, z))
        @test isequal(prod(i for i in [y, t, z]), new)
        @test isequal(prod([y, t, z]), new)
        # test normal products 
        @test prod(i for i in 1:3) == 6
    end
    # prepare the test grid 
    aff = 2z - 2
    quad = y^2 + y
    p = Node(NodeData(:sin))
    addchild(p, NodeData(y))
    nlp = NLPExpr(p)
    # test the multiplication operator
    @testset "Multiplication Operator" begin
        for (i, iorder) in [(42, 0), (y, 1), (aff, 1), (quad, 2), (nlp, 3)]
            for (j, jorder) in [(42, 0), (y, 1), (aff, 1), (quad, 2), (nlp, 3)]
                if iorder + jorder >= 3
                    expected = NLPExpr(InfiniteOpt._call_graph(:*, i, j))
                    @test isequal(i * j, expected)
                end
            end
        end
        expected = NLPExpr(InfiniteOpt._call_graph(:*, y, t, z,nlp))
        @test isequal(y * t * z * nlp, expected)
        @test isequal(*(nlp), nlp)
    end
    # test division operator 
    @testset "Division Operator" begin
        for i in [42, y, aff, quad, nlp]
            for j in [y, aff, quad, nlp]
                expected = NLPExpr(InfiniteOpt._call_graph(:/, i, j))
                @test isequal(i / j, expected)
            end
        end
        expected = NLPExpr(InfiniteOpt._call_graph(:/, nlp, 42))
        @test isequal(nlp / 42, expected)
        @test isequal(nlp / 1, nlp)
        @test_throws ErrorException nlp / 0
    end
    # test the power operator 
    @testset "Power Operator" begin
        for i in [42, y, aff, quad, nlp]
            for j in [y, aff, quad, nlp]
                expected = NLPExpr(InfiniteOpt._call_graph(:^, i, j))
                @test isequal(i ^ j, expected)
            end
        end
        one_aff = one(JuMP.GenericAffExpr{Float64, GeneralVariableRef})
        for i in [y, aff, quad, nlp]
            for f in [Float64, Int]
                @test isequal(i^zero(f), one_aff)
                @test isequal(i^one(f), i)
                if i isa NLPExpr
                    expected = NLPExpr(InfiniteOpt._call_graph(:^, i, f(2)))
                    @test isequal(i ^ f(2), expected)
                else
                    @test isequal(i^f(2), i * i)
                end
                expected = NLPExpr(InfiniteOpt._call_graph(:^, i, f(3)))
                @test isequal(i ^ f(3), expected)
            end
        end
        # extra tests 
        @test isequal(y^0, one_aff)
        @test isequal(y^1, y)
        @test isequal(y^2, y * y)
        @test isequal(y^0.0, one_aff)
        @test isequal(y^1.0, y)
        @test isequal(y^2.0, y * y)
    end
    # test the subtraction operator
    @testset "Subtraction Operator" begin
        for i in [42, y, aff, quad, nlp]
            expected = NLPExpr(InfiniteOpt._call_graph(:-, nlp, i))
            @test isequal(nlp - i, expected)
            expected = NLPExpr(InfiniteOpt._call_graph(:-, i, nlp))
            @test isequal(i - nlp, expected)
        end
        expected = NLPExpr(InfiniteOpt._call_graph(:-, nlp))
        @test isequal(-nlp, expected)
        @test isequal(y - y, zero(JuMP.GenericAffExpr{Float64, GeneralVariableRef})) 
        @test (y - z).constant == 0 
        @test (y - z).terms[y] == 1
        @test (y - z).terms[z] == -1
    end
    # test the addition operator 
    @testset "Addition Operator" begin
        for i in [42, y, aff, quad, nlp]
            expected = NLPExpr(InfiniteOpt._call_graph(:+, nlp, i))
            @test isequal(nlp + i, expected)
            expected = NLPExpr(InfiniteOpt._call_graph(:+, i, nlp))
            @test isequal(i + nlp, expected)
        end
        @test isequal(+nlp, nlp)
    end
    # test the comparison operators 
    @testset "Comparison Operators" begin 
        for i in [42, y, aff, quad, nlp]
            for j in [42, y, aff, quad, nlp]
                for (n, f) in (:< => Base.:(<), :(==) => Base.:(==), 
                               :> => Base.:(>), :<= => Base.:(<=), 
                               :>= => Base.:(>=))
                    if !(isequal(i, 42) && isequal(j, 42)) && !(isequal(i, y) && isequal(j, y))
                        expected = NLPExpr(InfiniteOpt._call_graph(n, i, j))
                        @test isequal(f(i, j), expected)
                    end
                end
            end
        end
        for (n, f) in (:(==) => Base.:(==), :<= => Base.:(<=), :>= => Base.:(>=))
            expected = NLPExpr(InfiniteOpt._call_graph(n, y, z))
            @test isequal(f(y, z), expected)
            @test f(y, y)
        end
        for (n, f) in (:< => Base.:(<), :> => Base.:(>))
            expected = NLPExpr(InfiniteOpt._call_graph(n, y, z))
            @test isequal(f(y, z), expected)
            @test !f(y, y)
        end
        # extra tests
        @test (y == 0) isa NLPExpr
        @test (y <= 0) isa NLPExpr
        @test (y >= 0) isa NLPExpr
        @test (y > 0) isa NLPExpr
        @test (y < 0) isa NLPExpr
        @test (0 == y) isa NLPExpr
        @test (0 <= y) isa NLPExpr
        @test (0 >= y) isa NLPExpr
        @test (0 > y) isa NLPExpr
        @test (0 < y) isa NLPExpr
    end
    # test the logic operators 
    @testset "Logic Operators" begin 
        for i in [y, nlp]
            for j in [y, nlp]
                for (n, f) in (:&& => Base.:&, :|| => Base.:|)
                    expected = NLPExpr(InfiniteOpt._call_graph(n, i, j))
                    @test isequal(f(i, j), expected)
                end
            end
        end
        for i in [y, nlp]
            @test !(i & false)
            @test !(false & i)
            @test isequal(i & true, i)
            @test isequal(true & i, i)
            @test (i | true)
            @test (true | i)
            @test isequal(i | false, i)
            @test isequal(false | i, i)
        end
    end
    # test ifelse 
    @testset "ifelse" begin
        for i in [42, y, aff, quad, nlp]
            for j in [42, y, aff, quad, nlp]
                expected = NLPExpr(InfiniteOpt._call_graph(:ifelse, nlp, i, j))
                @test isequal(InfiniteOpt.ifelse(nlp, i, j), expected)
            end
        end
        @test isequal(InfiniteOpt.ifelse(true, y, z), y)
        @test isequal(InfiniteOpt.ifelse(false, y, z), z)
    end
    # test the default functions 
    @testset "Default Registered Functions" begin 
        for i in [y, aff, quad, nlp]
            for (n, f) in InfiniteOpt._Base1ArgFuncList
            expected = NLPExpr(InfiniteOpt._call_graph(n, i))
            @test isequal(f(i), expected)
            end
        end
        for i in [y, aff, quad, nlp]
            for (n, f) in InfiniteOpt._Special1ArgFuncList
            expected = NLPExpr(InfiniteOpt._call_graph(n, i))
            @test isequal(f(i), expected)
            end
        end
    end
end

# Test the MutableArithmetics stuff
@testset "MutableArithmetics" begin 
    # setup model data 
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @variable(m, y, Infinite(t))
    @variable(m, z)
    aff = 2z - 2
    quad = y^2 + y
    nlp = sin(y)
    # test MA.mutability
    @testset "MA.mutability" begin 
        @test MA.mutability(NLPExpr) == MA.IsMutable()
    end
    # test MA.promote_operation
    @testset "MA.promote_operation" begin 
        for i in (2.0, 2, y, aff, quad, nlp)
            for f in (+, -, *, /, ^)
                @test MA.promote_operation(f, typeof(i), NLPExpr) == NLPExpr
                @test MA.promote_operation(f, NLPExpr, typeof(i)) == NLPExpr
            end
        end
        for i in (y, aff, quad)
            for f in (*, /, ^)
                @test MA.promote_operation(f, typeof(i), typeof(quad)) == NLPExpr
                @test MA.promote_operation(f, typeof(quad), typeof(i)) == NLPExpr
            end
        end
        for i in (2, 2.0, y, aff, quad, nlp)
            for j in (y, aff, quad, nlp)
                for f in (/, ^)
                    @test MA.promote_operation(f, typeof(i), typeof(j)) == NLPExpr
                end
            end
        end
    end
    # test MA.scaling 
    @testset "MA.scaling" begin
        @test_throws ErrorException MA.scaling(nlp)
        @test MA.scaling(one(NLPExpr)) == 1
    end
    # test mutable_copy
    @testset "MA.mutable_copy" begin
        @test MA.mutable_copy(nlp) === nlp
    end
    # test operate!
    @testset "MA.operate!" begin 
        # test zero and one 
        @test MA.operate!(zero, nlp) !== nlp
        @test isequal(MA.operate!(zero, nlp), zero(NLPExpr))
        @test MA.operate!(one, nlp) !== nlp
        @test isequal(MA.operate!(one, nlp), one(NLPExpr))
        # test operators 
        for i in (2.0, 2, y, aff, quad, nlp)
            for f in (+, -, *, /, ^)
                @test isequal(MA.operate!(f, nlp, i), f(nlp, i))
                @test isequal(MA.operate!(f, i, nlp), f(i, nlp))
            end
        end
        # test AddSubMul
        @test isequal(MA.operate!(MA.add_mul, nlp, y, 2), nlp + y * 2)
        @test isequal(MA.operate!(MA.sub_mul, nlp, y, 2), nlp - y * 2)
    end
end

# Test LinearAlgebra stuff 
@testset "Linear Algebra" begin 
    # setup the model data
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @variable(m, y, Infinite(t))
    @variable(m, z)
    aff = 2z - 2
    quad = y^2 + y
    nlp = sin(y)
    # test promotions 
    @testset "Base.promote_rule" begin 
        for i in (3, y, aff, quad)
            @test [nlp, i] isa Vector{NLPExpr}
        end
        # extra tests 
        @test promote_rule(NLPExpr, typeof(quad)) == NLPExpr
    end
    # test dot 
    @testset "LinearAlgebra.dot" begin 
        for i in (2, y, aff, quad, nlp)
            for j in (2, y, aff, quad, nlp)
                @test isequal(dot(i, j), i * j)
            end
        end
    end
    # test linear algebra operations 
    @testset "Operations" begin
        # make variables
        @variable(m, A[1:2, 1:3])
        @variable(m, x[1:2])
        @variable(m, w[1:3])
        # test expression
        @test x' * A * w isa NLPExpr # TODO fix with improved MA extension
    end
end

# Test registration utilities
@testset "Registration Methods" begin
    # setup model data 
    m = InfiniteModel()
    @variable(m, y)
    # test name_to_function
    @testset "name_to_function" begin 
        @test name_to_function(m, :tan, 1) == tan
        @test name_to_function(m, :+, 10) == +
        @test name_to_function(m, :*, 2) == *
    end
    # test all_registered_functions
    @testset "all_registered_functions" begin 
        @test all_registered_functions(m) isa Vector{Function}
    end
    # test user_registered_functions
    @testset "user_registered_functions" begin 
        @test user_registered_functions(m) == RegisteredFunction[]
    end
    # define functions for tests 
    f(a) = a^3
    g(a::Int) = 42
    h(a, b) = 42
    f1(a) = 32
    f2(a) = 10
    h1(a, b) = 13
    function hg(v, a, b)
        v[1] = 1
        v[2] = 2
        return 
    end
    # test creation helper errors 
    @testset "Registration Helpers" begin 
        @test_throws ErrorException RegisteredFunction(:a, 1, f, g)
        @test_throws ErrorException RegisteredFunction(:a, 2, f, g)
        @test_throws ErrorException RegisteredFunction(:a, 1, f, g, f)
        @test_throws ErrorException RegisteredFunction(:a, 1, f, f, g)
        @test_throws ErrorException InfiniteOpt._register(error, Main, m, :f, 1, 1)
        @test_throws ErrorException InfiniteOpt._register(error, Main, m, :sin, 1, sin)
        @test_throws ErrorException InfiniteOpt._register(error, Main, m, :g, 1, g)
        @test_throws ErrorException InfiniteOpt._register(error, Main, m, :eta, 1, eta)
    end
    # test @register
    @testset "@register" begin 
        # test errors 
        @test_macro_throws ErrorException @register(Model(), f(a))
        @test_macro_throws ErrorException @register(m, f(a), bad = 42)
        @test_macro_throws ErrorException @register(m, f(a), g, h, 3)
        @test_macro_throws ErrorException @register(m, f[i](a))
        @test_macro_throws ErrorException @register(m, f(a::Int))
        @test_macro_throws ErrorException @register(m, f(a), 2)
        @test_macro_throws ErrorException @register(m, g(a))
        @test_macro_throws ErrorException @register(m, eta(a))
        @test_macro_throws ErrorException @register(m, f(a), g)
        @test_macro_throws ErrorException @register(m, f(a), g, g)
        @test_macro_throws ErrorException @register(m, f(a), f, g)
        @test_macro_throws ErrorException @register(m, h(a, b), h)
        @test_macro_throws ErrorException @register(m, sin(a))
        # test univariate function with no gradient or hessian
        @test @register(m, f(a)) isa Function
        @test f(y) isa NLPExpr
        @test f(y).tree_root.data.value == :f
        @test isequal(f(y).tree_root.child.data.value, y)
        @test f(2) == 8
        @test length(user_registered_functions(m)) == 1
        @test user_registered_functions(m)[1].name == :f
        @test user_registered_functions(m)[1].func == f
        @test user_registered_functions(m)[1].num_args == 1
        @test user_registered_functions(m)[1].gradient isa Nothing 
        @test user_registered_functions(m)[1].hessian isa Nothing 
        @test name_to_function(m, :f, 1) == f
        # test univariate function with gradient
        @test @register(m, f1(a), f) isa Function
        @test f1(y) isa NLPExpr
        @test f1(y).tree_root.data.value == :f1
        @test isequal(f1(y).tree_root.child.data.value, y)
        @test f1(2) == 32
        @test length(user_registered_functions(m)) == 2
        @test user_registered_functions(m)[2].name == :f1
        @test user_registered_functions(m)[2].func == f1
        @test user_registered_functions(m)[2].num_args == 1
        @test user_registered_functions(m)[2].gradient == f
        @test user_registered_functions(m)[2].hessian isa Nothing 
        @test name_to_function(m, :f1, 1) == f1
        # test univariate function with gradient and hessian
        @test @register(m, f2(a), f, f1) isa Function
        @test f2(y) isa NLPExpr
        @test f2(y).tree_root.data.value == :f2
        @test isequal(f2(y).tree_root.child.data.value, y)
        @test f2(2) == 10
        @test length(user_registered_functions(m)) == 3
        @test user_registered_functions(m)[3].name == :f2
        @test user_registered_functions(m)[3].func == f2
        @test user_registered_functions(m)[3].num_args == 1
        @test user_registered_functions(m)[3].gradient == f
        @test user_registered_functions(m)[3].hessian == f1
        @test name_to_function(m, :f2, 1) == f2
        # test multivariate function with no gradient
        @test @register(m, h(a, b)) isa Function
        @test h(2, y) isa NLPExpr
        @test h(2, y).tree_root.data.value == :h
        @test isequal(h(2, y).tree_root.child.data.value, 2)
        @test isequal(h(2, y).tree_root.child.sibling.data.value, y)
        @test h(y, y) isa NLPExpr
        @test h(y, 2) isa NLPExpr
        @test h(2, 2) == 42
        @test length(user_registered_functions(m)) == 4
        @test user_registered_functions(m)[4].name == :h
        @test user_registered_functions(m)[4].func == h
        @test user_registered_functions(m)[4].num_args == 2
        @test user_registered_functions(m)[4].gradient isa Nothing
        @test user_registered_functions(m)[4].hessian isa Nothing
        @test name_to_function(m, :h, 2) == h
        # test multivariate function with gradient
        @test @register(m, h1(a, b), hg) isa Function
        @test h1(2, y) isa NLPExpr
        @test h1(2, y).tree_root.data.value == :h1
        @test isequal(h1(2, y).tree_root.child.data.value, 2)
        @test isequal(h1(2, y).tree_root.child.sibling.data.value, y)
        @test h1(y, y) isa NLPExpr
        @test h1(y, 2) isa NLPExpr
        @test h1(2, 2) == 13
        @test length(user_registered_functions(m)) == 5
        @test user_registered_functions(m)[5].name == :h1
        @test user_registered_functions(m)[5].func == h1
        @test user_registered_functions(m)[5].num_args == 2
        @test user_registered_functions(m)[5].gradient == hg
        @test user_registered_functions(m)[5].hessian isa Nothing
        @test name_to_function(m, :h1, 2) == h1
        # test wrong model error 
        @test_throws ErrorException f(@variable(InfiniteModel()))
        # test functional registration
        function registration_test()
            mt = InfiniteModel()
            @variable(mt, x)
            q(a) = 1
            @test @register(mt, q(a)) isa Function 
            @test q(x) isa NLPExpr
            return 
        end
        @test registration_test() isa Nothing 
        @test registration_test() isa Nothing 
    end
    # test add_registered_to_jump
    @testset "add_registered_to_jump" begin 
        # test normal 
        m1 = Model()
        @test add_registered_to_jump(m1, m) isa Nothing 
        @test m1.nlp_data.largest_user_input_dimension == 2
        r1 = m1.nlp_data.user_operators
        @test length(keys(r1.univariate_operator_to_id)) == 3
        @test r1.univariate_operator_f == [f, f1, f2]
        @test r1.univariate_operator_fprime[2:end] == [f, f]
        @test r1.univariate_operator_fprimeprime[3] == f1
        @test length(r1.multivariate_operator_to_id) == 2
        # test error 
        m2 = Model()
        h2(a, b) = 3
        @test @register(m, h2(a, b), hg, f) isa Function
        @test_throws ErrorException add_registered_to_jump(m2, m) isa Nothing 
        r2 = m2.nlp_data.user_operators
        @test length(keys(r2.univariate_operator_to_id)) == 3
        @test r2.univariate_operator_f == [f, f1, f2]
        @test r2.univariate_operator_fprime[2:end] == [f, f]
        @test r2.univariate_operator_fprimeprime[3] == f1
        @test length(r2.multivariate_operator_to_id) == 2
    end
end 

# Test string methods
@testset "String Methods" begin
    # setup the model data
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1])
    @variable(m, y, Infinite(t))
    @variable(m, z)
    # test making strings
    @testset "String Creation" begin
        # test some simple ones
        @test string(sin(y) + 2) == "sin(y(t)) + 2"
        @test string((z*z) ^ 4) == "(z²)^4"
        @test string((cos(z) + sin(z)) / y) == "(cos(z) + sin(z)) / y(t)"
        @test string(-cos(z + y) * z^2.3) == "-cos(z + y(t)) * z^2.3"
        @test string((-InfiniteOpt.ifelse(z == 0, 0, z)) ^ 3) == "(-(ifelse(z == 0, 0, z))^3"
        # test AffExpr cases 
        aff0 = zero(GenericAffExpr{Float64, GeneralVariableRef})
        @test string(aff0^4 + (2z - 3y + 42) ^ (-1z)) == "0^4 + (2 z - 3 y(t) + 42)^(-z)"
        @test string((1z) / (2.3 * y)) == "z / (2.3 y(t))"
        # test QuadExpr cases 
        quad0 = zero(GenericQuadExpr{Float64, GeneralVariableRef})
        @test string(quad0 / (z^2 + y * z)) == "0 / (z² + y(t)*z)"
        @test string((0z^2 + 42) * sin(y)) == "42 * sin(y(t))"
        @test string((z^2) / y) == "(z²) / y(t)"
    end
end
