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
    # TODO finish
end

# TODO carry on...
