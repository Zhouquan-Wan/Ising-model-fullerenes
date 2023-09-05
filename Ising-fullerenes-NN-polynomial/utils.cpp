#include "utils.h"

void graphEdgeIndexInit(Graph& g) {
    // initialize edge index
    edgeMap e_index = get(edge_index, g);
    graph_traits< Graph >::edges_size_type edge_count = 0;
    graph_traits< Graph >::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
        put(e_index, *ei, edge_count++);
}

dualGraphConstruction::dualGraphConstruction(Graph* _g, Graph* _dg, std::set<vPair>* _outerE) { 
    g = _g;
    dualg = _dg;
    outerE = _outerE;
    eInd = get(edge_index, *g);
    numE = num_edges(*g);
    edgeStart = std::vector<int>(numE, -1);
    edgeInclude = std::vector<bool>(numE, 1);
    treeVisited = std::vector<bool>(num_vertices(*dualg), 0);
    representingEdge = std::vector<vPair>(num_vertices(*dualg), vPair(0, 0));
    preV = -1;
    faceCount = 0;
    head = -1;
}

void dualGraphConstruction::meetEdge(Vertex v1, Vertex v2) {
    int ei = eInd[ edge(v1, v2, *g).first ];
    if (!edgeInclude[ei]) return;

    int curE = eInd[edge(v1, v2, *g).first];
    // std::cout << "meet edge " << curE << " " << v1 << " " << v2 << "\n";
    if(edgeStart[curE] == -1) { // unmet edge
        edgeStart[curE] = faceCount;
    } else {
        add_edge(faceCount, edgeStart[curE], *dualg);
        // std::cout << "new face edge " << faceCount << " " << edgeStart[curE] << "\n";
    }
}

void dualGraphConstruction::end_face() {
    meetEdge(head, preV);
    representingEdge[faceCount] = vPair(preV, head);
    if (outerE != NULL) {
        if(outerE->count(vPair(preV, head))) {
            root = faceCount;
        }
    }
    // std::cout << "face " << faceCount << " with edge " << preV << " " << head << "\n";
    faceCount ++;
    head = -1;
}

dualGraphTransform::dualGraphTransform(Graph& _g) {
    g = _g;
    planarize();
    eInd = get(edge_index, g);
    edgePool = std::vector<bool>(numE, 1);
    vertexPool = std::vector<bool>(numV, 1);
    vertexParity = std::vector<bool>(numV, 0);
    edgeWeight = std::vector<int>(numE, 1);
}

void dualGraphTransform::vertexSplit() {
    graph_traits<Graph>::vertex_iterator vi, vend;
    int outd, parity, vcur, vnew, vprev;
    for(tie(vi, vend) = vertices(g); vi != vend; vi++) {
        // std::cout << "visiting " << *vi << "\n";
        outd = out_degree(*vi, g);
        parity =  outd % 2;
        vertexParity.at(*vi) = parity;


        // if degree > 3, (degree - 3) vertices should be added,
        // and the original vertex be splitted
        if (outd > 3) {
            auto e = embedding[*vi].begin(); e+=2;
            vprev = *vi; 
            for (int j = 0; j < outd-3; j++) {
                // std::cout << j <<" " << "removing e " << *e << "\n";
                removeEdgeAndUpdate(*e);
                // remove_edge(*e, g);
                // edgePool.at(eInd[*e]) = 0;
                vcur = getAnotherVertex(*e, *vi);
                // std::cout << "vcur " << vcur << "\n";
                vnew = allocInd(vertexPool);
                vertexParity.push_back(0);
                addEdgeAndUpdate(vnew, vcur, 1);
                // auto newE = add_edge(vnew, vcur, g);
                // put(eInd, newE.first, allocInd(edgePool));

                addEdgeAndUpdate(vnew, vprev);
                // newE = add_edge(vnew, vprev, g);
                // put(eInd, newE.first, allocInd(edgePool));
                vprev = vnew;
                e++;
            }
            vcur = getAnotherVertex(*e, *vi);
            // std::cout << " " << "removing e " << *e << "\n";
            removeEdgeAndUpdate(*e);
            // remove_edge(*e, g);
            // edgePool.at(eInd[*e]) = 0;
            addEdgeAndUpdate(vnew, vcur, 1);
            // auto newE = add_edge(vnew, vcur, g);
            // put(eInd, newE.first, allocInd(edgePool));

        }
        // std::cout << "finish visiting " << *vi << "\n";
        // graph_traits< Graph >::edge_iterator ei, e_end;
        // for (tie(ei, e_end) = edges(g); ei != e_end; ei++) {
        //     std::cout << *ei << " index: " << eInd[*ei] << " weight: " << edgeWeight.at(eInd[*ei]) << "\n";
        // }
        // std::cout << " ========== \n";
        planarize();
    }

    vertexDeform();

    // for(tie(vi, vend) = vertices(g); vi != vend; vi++) {
    //     std::cout << *vi << "parity: " << vertexParity.at(*vi) << "\n";
    // }

}

int dualGraphTransform::allocInd(std::vector<bool>& pool) {
    for (auto cur = pool.begin(); cur != pool.end(); cur ++ ) {
        if (!(*cur)) {
            *cur = 1;
            return std::distance(pool.begin(), cur);
        }
    }
    pool.push_back(1);
    return pool.size() - 1;
}

void dualGraphTransform::planarize() {
    numV = num_vertices(g);
    numE = num_edges(g);
    embedding_storage = embedding_storage_t(numV);
    embedding = embedding_t(embedding_storage.begin(), get(vertex_index, g));
    boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
        boyer_myrvold_params::embedding = embedding);
}

int dualGraphTransform::getAnotherVertex(Edge e, int v) {
    if (source(e, g) != v) {
        return source(e, g);
    } else {
        return target(e, g);
    }
}

void dualGraphTransform::addEdgeAndUpdate(int vi, int v_end, int weight) {
    auto newE = add_edge(vi, v_end, g);
    int id = allocInd(edgePool);
    put(eInd, newE.first, id);
    // if new edge id, alloc new edgeWeight
    if (id == edgePool.size() - 1) {
        edgeWeight.push_back(weight);
    } else {
        edgeWeight.at(id) = weight;
    }
}

void dualGraphTransform::removeEdgeAndUpdate(Edge e) {
    remove_edge(e, g);
    // printf("here2!\n");
    edgePool.at(eInd[e]) = 0;
    edgeWeight.at(eInd[e]) = 0;
}

// int dualGraphTransform::removeEdgeGetWeight(int v1, int v2) {

// }

void dualGraphTransform::vertexDeform() {
    graph_traits<Graph>::vertex_iterator vi, vend;
    int adjV[3];
    int weight[3];
    int count;
    for(tie(vi, vend) = vertices(g); vi != vend; vi++) {
        // std::cout << "visiting " << *vi << "\n";
        count = 0;

        // for debug purpose
        // current program should only support those
        // Ising problems with no parallel edge
        assert(embedding[*vi].size() == 3);
        auto& edge3 = embedding[*vi];
        for (auto e = edge3.begin(); e != edge3.end(); e++) {
            // std::cout << "here! " << count << " " << *e << " " << eInd[*e] << "\n";
            weight[count] = edgeWeight.at(eInd[*e]);
            adjV[count] = getAnotherVertex(*e, *vi);
            if (count > 0) {
                removeEdgeAndUpdate(*e);
            }
            count ++;
        }

        if ( vertexParity.at(*vi) ) {
            int newv1 = allocInd(vertexPool);
            int newv2 = allocInd(vertexPool);
            addEdgeAndUpdate(adjV[1], newv1, weight[1]);
            addEdgeAndUpdate(adjV[2], newv2, weight[2]);
            addEdgeAndUpdate(newv1, newv2);
            addEdgeAndUpdate(*vi, newv1);
            addEdgeAndUpdate(*vi, newv2);
        } else {
            int newv[5];
            for(int i=0; i<5; i++) newv[i] = allocInd(vertexPool);
            addEdgeAndUpdate(adjV[1], newv[3], weight[1]);
            addEdgeAndUpdate(adjV[2], newv[4], weight[2]);
            addEdgeAndUpdate(newv[0], newv[1]); addEdgeAndUpdate(newv[1], newv[2]); addEdgeAndUpdate(newv[2], newv[0]);
            addEdgeAndUpdate(*vi, newv[0]);
            addEdgeAndUpdate(newv[3], newv[1]);
            addEdgeAndUpdate(newv[4], newv[2]);
        }

        planarize();

        // std::cout << "finish visiting " << *vi << "\n";
        // graph_traits< Graph >::edge_iterator ei, e_end;
        // for (tie(ei, e_end) = edges(g); ei != e_end; ei++) {
        //     std::cout << *ei << " index: " << eInd[*ei] << " weight: " << edgeWeight.at(eInd[*ei]) << "\n";
        // }
        // std::cout << num_edges(g) << " " << num_vertices(g) << " ========== \n";

    }
}

void alignPair(vPair& e) {
    if (e.first > e.second) {
        Vertex t = e.second;
        e.second = e.first;
        e.first = t;
    }
}

Vertex dualGraphTransform::getNext(Vertex v1, Vertex v2) {
    // std::cout << "query next of " << v1 << " " << v2 << "\n";
    auto e = embedding[v2].begin();
    Vertex v_start = getAnotherVertex(*e, v2);
    Vertex vcur;
    bool flag = false;
    for(; e != embedding[v2].end(); e++) {
        vcur = getAnotherVertex(*e, v2);
        if(flag) {
            return vcur;
        }

        if(vcur == v1) {
            flag = true;
            continue;
        }
    }
    // for debug purpose
    assert(flag);

    return v_start;
}

struct coord_t
{
    std::size_t x;
    std::size_t y;
};

void dualGraphTransform::getOuterFace() {
    gMaximal = g;
    embedding_storage_t emb_s = embedding_storage;
    embedding_t emb = embedding_t(emb_s.begin(), get(vertex_index, gMaximal));
    boyer_myrvold_planarity_test(boyer_myrvold_params::graph = gMaximal,
        boyer_myrvold_params::embedding = emb);
    make_biconnected_planar(gMaximal, emb);
    graphEdgeIndexInit(gMaximal);
    make_maximal_planar(gMaximal, emb);
    boyer_myrvold_planarity_test(boyer_myrvold_params::graph = gMaximal,
        boyer_myrvold_params::embedding = emb);
    graphEdgeIndexInit(gMaximal);
    // g_output(gMaximal);
    std::vector< graph_traits< Graph >::vertex_descriptor > ordering;
    planar_canonical_ordering(gMaximal, emb, std::back_inserter(ordering));

    // Set up a property map to hold the mapping from vertices to coord_t's
    typedef std::vector< coord_t > straight_line_drawing_storage_t;
    typedef boost::iterator_property_map<
        straight_line_drawing_storage_t::iterator,
        property_map< Graph, vertex_index_t >::type >
        straight_line_drawing_t;

    straight_line_drawing_storage_t straight_line_drawing_storage(
        num_vertices(gMaximal));
    straight_line_drawing_t straight_line_drawing(
        straight_line_drawing_storage.begin(), get(vertex_index, gMaximal));

    // Compute the straight line drawing
    chrobak_payne_straight_line_drawing(
        gMaximal, emb, ordering.begin(), ordering.end(), straight_line_drawing);

    // std::cout << "The straight line drawing is: " << std::endl;
    // graph_traits< Graph >::vertex_iterator vi, vi_end;
    // for (boost::tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
    // {
    //     coord_t coord(get(straight_line_drawing, *vi));
    //     std::cout << *vi << " -> (" << coord.x << ", " << coord.y << ")"
    //               << std::endl;
    // }

    graph_traits< Graph >::vertex_iterator vi, vi_end;
    Vertex vLeftmost;
    coord_t left;
    left.x = 1e9;
    for (tie(vi, vi_end) = vertices(gMaximal); vi != vi_end; ++vi)
    {
        coord_t coord(get(straight_line_drawing, *vi));
        if (coord.x < left.x) {
            left = coord;
            vLeftmost = *vi;
        }
    }

    // std::cout << "Left most vertex " << vLeftmost << "\n" ;

    graph_traits<Graph>::adjacency_iterator ai, a_end;
    Vertex o;
    float slope = -1e9;
    for(tie(ai, a_end) = adjacent_vertices(vLeftmost, g); ai != a_end; ai++) {
        coord_t coord(get(straight_line_drawing, *ai));
        // std::cout << *ai << "\n";
        if(coord.x == left.x) {
            o = *ai;
            break;
        } else {
            float tempSlope = float(coord.y - left.y) / float(coord.x - left.x);
            if (tempSlope > slope) {
                slope = tempSlope;
                o = *ai;
            }
        }
    }
    // std::cout << o << "\n";

    int v1 = vLeftmost;
    int v2 = o;
    int v_next = -1;
    for(; v_next != o; ) {
        // Edge now = edge(v1, v2, g).first;
        // std::cout << "visiting" << v1 << " " << v2 << "\n";
        outerEdges.insert(vPair(v1, v2));
        v_next = getNext(v1, v2);
        // printf("get %d \n", v_next);
        v1 = v2;
        v2 = v_next;
    }
}

void g_output(Graph& g) {
    std::cout << num_vertices(g) << " " << num_edges(g) << "\n";
    graph_traits< Graph >::edge_iterator ei, e_end;
    for (tie(ei, e_end) = edges( g ); ei != e_end; ei++) {
        std::cout << *ei << "\n";
    }
    std::cout <<"======\n";
}

int getTreeLeaf(dualGraphConstruction* dg, int root) {
    Graph* g = dg->dualg;
    int preV = -1;
    int cur = root;
    graph_traits<Graph>::adjacency_iterator ai, a_end;
    // int count = 0;
    while(true) {
        // count ++;
        // if (count > 10) return 0;
        // std::cout << "visiting" << cur << "\n";
        bool has_child = false;
        for(tie(ai, a_end) = adjacent_vertices(cur, *g); ai != a_end; ai++) {
            // std::cout << "child " << *ai << "\n";
            if ((*ai != preV) && (! dg->treeVisited[*ai] )) {
                preV = cur;
                cur = *ai;
                has_child = true;
                break;
            }
        }
        if(!has_child) {
            dg->treeVisited[cur] = true;
            // std::cout << "vertex removed " << cur << "\n";
            return cur;
        }
    }
}

void dualGraphTransform::directionAssignment() {
    std::map<Edge, int> eWeights;
    associative_property_map<std::map<Edge, int>> eWMap(eWeights);
    std::vector<Edge> spanning_tree;

    graph_traits< Graph >::edge_iterator ei, e_end;
    for (tie(ei, e_end) = edges( g ); ei != e_end; ei++) {
        eWeights.insert(std::pair<Edge, int>(*ei, 1));
    }

    kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree), weight_map(eWMap));

    // std::cout << "Print the edges in the MST:" << std::endl;
    // for (std::vector< Edge >::iterator ei = spanning_tree.begin();
    //      ei != spanning_tree.end(); ++ei)
    // {
    //     std::cout << source(*ei, g) << " <--> " << target(*ei, g)
    //               << " with weight of " << eWeights[*ei] << std::endl;
    // }

    faceCount mFaceCount;
    planar_face_traversal(g, embedding, mFaceCount);
    Graph tempdg(mFaceCount.faceNum);
    dualGraphConstruction dg(&g, &tempdg, &outerEdges);

    for (auto ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei) {
        dg.edgeInclude[eInd[*ei]] = false;
    }
    planar_face_traversal(g, embedding, dg);

    edgeDirection = std::vector<bool>(num_edges(g), 0);
    // g_output(*(dg.dualg));
    int leaf;
    // std::cout << "root is " << dg.root << "\n";
    do {
        leaf = getTreeLeaf(&dg, dg.root);
        // g_output(*(dg.dualg));
        vPair e = dg.representingEdge[leaf];
        int coOrientCount = 0;
        int v_next = -1;
        int v1 = e.first;
        int v2 = e.second;
        Edge now;
        int nowid, uid;
        bool ord;
        // std::cout << "visiting " << v1 << " " << v2 << " " << leaf << "\n";
        int notOrientCount = 0;
        for(; v_next != e.second;) {
            now = edge(v1, v2, g).first;
            nowid = eInd[now];
            // std::cout << nowid << "\n";
            if( ! dg.edgeInclude[nowid] ) {
                coOrientCount += int( (v1 < v2) != ( edgeDirection[nowid] ) );
            } else {
                // std::cout << now << "is to be determined\n";
                dg.edgeInclude[nowid] = false;
                uid = nowid;
                ord = v1 < v2;
                notOrientCount ++;
            }
            v_next = getNext(v1, v2);
            v1 = v2;
            v2 = v_next;
        }
        if (notOrientCount == 1) {
            if(coOrientCount % 2) {
                edgeDirection[uid] = ord;
            } else {
                edgeDirection[uid] = !ord;
            }
        }
        // std::cout << "=======\n";
    } while(leaf != dg.root);

    // graph_traits< Graph >::edge_iterator ei, e_end;
    // for (tie(ei, e_end) = edges(g); ei != e_end; ei++) {
    //         std::cout << *ei << " visited: " << eVisit[eInd[*ei]] << " orient: " << edgeDirection.at(eInd[*ei]) << "\n";
    // }
        
}