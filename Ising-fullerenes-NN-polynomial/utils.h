#ifndef UTILS_H
#define UTILS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <vector>

#include <boost/graph/planar_canonical_ordering.hpp>
// #include <boost/graph/is_straight_line_drawing.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <set>
#include <queue>

#include <iostream>

using namespace boost;

typedef adjacency_list< vecS, vecS, undirectedS,
        property< vertex_index_t, int >, property< edge_index_t, int > >
        Graph;

typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;

typedef property_map< Graph, edge_index_t >::type edgeMap;

// Define the storage type for the planar embedding
typedef std::vector< std::vector< graph_traits< Graph >::edge_descriptor > >
    embedding_storage_t;
typedef boost::iterator_property_map< embedding_storage_t::iterator,
    property_map< Graph, vertex_index_t >::type >
    embedding_t;

typedef std::pair<Vertex, Vertex> vPair;

void graphEdgeIndexInit(Graph& g);

void g_output(Graph& g);

// count faces, utils for constructing a dual graph
struct faceCount: public planar_face_traversal_visitor {
    int faceNum = 0;    
    void begin_face() {faceNum ++;}
};

struct dualGraphConstruction: public planar_face_traversal_visitor {
    Graph* g;
    Graph* dualg;
    edgeMap eInd;
    int numE, faceCount;
    std::vector<int> edgeStart;
    std::vector<bool> edgeInclude;
    std::vector<vPair> representingEdge;
    std::vector<bool> treeVisited;
    Vertex preV, head;
    int root;
    std::set<vPair>* outerE;
    
    dualGraphConstruction(Graph* _g, Graph* _dg, std::set<vPair>* _outerE);

    void meetEdge(Vertex v1, Vertex v2);

    template <typename Vertex> void next_vertex(Vertex v) {
        if (head == -1) {
            head = v;
        } else {
            meetEdge(preV, v);
        }
        preV = v;
    }

    void end_face();
};

struct dualGraphTransform {
    Graph g;
    Graph gMaximal;
    int numV, numE;
    std::vector<bool> edgePool;
    std::vector<bool> vertexPool;
    // std::vector<bool> vertexColor;
    std::vector<bool> vertexParity;
    std::vector<int> edgeWeight;
    std::vector<bool> edgeDirection;
    std::set<vPair> outerEdges;

    // enum eStat {NOT, POS, NEG, BOTH};
    std::vector<int> eVisit;

    embedding_storage_t embedding_storage;
    embedding_t embedding;
    edgeMap eInd;
    dualGraphTransform(Graph& _g);

    void vertexSplit();
    void vertexDeform();

    int allocInd(std::vector<bool>& pool);

    void planarize();

    int getAnotherVertex(Edge e, int v);

    void addEdgeAndUpdate(int vi, int v_end, int weight=0);
    void removeEdgeAndUpdate(Edge e);

    void directionAssignment();

    void getOuterFace();

    // void graphBFS();
    Vertex getNext(Vertex v1, Vertex v2);
    // int removeEdgeGetWeight(int v1, int v2);

};

#endif