#include <iostream>
#include <fstream>
#include "utils.h"

Graph constructIsingGraph(std::ifstream& f) {
    int vNum, eNum, v1, v2;
    f >> vNum >> eNum;
    Graph g(vNum);
    for(int i=0; i<eNum; i++) {
        f >> v1 >> v2;
        add_edge(v1, v2, g);
    }
    return g;
}

void graph_output(dualGraphTransform& d) {
    std::cout << num_vertices(d.g) << " " << num_edges(d.g) << "\n";
    graph_traits< Graph >::edge_iterator ei, e_end;
    for (tie(ei, e_end) = edges( d.g ); ei != e_end; ei++) {
        std::cout << source(*ei, d.g) << " " << target(*ei, d.g) << " " << d.edgeDirection.at(d.eInd[*ei]) << " " << d.edgeWeight[ d.eInd[*ei] ] << "\n";
    }
}

int main(int argc, char** argv)
{
    using namespace std;
    // cout << "reading " << argv[1] << "\n";
    ifstream f(argv[1]);
    Graph g = constructIsingGraph(f);
    graphEdgeIndexInit(g);
    // g_output(g);

    // planarize the original latice
    embedding_storage_t embedding_storage(num_vertices(g));
    embedding_t embedding(embedding_storage.begin(), get(vertex_index, g));
    bool is_planar = boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
        boyer_myrvold_params::embedding = embedding);
    assert(is_planar == true);
    
    // dual graph construction
    faceCount mFaceCount;
    planar_face_traversal(g, embedding, mFaceCount);
    // cout << mFaceCount.faceNum << endl;

    Graph dualg(mFaceCount.faceNum);
    dualGraphConstruction dC(&g, &dualg, NULL);
    planar_face_traversal(g, embedding, dC);
    graphEdgeIndexInit(dualg);

    // g_output(dualg);

    dualGraphTransform dGT(dualg);
    // cout << "here! \n";
    dGT.vertexSplit();
    // g_output(dGT.g);

    dGT.getOuterFace();

    dGT.directionAssignment();
    // cout << "here! \n";
    // dGT.graphBFS();
    // cout << "here! \n";
    graph_output(dGT);

    return 0;
}