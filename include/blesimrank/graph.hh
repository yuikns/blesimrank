#ifndef BLESIMRANK_GRAPH_HH
#define BLESIMRANK_GRAPH_HH

#include <cstring>

#include <vector>
#include <set>
#include <map>
#include <queue>  // queue
#include <string>

//#include "moodycamel/concurrentqueue.h"

// using moodycamel::ConcurrentQueue;

class GEdge {
public:
    size_t id;
    double weight;
    double accum_weight;
    GEdge(size_t id, double weight, double accum_weight)
        : id(id), weight(weight), accum_weight(accum_weight) {}
    // GEdge(double weight, double accum_weight) : weight(weight), accum_weight(accum_weight) {}
    ~GEdge() {}
};

class GNode {
public:
    GNode() {
        weight = 0.0;
        gamma = 0.0;
    }
    //std::map<size_t, GEdge> edges;
    std::vector<GEdge> edges;
    std::vector<size_t> path;
    // ConcurrentQueue<int> path_cq;
    // std::queue<int> path_q;
    double weight;
    //std::vector<double> gamma;
    double gamma;
};


class BGrapgh {
public:
    BGrapgh(BGrapgh *_g) : BGrapgh(_g->nsize()) {
        //Do I Need Copy The Edges? 
        //memcpy(_nodes, _g->_nodes, sizeof(GNode) * node_size);
    }

    BGrapgh(size_t node_size)
        : node_size(node_size),
          _nodes(new GNode[node_size]) {}

    ~BGrapgh() {
        delete[] _nodes;
    }

    GNode *get_node(size_t offset) {
        if (offset >= node_size) {
            return NULL;
        } else {
            return &(_nodes[offset]);
        }
    }

    size_t nsize() const { return node_size; }

    size_t esize() const {
        size_t edge_size = 0;
        for (int i = 0; i < node_size; i++) {
            // printf("[Graph] %d : %zu \n", i, _nodes[i].neighbors.size());
            // for (int j = 0; j < _nodes[i].neighbors.size(); j++) {
            //    printf("\t%d -> %d \n",i,_nodes[i].neighbors[j].id);
            //}
            edge_size += _nodes[i].edges.size();
        }
        return edge_size;
    }

    void copy(BGrapgh *_g) {}

private:
    size_t node_size;
    GNode *_nodes;
};

#endif  // BLESIMRANK_GRAPH_HH
