#ifndef GRAFO_BASICO_EDGE_H
#define GRAFO_BASICO_EDGE_H

#include "../Headers/Defines.h"
#include "../Headers/Node.h"


struct Edge
{
    Edge  *_next_edge;
    float  _weight;
    size_t _target_id; // id de destino
    size_t _source_id; // id de origem
    //Node *_next_node;

    //Edge(size_t target_id, float weight) : _target_id(target_id), _weight(weight), _next_edge(nullptr) {}
     
    bool operator<(const Edge& other) const {
        return _weight < other._weight;
    }
};

#endif /* GRAFO_BASICO_EDGE_H */
