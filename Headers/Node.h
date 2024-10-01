#ifndef GRAFO_BASICO_NODE_H
#define GRAFO_BASICO_NODE_H

#include "../Headers/Edge.h"
#include "../Headers/Defines.h"

struct Node
{
    size_t _number_of_edges;
    size_t _id;
    float  _weight;
    Edge  *_first_edge;
    Node  *_next_node;
    Node  *_previous_node;

    std::unordered_map<size_t, Edge*> edges;
    int _cluster_id;
    int demanda;

    Edge* get_edge_to(size_t target_id) {
        Edge* current_edge = _first_edge;
        while (current_edge != nullptr) {
            if (current_edge->_target_id == target_id) {
                return current_edge; // Retorna a aresta encontrada
            }
            current_edge = current_edge->_next_edge;
        }
        return nullptr; // Retorna nullptr se não encontrar a aresta
    }

    //  GET

    size_t get_number_of_edges(){return this->_number_of_edges;};
    size_t get_id(){return this->_id;};
    size_t get_cluster_id(){return this->_cluster_id;};
    float get_weight(){return this->_weight;};
    Edge *get_first_edge(){return this->_first_edge;};
    Node *get_next_node(){return this->_next_node;};
    Node *get_previous_node(){return this->_previous_node;};

    size_t getDemanda() { return this->demanda; };

    std::unordered_map<size_t, Edge*> getEdges() {
        std::unordered_map<size_t, Edge*> edges;
        Edge* current = _first_edge; // Supondo que `_first_edge` seja o ponteiro para a primeira aresta
        while (current != nullptr) {
        edges[current->_target_id] = current;
        current = current->_next_edge; // Avança para a próxima aresta
    }
    return edges;
    }

    //SET

    void set_number_of_edges(size_t numEdges){this->_number_of_edges = numEdges;};
    void set_id(size_t id){this->_id = id;};
    void set_cluster_id(size_t cluster_id){this->_cluster_id = cluster_id;};
    void set_weight(float weight){this->_weight = weight;};
    void set_first_edge(Edge *firstEdge){this->_first_edge = firstEdge;};
    void set_next_node(Node *nextNode){this->_next_node = nextNode;};
    void set_previous_node(Node *previousNode){this->_previous_node = previousNode;};
    void setDemanda(size_t demanda) { this->demanda = demanda; };    

    void setEdges(const std::unordered_map<size_t, Edge*>& newEdges) {
        this->edges = newEdges;
    }

    void insertEdge(int target_id, float weight) 
    {
    weight = 0;
    if(edges.count(target_id) == 0) {
        Edge* e = new Edge();
        edges.insert({target_id, e});
    }
    else
        //std::cout << "Aviso: " << _id << " -> " << target_id << " ja existe." << std::endl;
        return;
    }
    };

#endif  //GRAFO_BASICO_NODE_H
