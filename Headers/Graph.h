#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "../Headers/Node.h"
#include "../Headers/Edge.h"
#include "../Headers/Defines.h"

using namespace std;

class Graph
{
public:

    // CONSTRUTOR E DESTRUTOR
    Graph(size_t _number_of_nodes,size_t _number_of_edges, size_t clusters, bool directed, bool weighted_edges, bool weighted_nodes);
    ~Graph();

    //GET
    size_t get_number_of_nodes();
    size_t get_number_of_edges();
    bool getDirected();
    bool getWeightedEdges();
    bool getWeightedNodes();

    //SET
    void set_number_of_edges(size_t numEdges);
    void setDirected(bool directed);
    void setWeightedEdges(bool weighted_edges);
    void setWeightedNodes(bool weighted_nodes);

    //FUNCOES DO GRAFO
    void remove_node(size_t node_id);
    void remove_edge(size_t node_id_1, size_t node_id_2);
    void add_node(size_t node_id, float weight = 0);
    void add_edge(size_t node_id_1, size_t node_id_2, float weight);
    void print_graph(std::ofstream& output_file);
    void print_graph();
    bool conected(size_t node_id_1, size_t node_id_2);
    set<size_t> getVertices();
    vector<Edge> edges; // Lista de arestas
    vector<Edge> getEdges();

    void setOutfileName(std::string outfile_name);

    int algoritmoGuloso(float alpha, ofstream& output_file);
    int algoritmoGulosoRandomizadoAdaptativo(float alpha, int iterations, ofstream& output_file);
    void imprimeGulosoRandomizadoAdaptativo(float best_cost, float alfa, int best_it, int seed);
    int algoritmoGulosoRandomizadoAdaptativoReativo(float alphas[], int tam_alpha, int iterations, int stack,ofstream& output_file);
    void imprimeAlgoritmoGulosoRandomizadoAdaptativoReativo(float best_cost, float best_alfa, int best_it, int seed);

    void printVertices(const std::unordered_map<int, Node*>& vertices);
    void printEdges(const std::unordered_map<size_t, Edge*>& edges);
    void insertEdge(int id, int target_id, float edge_weight = 1, float source_vertex_weight = 1, float target_vertex_weight = 1);

private:
    size_t _number_of_nodes;
    size_t _number_of_edges;
    size_t _clusters; //numero de subgrafos
    bool   _directed;
    bool   _weighted_edges;
    bool   _weighted_nodes;
    Node  *_first;
    Node  *_last;

    Node* get_by_id(size_t node_id);
    bool node_id_in_graph(size_t node_id);

    std::unordered_map<size_t, std::unordered_map<size_t, size_t>> mat_adj;

    std::unordered_map<int, Node*> vertices;
    std::string outfile_name;

};

#endif  //GRAPH_HPP
