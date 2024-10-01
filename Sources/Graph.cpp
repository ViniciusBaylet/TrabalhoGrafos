#include "../Headers/Defines.h"
#include "../Headers/Graph.h" 
#include "../Headers/Node.h"
#include "../Headers/Edge.h"

#define INT_INFINITY std::numeric_limits<int>::max()

using namespace std;

// Construtor da classe Graph
Graph::Graph(size_t _number_of_nodes, size_t _number_of_edge, size_t clusters, bool directed, bool weighted_edges, bool weighted_nodes)
{
    this->_number_of_nodes = _number_of_nodes;
    this->_clusters = clusters;
    this->_directed = directed;
    this->_weighted_edges = weighted_edges;
    this->_weighted_nodes = weighted_nodes;
    this->_number_of_edges = _number_of_edges;
}

// Destrutor da classe Graph
Graph::~Graph() {
    // Libera a memoria alocada para os nos e arestas
    Node* current_node = _first;
    while (current_node != nullptr) {
        Edge* current_edge = current_node->_first_edge;
        while (current_edge != nullptr) {
            Edge* next_edge = current_edge->_next_edge;
            delete current_edge;
            current_edge = next_edge;
        }
        Node* next_node = current_node->_next_node;
        delete current_node;
        current_node = next_node;
    }
}

// Getters para os atributos privados

size_t Graph::get_number_of_nodes()
{
    return this->_number_of_nodes;
}

size_t Graph::get_number_of_edges()
{
    return this->_number_of_edges;
}

bool Graph::getDirected(){
    return this->_directed;
}

bool Graph::getWeightedEdges(){
    return this->_weighted_edges;
}

bool Graph::getWeightedNodes(){
    return this->_weighted_nodes;
}

// Setters para os atributos privados

void Graph::set_number_of_edges(size_t numEdges){
   this->_number_of_edges = numEdges;
}

void Graph::setDirected(bool directed) {
   this->_directed = directed;
}

void Graph::setWeightedEdges(bool weighted_edges) {
    this->_weighted_edges = weighted_edges;
}

void Graph::setWeightedNodes(bool weighted_nodes) {
    this->_weighted_nodes = weighted_nodes;
}

void Graph::setOutfileName(std::string outfile_name)
{
    this->outfile_name = "output/" + outfile_name;
}

bool Graph::node_id_in_graph(size_t node_id) {
    Node* current_node = _first;

    while (current_node != nullptr) {
        if(current_node->_id == node_id) { 
            return true;
        }
        current_node = current_node->_next_node;
    }
    return false;
}

// Adiciona um no ao grafo
void Graph::add_node(size_t node_id, float weight) {
    // verifica se o no ja existe
    if(this->node_id_in_graph(node_id)) {
        cout << "O no node_id " << node_id << " ja pertence ao grafo. " << endl;
        return;
    };

    // Cria um novo no
    Node* new_node = new Node();
    new_node->_id = node_id;
    new_node->_weight = weight;
    new_node->_first_edge = nullptr;
    
    vertices[node_id] = new_node;

    // Se a lista de nos estiver vazia
    if (_first == nullptr) {
        _first = new_node;
        _last = new_node;
        new_node->_next_node = nullptr;
        new_node->_previous_node = nullptr;
    } else {
        // Insere o novo no no final da lista
        _last->_next_node = new_node;
        new_node->_previous_node = _last;
        _last = new_node;
        new_node->_next_node = nullptr;
    }
    mat_adj[node_id] = std::unordered_map<size_t, size_t>();
    //cout << "No com id " << node_id << " e peso " << weight << " adicionado com sucesso!" << endl;
}

Node* Graph::get_by_id(size_t node_id_1) {
    for(Node* node = this->_first; node != nullptr; node = node->_next_node) {
        if(node->_id == node_id_1) {
            return node;
        }
    }
    return nullptr;
}

// Adiciona uma aresta ao grafo
void Graph::add_edge(size_t node_id_1, size_t node_id_2, float weight) {
    // Encontra os nos que representam os vertices
    Node* node_1 = this->get_by_id(node_id_1);
    Node* node_2 = this->get_by_id(node_id_2);

    // Se os nos foram encontrados, cria uma nova aresta
    if (node_1 != nullptr && node_2 != nullptr) {
        // Cria uma nova aresta
        Edge* new_edge = new Edge();
        new_edge->_target_id = node_id_2;
        new_edge->_weight = weight;

        // Insere a nova aresta na lista de adjacencia do primeiro no
        new_edge->_next_edge = node_1->_first_edge;
        node_1->_first_edge = new_edge;

        // Atualiza a matriz de adjacência
        mat_adj[node_id_1][node_id_2] = weight;

        vertices[node_id_1]->insertEdge(node_id_2, weight);
        
        _number_of_edges++;

        // Se o grafo nao e direcionado, tambem adiciona a aresta inversa
        if(!_directed) {
            Edge* new_edge_inverse = new Edge();
            new_edge_inverse->_target_id = node_id_1;
            new_edge_inverse->_weight = weight;
            new_edge_inverse->_next_edge = node_2->_first_edge;
            node_2->_first_edge = new_edge_inverse;

            // Atualiza a matriz de adjacência para a aresta inversa
            mat_adj[node_id_2][node_id_1] = weight;

            vertices[node_id_1]->insertEdge(node_id_2, weight);
            //vertices[node_id_2]->insertEdge(node_id_1, weight);

            _number_of_edges++;

        }

        node_1->_number_of_edges++;
        if (!_directed) {
            node_2->_number_of_edges++;
        }
        
        //cout << "Aresta do no " << node_id_1 << " para o " << node_id_2 << " com peso " << weight << " adicionada!" << endl;
    }
    vector<Edge> edges = getEdges();
    edges.push_back({ nullptr, weight, node_id_2, node_id_1 });
}

// Remove um no do grafo
void Graph::remove_node(size_t node_id) {
    // Encontra o no a ser removido
    Node* current_node = _first;
    while (current_node != nullptr && current_node->_id != node_id) {
        current_node = current_node->_next_node;
    }

    // Se o no foi encontrado, remove-o da lista de nos
    if (current_node != nullptr) {
        // Remove as arestas do no
        Edge* current_edge = current_node->_first_edge;
        while (current_edge != nullptr) {
            Edge* next_edge = current_edge->_next_edge;
            delete current_edge;
            current_edge = next_edge;
        }

        // Remove o no da lista de nos
        if (current_node == _first) {
            _first = current_node->_next_node;
            if (_first != nullptr) {
                _first->_previous_node = nullptr;
            }
        } else if (current_node == _last) {
            _last = current_node->_previous_node;
            _last->_next_node = nullptr;
        } else {
            current_node->_previous_node->_next_node = current_node->_next_node;
            current_node->_next_node->_previous_node = current_node->_previous_node;
        }
        // Atualiza o numero de nos e de arestas
        _number_of_nodes--;
        _number_of_edges -= current_node->_number_of_edges;
        delete current_node;
    }
}

// Remove uma aresta do grafo
void Graph::remove_edge(size_t node_id_1, size_t node_id_2) {
    // Encontra os nos que representam os vertices
    Node* node_1 = _first;
    while (node_1 != nullptr && node_1->_id != node_id_1) {
        node_1 = node_1->_next_node;
    }

    Node* node_2 = _first;
    while (node_2 != nullptr && node_2->_id != node_id_2) {
        node_2 = node_2->_next_node;
    }

    // Se os nos foram encontrados, remove a aresta
    if (node_1 != nullptr && node_2 != nullptr) {
        // Remove a aresta da lista de adjacencia do primeiro no
        Edge* previous_edge = nullptr;
        Edge* current_edge = node_1->_first_edge;
        while (current_edge != nullptr && current_edge->_target_id != node_id_2) {
            previous_edge = current_edge;
            current_edge = current_edge->_next_edge;
        }

        // Se a aresta foi encontrada, remove-a
        if (current_edge != nullptr) {
            if (previous_edge == nullptr) {
                node_1->_first_edge = current_edge->_next_edge;
            } else {
                previous_edge->_next_edge = current_edge->_next_edge;
            }
            delete current_edge;
            _number_of_edges--;
            node_1->_number_of_edges--;

            // Se o grafo nao e direcionado, remove a aresta inversa
            if (!_directed) {
                previous_edge = nullptr;
                current_edge = node_2->_first_edge;
                while (current_edge != nullptr && current_edge->_target_id != node_id_1) {
                    previous_edge = current_edge;
                    current_edge = current_edge->_next_edge;
                }

                if (current_edge != nullptr) {
                    if (previous_edge == nullptr) {
                        node_2->_first_edge = current_edge->_next_edge;
                    } else {
                        previous_edge->_next_edge = current_edge->_next_edge;
                    }
                    delete current_edge;
                    node_2->_number_of_edges--;
                }
            }
        }
    }
}

// Imprime o grafo na tela
void Graph::print_graph() {
    cout << " ----- IMPRIMINDO O GRAFO -----" << endl;
    cout << endl;
    cout << "Numero de vertices: " << _number_of_nodes << endl;
    cout << "Numero de arestas: " << _number_of_edges << endl;

    // Percorre a lista de nos
    Node* current_node = _first;
    while (current_node != nullptr) {
        // Imprime o ID do no
        cout << "Vertice: " << current_node->_id;
        if (_weighted_nodes) {
            cout << " (Peso: " << current_node->_weight << ")";
        }
        cout << endl;

        // Percorre a lista de arestas do no
        Edge* current_edge = current_node->_first_edge;
        while (current_edge != nullptr) {
            // Imprime o ID do vertice de destino
            cout << "  - Destino: " << current_edge->_target_id;
            if (_weighted_edges) {
                cout << " (Peso: " << current_edge->_weight << ")";
            }
            cout << endl;

            current_edge = current_edge->_next_edge;
        }
        current_node = current_node->_next_node;
    }
    cout << endl;
}

// Imprime o grafo em um arquivo
void Graph::print_graph(ofstream& output_file) {
    output_file << " ----- IMPRIMINDO O GRAFO -----" << endl;
    output_file << endl;
    output_file << "Numero de vertices: " << _number_of_nodes << endl;
    output_file << "Numero de arestas: " << _number_of_edges << endl;

    // Percorre a lista de nos
    Node* current_node = _first;
    while (current_node != nullptr) {
        // Imprime o ID do no
        output_file << "Vertice: " << current_node->_id;
        if (_weighted_nodes) {
            output_file << " (Peso: " << current_node->_weight << ")";
        }
        output_file << endl;

        // Percorre a lista de arestas do no
        Edge* current_edge = current_node->_first_edge;
        while (current_edge != nullptr) {
            // Imprime o ID do vertice de destino
            output_file << "  - Destino: " << current_edge->_target_id;
            if (_weighted_edges) {
                output_file << " (Peso: " << current_edge->_weight << ")";
            }
            output_file << endl;

            current_edge = current_edge->_next_edge;
        }
        current_node = current_node->_next_node;
    }
    output_file << endl;
}

// Verifica se dois nos sao conectados
bool Graph::conected(size_t node_id_1, size_t node_id_2) {
    // Encontra o no que representa o primeiro vertice
    Node* node_1 = _first;
    while (node_1 != nullptr && node_1->_id != node_id_1) {
        node_1 = node_1->_next_node;
    }

    // Se o primeiro no foi encontrado, percorre a lista de arestas
    if (node_1 != nullptr) {
        Edge* current_edge = node_1->_first_edge;
        while (current_edge != nullptr) {
            // Se o ID do vertice de destino e igual ao ID do segundo vertice, os nos sao conectados
            if (current_edge->_target_id == node_id_2) {
                return true;
            }
            current_edge = current_edge->_next_edge;
        }
    }

    // Se a aresta noo foi encontrada, os nos nao sao conectados
    return false;
}

std::set<size_t> Graph::getVertices() {
    std::set<size_t> vertices;
    Node* current = _first;
    while (current != nullptr) {
        vertices.insert(current->_id);
        current = current->_next_node;
    }
    return vertices;
}

std::vector<Edge> Graph::getEdges() {
    std::vector<Edge> edges;
    std::set<Edge> edge_set; // Usado para evitar duplicatas em grafos não direcionados

    Node* current_node = _first;
    while (current_node != nullptr) {
        Edge* current_edge = current_node->_first_edge;
        while (current_edge != nullptr) {
            if (_directed) {
                edges.push_back(*current_edge);
            } else {
                // Para grafos não direcionados, garantimos que cada aresta seja adicionada apenas uma vez
                Edge edge_to_add = *current_edge;
                if (edge_to_add._source_id > edge_to_add._target_id) {
                    std::swap(edge_to_add._source_id, edge_to_add._target_id);
                }
                edge_set.insert(edge_to_add);
            }
            current_edge = current_edge->_next_edge;
        }
        current_node = current_node->_next_node;
    }

    if (!_directed) {
        edges.assign(edge_set.begin(), edge_set.end());
    }

    return edges;
}


struct Candidate
{
    int cluster;
    int source_id;
    int target_id;
    int increase_gap;

    friend bool operator<(Candidate c1, Candidate c2) 
    {
        return c1.increase_gap < c2.increase_gap;
    }

    friend bool operator>(Candidate c1, Candidate c2) 
    {
        return c1.increase_gap > c2.increase_gap;
    }

    friend bool operator==(Candidate c1, Candidate c2)
    {
        return c1.target_id == c2.target_id;
    }
}; 

typedef struct
{
    std::vector<int> vertices_ids;
    int higher;
    int lower;
} Solution;

struct AuxEdge
{
    int source;
    int target;
    int gap;

    friend bool operator<(AuxEdge e1, AuxEdge e2) 
    {
        return e1.gap < e2.gap;
    }

    friend bool operator>(AuxEdge e1, AuxEdge e2) 
    {
        return e1.gap > e2.gap;
    }
};

void removeByValue(std::vector<int> &v, int value) 
{
    int i;
    for(i = 0; i < v.size(); i++) {
        if(v[i] == value)
            break;
    }
    v.erase(v.begin() + i);
}

bool checkAvailable(std::vector<int> v, int a, int b) {
    int found = 0;
    for(int i = 0; i < v.size(); i++) {
        if(v[i] == a || v[i] == b)
            found++;
        if(found == 2) 
            return true;
    }
    return false;
}

void init_genrand(unsigned long s)
{
    static unsigned long mt[N];
    static int mti = N + 1;	

	mt[0] = s & 0xffffffffUL;
	for (mti = 1; mti < N; mti++)
	{
		mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
		mt[mti] &= 0xffffffffUL;
	}
}

unsigned long genrand_int32(void)
{
	unsigned long y;
	static unsigned long mag01[2] = {0x0UL, MATRIX_A};
	/* mag01[x] = x * MATRIX_A  for x=0,1 */

	if (mti >= N)
	{ 
		int kk;

		if (mti == N + 1)
			init_genrand(5489UL); 

		for (kk = 0; kk < N - M; kk++)
		{
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (; kk < N - 1; kk++)
		{
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
	}

	y = mt[mti++];

	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680UL;
	y ^= (y << 15) & 0xefc60000UL;
	y ^= (y >> 18);

	return y;
}


double genrand_real2(void) {
  return genrand_int32()*(1.0/4294967296.0);
  /* divided by 2^32 */
}

unsigned int intRandom(const unsigned int maxValue)
{
	#ifdef MTWISTER
		return ((double)genrand_real2()) * ((double)maxValue);
	#else
		static unsigned int res;
		static unsigned int rgen;
		static double factor;
		rgen = rand();
		factor = ((double)rgen / (double)INT_MAX);

		res = (unsigned int)(maxValue * factor);
		if (res == maxValue)
			res--;
		return res;
	#endif
}

void printAvailable(std::vector<int> v) 
{
    std::cout << "Indexes: \n   ";
    for(int i = 1; i < v.size(); i++)
    {
        std::cout << v[i] << " ";
    }
    std::cout << "\n";
}

void printSolutions(Solution solutions[], int clusters) 
{
    std::cout << "\nSolutions" << std::endl;
    for(int i = 0; i < clusters; i++) {
        printf("    %c = {", 'a'+i);   
        for(int j = 0; j < solutions[i].vertices_ids.size(); j++) {
            std::cout << solutions[i].vertices_ids[j];        
            if(j != solutions[i].vertices_ids.size() - 1)
                std::cout << ", ";
        }
        std::cout << "}" << "   Cost = " << (solutions[i].higher - solutions[i].lower)/2 << "\n\n";
    }
}

void printCandidates(std::list<Candidate> v){
    std::cout << "Candidates: \n";
    for(std::list<Candidate>::iterator it = v.begin(); it != v.end(); ++it)
    {
        std::cout << "  Candidate: " << it->target_id << ", cluster: " << it->cluster << ", increase in gap: " << it->increase_gap << "\n";
    }
    std::cout << "\n";
}

int Graph::algoritmoGuloso(float alpha, ofstream& output_file){

    std::vector<int> available;
    auto start = std::chrono::steady_clock::now();
    for(int i = 1; i <= _number_of_nodes; i++)
        available.push_back(i);
    
    Solution solutions[_clusters];

    std::list<Candidate> candidates;

    std::list<AuxEdge> all_edges;

    for (std::unordered_map<int, Node*>::iterator itV = vertices.begin(); itV != vertices.end(); ++itV) {
   
    Node *v = itV->second;
    
    // Certifique-se de que getEdges() está retornando um mapa válido
    std::unordered_map<size_t, Edge *> edges = v->getEdges();

    for (std::unordered_map<size_t, Edge*>::iterator itE = edges.begin(); itE != edges.end(); ++itE) {
        Edge *e = itE->second;
        AuxEdge aux;
        aux.source = v->get_id();
        aux.target = e->_target_id;
        aux.gap = abs(v->get_weight() - vertices[e->_target_id]->get_weight());
        all_edges.push_back(aux);
    }
}

// Certifique-se de que all_edges está sendo manipulado corretamente
all_edges.sort();

    int added = 0;
    std::list<AuxEdge>::iterator it;
    for (it = all_edges.begin(); it != all_edges.end(); ++it)
    {
        if(added < _clusters) 
        {
            if(checkAvailable(available, it->source, it->target)) 
            {
                solutions[added].vertices_ids.push_back(it->source);
                solutions[added].vertices_ids.push_back(it->target);
                if(vertices[it->source]->get_weight() > vertices[it->target]->get_weight()) {
                    solutions[added].higher = vertices[it->source]->get_weight();
                    solutions[added].lower = vertices[it->target]->get_weight();
                }
                else {
                    solutions[added].lower = vertices[it->source]->get_weight();
                    solutions[added].higher = vertices[it->target]->get_weight();
                }
                removeByValue(available,it->source);
                removeByValue(available,it->target);
                added++;
            }
        }
    }

     //printSolutions(solutions, _clusters);
     //printAvailable(available);
    for(int i = 0; i < _clusters; i++) 
    {       
        for(int j = 0; j < solutions[i].vertices_ids.size(); j++) 
        {
            int current_node_id = solutions[i].vertices_ids[j];
            std::unordered_map<size_t, Edge *> edges = vertices[current_node_id]->getEdges();
            for (std::unordered_map<size_t, Edge *>::iterator it = edges.begin(); it != edges.end(); ++it)
            {
                Edge *e = it->second;
                
                bool id_available = false;
                for(int i = 0; i < available.size(); i++) {
                    if(available[i] == e->_target_id)
                        id_available = true;
                }
                if(id_available) 
                {
                    Candidate c;
                    c.cluster = i;
                    c.source_id = current_node_id;
                    c.target_id = e->_target_id;
                    
                    int weigth = vertices[e->_target_id]->get_weight();
                    int lower_weigth = solutions[i].lower;
                    int higher_weigth = solutions[i].higher;
                    if(weigth < lower_weigth) {
                        c.increase_gap = abs(weigth - lower_weigth);
                    }
                    else if(weigth > higher_weigth) {
                        c.increase_gap = abs(weigth - higher_weigth);
                    }
                    else
                        c.increase_gap = 0;
                                        
                    candidates.push_back(c);
                }
            }
        }        
    }
     candidates.sort();

for (const auto& c : candidates) {
    //std::cout << "Candidato: " << c.target_id << " em cluster: " << c.cluster << std::endl;
}
    
    int n = 1;
    while(!available.empty())
    {
        // std::cout << "\n\nIteration: " << n << "\n\n";
        if (candidates.empty()) {
        //std::cerr << "Candidatos vazios, nao e possivel continuar." << std::endl;
        break;  // Saia do loop
        }

        candidates.sort();

        int candidates_index;
        if(alpha != 0) {
           candidates_index = rand() % candidates.size();
        }
        else
            candidates_index = 0;
    
        int chosen_vertex;
        int chosen_cluster;
        if(alpha != 0) {
            std::list<Candidate>::iterator it;
            int n = 0;
            for (it = candidates.begin(); it != candidates.end(); ++it){
                if(n == candidates_index) {
                    chosen_vertex = it->target_id;
                    chosen_cluster = it->cluster;
                    break;
                }
                n++;
            }
        }
        else {
            chosen_vertex = candidates.front().target_id;
            chosen_cluster = candidates.front().cluster;
        }

        solutions[chosen_cluster].vertices_ids.push_back(chosen_vertex);
        if(vertices[chosen_vertex]->get_weight() > solutions[chosen_cluster].higher)
            solutions[chosen_cluster].higher = vertices[chosen_vertex]->get_weight();
        else if(vertices[chosen_vertex]->get_weight() < solutions[chosen_cluster].lower)
            solutions[chosen_cluster].lower = vertices[chosen_vertex]->get_weight();
       
        removeByValue(available, chosen_vertex);
        Candidate aux;
        aux.target_id = chosen_vertex;
        candidates.remove(aux);

        std::unordered_map<size_t, Edge *> edges = vertices[chosen_vertex]->getEdges();
        for (std::unordered_map<size_t, Edge *>::iterator it = edges.begin(); it != edges.end(); ++it)
        {
            Edge *e = it->second;
            
            bool id_available = false;
            for(int i = 0; i < available.size(); i++) {
                if(available[i] == e->_target_id)
                    id_available = true;
            }
  
            if(id_available) 
            {
                Candidate c;
                c.cluster = chosen_cluster;
                c.source_id = chosen_vertex;
                c.target_id = e->_target_id;
                
                int weigth = vertices[e->_target_id]->get_weight();
                int lower_weigth = solutions[chosen_cluster].lower;
                int higher_weigth = solutions[chosen_cluster].higher;
                
                if(weigth < lower_weigth) {
                     c.increase_gap = abs(weigth - lower_weigth) + rand() % 5;
                }
                else if(weigth > higher_weigth) {
                    c.increase_gap = abs(weigth - higher_weigth);
                }
                else
                    c.increase_gap = 0;
            
                std::list<Candidate>::iterator it;
                int n = 0;
                for (it = candidates.begin(); it != candidates.end(); ++it){
                   if(it->increase_gap > c.increase_gap)
                      break;
                 }
                candidates.insert(it, c);
                candidates.push_back(c);
            }
        }
        n++;
    }   
    int cost = 0;
    for(int i = 0; i < _clusters; i++) {
        cost += (solutions[i].higher - solutions[i].lower);
    }
    //printSolutions(solutions, _clusters);
    auto end = std::chrono::steady_clock::now();
    std::cout << endl << "Custo total: " << cost << endl;
    std::cout << "Demorou " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() 
    << " ms para processar o algoritmo." << std::endl;
    return cost;
}

int Graph::algoritmoGulosoRandomizadoAdaptativo(float alpha, int iterations, ofstream& output_file){

    int seed = time(0);
    init_genrand(seed);
    int best = INT_INFINITY;
    int best_it;
    int aux;

    for(int i = 0; i < iterations; i++) {
        std::cout << "Iteracao " << i + 1 << std::endl;
        aux = algoritmoGuloso(alpha, output_file);
        if(aux < best) {
            best = aux;
            best_it = i;
        }
        break;
    }
    //std::cout << "Melhor custo: " << best << std::endl;
    imprimeGulosoRandomizadoAdaptativo(best, alpha, best_it, seed);
    return best;
}

void Graph::imprimeGulosoRandomizadoAdaptativo(float best_cost, float alfa, int best_it, int seed){
    ofstream outfile;
    outfile.open(outfile_name, std::ios::app);
    cout << "Melhor custo: " << best_cost << " alfa: " << alfa << " melhor iteracao: " << best_it << endl;
    outfile << best_cost << "," << alfa << "," << best_it << endl;
    outfile.close();
}

int chooseAlfa(float prob_alfa[]) 
{
    float rand = genrand_real2();
    float sum_prob = 0;
    for(int i = 0; i < 10; i++) {
        if(rand >= sum_prob && rand < sum_prob + prob_alfa[i])
            return i;
        if(i == 9)
            return 9;
        sum_prob += prob_alfa[i];
    }
    std::cout << "ERROR: Didn't found probability" << std::endl;
    return -1;
}

void updateProbabilities(unsigned long int V[], unsigned short int W[], float prob_alfa[], int best_cost, int delta) 
{
    float q[10];
    float sum_q = 0;
    for(int i = 0; i < 10; i++)
    {
        q[i] = pow((float)best_cost / ((float)V[i] / W[i]), delta);
        sum_q += q[i];
    }
    for (int i = 0; i < 10; i++)
    {
        prob_alfa[i] = q[i] / sum_q;
    }    
}

int Graph::algoritmoGulosoRandomizadoAdaptativoReativo(float alphas[], int tam_alpha, int iterations, int stack,ofstream& output_file){
    int seed = time(0);
    init_genrand(seed);

    int delta = 10;

    float* prob_alfas = new float [tam_alpha];
    unsigned long int* V = new unsigned long int [tam_alpha];
    unsigned short int* W = new unsigned short int [tam_alpha];

    for(int i = 0; i < tam_alpha; i++) {
        prob_alfas[i] = 1.0 / tam_alpha;
        V[i] = 0;
        W[i] = 0;
    }
    
    int best        = INT_INFINITY;
    int best_it     = 0;
    float best_alfa = 0;

    for(int i = 0; i < tam_alpha; i++) {
        V[i] = algoritmoGuloso(alphas[i], output_file);
        W[i] = 1;
    }

    int aux;
    int index_alpha;

    for(int i = 0; i < iterations; i++) {
        std::cout << "Iteration " << i + 1 << std::endl;
        
        if(i % stack == 0) {
            updateProbabilities(V, W, prob_alfas, best, delta);
        }

        index_alpha = chooseAlfa(prob_alfas);
        
        aux = algoritmoGuloso(alphas[index_alpha], output_file);

        int novoaux = aux/2;
        
        V[index_alpha] += novoaux;
        W[index_alpha] += 1;
        
        if(novoaux < best) {
            best = novoaux;
            best_it = i;
            best_alfa = alphas[index_alpha];
        }
    }
    imprimeAlgoritmoGulosoRandomizadoAdaptativoReativo(best, best_alfa, best_it, seed);
    std::cout << "Melhor custo encontrado: " << best << std::endl;

    return best;
}

void Graph::imprimeAlgoritmoGulosoRandomizadoAdaptativoReativo(float best_cost, float best_alfa, int best_it, int seed){
    ofstream outfile;
    outfile.open(outfile_name, std::ios::app);
    cout << "Melhor custo: " << best_cost << " melhor alfa: " << best_alfa << " melhor iteracao: " << best_it << endl; 
    outfile << best_cost << "," << best_alfa << "," << best_it << "," << seed << ",";
    outfile.close();
}

void Graph::printVertices(const std::unordered_map<int, Node*>& vertices) {
     // Cria um vetor para armazenar as chaves
    std::vector<int> keys;
    // Preenche o vetor com as chaves do unordered_map
    for (const auto& pair : vertices) {
        keys.push_back(pair.first);
    }
    // Ordena o vetor de chaves
    std::sort(keys.begin(), keys.end());
    // Imprime os elementos em ordem crescente
    for (int key : keys) {
        Node* node = vertices.at(key); // Obtem o Node correspondente à chave
        std::cout << "Chave: " << key << ", Node ID: " << node->get_id() << std::endl;
    }
}

void Graph::printEdges(const std::unordered_map<size_t, Edge*>& edges) {
    for (const auto& pair : edges) {
        int key = pair.first; // id da aresta
        Edge* edge = pair.second; // ponteiro para o objeto Edge
        std::cout << "Chave: " << key << ", Target ID: " << edge->_target_id << std::endl;
        // Adicione mais informações para imprimir, se necessário
    }

    cout << "Nao ha arestas para contar." << endl;
}

void Graph::insertEdge(int id, int target_id, float edge_weight, float source_vertex_weight, float target_vertex_weight)
{
    if (id == target_id)
    {
        std::cerr << "ERROR: Vertices cannot have equal ID! (Selfloop not allowed)" << id << std::endl;
        return;
    }

    if (id > _number_of_nodes)
    {
        std::cout << "ERROR: Changed the order, this is wrong!!!" << std::endl;
        _number_of_nodes = id;
    }
    if (target_id > _number_of_nodes)
    {
        std::cout << "ERROR: Changed the order, this is wrong!!!" << std::endl;
        _number_of_nodes = target_id;
    }

    if (this->_directed)
    {
        // if(vertices.find(id) != vertices.end()){
        //     std::cout << "Vertex already present" << std::endl;
        //     if(vertices[id]->getEdges().find(target_id) != vertices[id]->getEdges().end())
        //         std::cout << "Edge already present" << std::endl;
        // }
        vertices[id]->insertEdge(target_id, edge_weight);
    }
    else
    {
        // if(vertices.find(id) != vertices.end()){
        //     std::cout << "Vertex already present" << std::endl;
        //     if(vertices[id]->getEdges().find(target_id) != vertices[id]->getEdges().end())
        //         std::cout << "Edge already present" << std::endl;
        // }
        vertices[id]->insertEdge(target_id, edge_weight);
        // if(vertices.find(target_id) != vertices.end()){
        //     std::cout << "Vertex already present" << std::endl;
        //     if(vertices[target_id]->getEdges().find(id) != vertices[target_id]->getEdges().end())
        //         std::cout << "Edge already present" << std::endl;
        // }
        vertices[target_id]->insertEdge(id, edge_weight);
    }
}

