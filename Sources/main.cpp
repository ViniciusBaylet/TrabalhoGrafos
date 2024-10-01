#include "../Headers/Graph.h"
#include "../Headers/Node.h"
#include "../Headers/Defines.h"

using namespace std;

Graph* readGraph(ifstream &input_file, bool directed, bool weighted_edge, bool weighted_node)
{
    Graph* graph;

    int order;
    int numEdges = 0;
    int node_id, node_weight;
    int source, target;        
    
    string line;
    
    // Skips lines
    for(int i = 0; i < 3; i++)
        getline(input_file, line);
    for(int i = 0; i < 4; i++)
        getline(input_file, line, ' ');

    // Reads number of clusters
    int clusters = stoi(line);    
    //cout << "Numero de subgrafos: " << clusters << endl;

    // Skips lines
    for(int i = 0; i < 2; i++)
        getline(input_file, line);
    for(int i = 0; i < 2; i++)
        getline(input_file, line, ' ');
    
    // Reads order of the graph
    order = stoi(line);
    //cout << "Numero de nos: " << order << endl;

    // Creates graph object
    graph = new Graph(order, numEdges, clusters, directed, weighted_edge, weighted_node);

    // Reads nodes and node weight
    for(int i = 0; i < 6; i++)
        getline(input_file, line);
    for(int i = 0; i < order; i++) 
    {
        input_file >> node_id >> node_weight;
        //graph->insertVertex(node_id, node_weight);
        graph->add_node(node_id, node_weight);
    }   

    int cont = 0;

    // Reads edges
    for(int i = 0; i < 5; i++)
        getline(input_file, line);
    for(int i = 0; i < order; i++) 
    {
        int start, end, pos = 0;
        getline(input_file, line);

        int edgesInLine = 0; 
                
        while(pos < line.size()) 
        {
            start = line.find('(', pos);
            if(start == string::npos)
                break;
            end = line.find(',', start);
            
            // Gets first node of the edge
            source = stoi(line.substr(start + 1, end));
            
            // Fix positions
            start = end;
            end = line.find(')', start);
            
            // Gets second node of the edge
            target = stoi(line.substr(start + 1, end));
            
            //graph->insertEdge(source, target);
            graph->add_edge(source, target, 0.0);
            numEdges++;
            cont++;

            // Fix position
            pos = end + 1;
        }  
        numEdges += edgesInLine;    
    }
    return graph;
}

int main(int argc, char* argv[]) {

    if (argc < 3) {
        cout
                << "ERROR: Expecting: ./<program_name> <input_file> <output_file> <Func> <param1> <maxIter>"
                << endl;
        return 1;
    }
    /*
    if (argc == 4 && stoi(argv[3]) != 1) {
        cout
                << "ERROR: Expecting: ./<program_name> <input_file> <output_file> <Func>"
                << endl;
        return 1;
    }
    
    // Le os argumentos de linha de comando
    string program_name(argv[1]);
    string input_file(argv[2]);

    string caminho_arquivo = "./Instances/" + input_file;
    */

    // Le os argumentos de linha de comando
    string input_file = argv[1];  // Nome do arquivo de entrada
    string output_file_name = argv[2];  // Nome do arquivo de saída

    string caminho_arquivo = "./Instances/" + input_file;

    // Abre o arquivo de entrada
    ifstream instance(caminho_arquivo);
    if (!instance.is_open()) {
        cerr << "Erro ao abrir o arquivo de entrada: " << caminho_arquivo << endl;
        return 1;
    }

    // Cria um objeto Graph
    Graph *graph;

    graph = readGraph(instance, false, false, true);

    // Fecha o arquivo de entrada
    instance.close();

    //* Guloso Randomizado Adaptativo
    float alfas_gra[3]        = {0.1, 0.2, 0.3};
    const size_t iterations_gra   = 1000;
    const size_t experiments_gra  = 30;

    //* Guloso Randomizado Adaptativo Reativo
    float alfas_grar[10]       = {0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500};
    //float alfas_grar[10]       = {0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.100};
    // float alfas_grar[10]       = {0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010};
    const size_t iterations_grar  = 4000;
    const size_t experiments_grar = 30;

    // Abre o arquivo de saída
    ofstream output_file;
    output_file.open(argv[2], ios::out | ios::trunc);

    // Menu de opcoes para o usuario
    int opcao;
    do {
        cout << "\nMenu de Opcoes:" << endl;
        cout << "1 - Algoritmo Guloso" << endl;
        cout << "2 - Algoritmo Guloso Randomizado Adaptativo" << endl;
        cout << "3 - Algoritmo Guloso Randomizado Adaptativo Reativo" << endl;
        cout << "0 - Sair" << endl;
        cout << "Digite a opcao desejada: ";
        cin >> opcao;

        switch (opcao) {
    case 1: {
        // Algoritmo Guloso
        cout << "Algoritmo Guloso para o grafo " << input_file << ": ";
        graph->algoritmoGuloso(0.0, output_file);
        cout << endl;
        
        break;
    }
    case 2: {
        // Algoritmo Gulos Randomizado Adaptativo
        cout << "Algoritmo Guloso Randomizado Adaptativo para o grafo " << input_file << ": ";
        for(int a = 0; a < experiments_gra; a++) {
        for(int i = 0 ; i < experiments_gra; i++) {
           graph->algoritmoGulosoRandomizadoAdaptativo(alfas_gra[a], iterations_gra, output_file);
         }
     }
        cout << endl;
        break;
    }
    case 3: {
        // Algoritmo Guloso Randomizado Adaptativo Reativo
        cout << "Algoritmo Guloso Randomizado Adaptativo Reativo para o grafo " << input_file << ": ";
        for(size_t i = 0 ; i < experiments_grar; i++) {
           graph->algoritmoGulosoRandomizadoAdaptativoReativo(alfas_grar, 10, iterations_grar, 40, output_file);
    }
        cout << endl;

        break;
    }
    case 0:
        cout << "Saindo do programa..." << endl;
        break;
    default:
        cout << "Opcao invalida!" << endl;
}
    } while (opcao != 0);

    // Pergunta ao usuário se deseja salvar a saída em um arquivo
    char choice;
    cout << "Deseja salvar a saida em um arquivo? (s/n): ";
    cin >> choice;

    if (choice == 's' || choice == 'S') {
        // Abre o arquivo e salva a saída
        if (output_file.is_open()) {
            cout << "Saida salva em arquivo_saida" << endl;
        } else {
            cout << "Erro ao abrir o arquivo_saida" << endl;
        }
    } else {
        cout << "Saida nao foi salva." << endl;
    }

    //Fecha o arquivo de saida
    output_file.close();

    return 0;
}
