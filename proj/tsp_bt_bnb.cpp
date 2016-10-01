/*******************************************************************************
 * MC658 - Projeto e Análise de Algoritmos III - 2s2016
 * Prof: Flavio Keidi Miyazawa
 * PED: Mauro Henrique Mulati
 * Usa ideias e código de Rafael Arakaki e Flávio Keidi Miyazawa 
 ******************************************************************************/

/*******************************************************************************
 * EDITE ESTE ARQUIVO APENAS ONDE INDICADO
 * DIGITE SEU RA: 148234
 * SUBMETA SOMENTE ESTE ARQUIVO
 ******************************************************************************/

#include <time.h>
#include <iostream>
#include <float.h>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "tsp_bt_bnb.h"

bool bfs (TSP_Data &tsp, int maxTime, Node u, int visit, double cost, NodeBoolMap &node, clock_t t, vector<Node> circuit);

//------------------------------------------------------------------------------

bool bt(TSP_Data &tsp, int maxTime) {

	clock_t t = clock();
	vector<Node> circuit;

	NodeBoolMap node(tsp.g);
	for (ListGraph::NodeIt n(tsp.g); n != INVALID; ++n)
		node[n] = false;
	
	for (ListGraph::NodeIt n(tsp.g); n != INVALID; ++n)
		return bfs(tsp, maxTime, n, 0, 0, node, t, circuit);
	
	return false;
	
}

bool bfs (TSP_Data &tsp, int maxTime, Node u, int visit, double cost, NodeBoolMap &node, clock_t t, vector<Node> circuit) {

	node[u] = true;
	circuit.push_back(u);
		
	// dado um vertice, passa por todos seus vizinhos	

	for (ListGraph::IncEdgeIt e(tsp.g, u); e != INVALID; ++e) {
		
		// verifica o tempo de execucao

		if (maxTime < (clock() - t) / CLOCKS_PER_SEC)
			return false;

		Node v = tsp.g.target(e);
		
//		cerr << "-------------------" << endl;		
//		cerr << "BFS : " << visit << endl;
//		cerr << "Node u : " << tsp.g.id(u) << endl;					
//		cerr << "Node v : " << tsp.g.id(v) << endl;
//		cerr << "Arc u->v : " << tsp.g.id(e) << endl;
//		cerr << "Weight : " << tsp.weight[e] << endl;
//		cerr << "Map v : " << node[v] << endl;		
				
		// se existe uma potencial solucao, continua a busca
		// se achar uma solucao melhor, atualiza a solucao
		
		if (node[v] == false && cost + tsp.weight[e] < tsp.BestCircuitValue) {
			if (bfs(tsp, maxTime, v, visit+1, cost + tsp.weight[e], node, t, circuit))
				return true;	
		}
		
		// se achou o ciclo, verifica se a solucao melhora
		// se achar uma solucao melhor, atualiza a solucao
		
		else if (v == circuit.front() && visit+1 == tsp.NNodes) {
			if (cost + tsp.weight[e] < tsp.BestCircuitValue) {					
				tsp.BestCircuitValue = cost + tsp.weight[e];
				tsp.BestCircuit = circuit;
				return true;
			}
		}					
	
	}

	circuit.pop_back();	
	node[u] = false;
	return false;

}

//------------------------------------------------------------------------------

bool bnb(TSP_Data &tsp,  int maxTime) {

	

	return false;

}
//------------------------------------------------------------------------------
