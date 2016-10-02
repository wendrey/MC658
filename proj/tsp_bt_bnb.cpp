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
#include <string.h>
#include <iostream>
#include <float.h>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "tsp_bt_bnb.h"

void updateSolution (TSP_Data &tsp, double cost, vector<Node> circuit);
bool bfs (TSP_Data &tsp, int maxTime, Node u, double cost, NodeBoolMap &node, clock_t t, vector<Node> circuit);

//		cerr << "-------------------" << endl;		
//		cerr << "BFS : " << visit << endl;
//		cerr << "Node u : " << tsp.g.id(u) << endl;					
//		cerr << "Node v : " << tsp.g.id(v) << endl;
//		cerr << "Arc u->v : " << tsp.g.id(e) << endl;
//		cerr << "Weight : " << tsp.weight[e] << endl;
//		cerr << "Map v : " << node[v] << endl;		

//------------------------------------------------------------------------------

bool bt(TSP_Data &tsp, int maxTime) {

	clock_t t = clock();
	vector<Node> circuit;
	circuit.clear();
	tsp.BestCircuit.clear();
	NodeBoolMap node(tsp.g);

	for (ListGraph::NodeIt n(tsp.g); n != INVALID; ++n)
		node[n] = false;
			
	for (ListGraph::NodeIt n(tsp.g); n != INVALID; ++n)
		return bfs(tsp, maxTime, n, 0, node, t, circuit);

	return false;
	
}

bool bfs (TSP_Data &tsp, int maxTime, Node u, double cost, NodeBoolMap &node, clock_t t, vector<Node> circuit) {

	// verifica o tempo de execucao

	if (maxTime < (clock() - t) / CLOCKS_PER_SEC)
		return false;
		
	// poe o vertice na solucao e passa por todos seus vizinhos	

	node[u] = true;
	circuit.push_back(u);

	for (ListGraph::IncEdgeIt e(tsp.g, u); e != INVALID; ++e) {
		
		Node v = tsp.g.target(e);
						
		// se existe uma potencial solucao, continua a busca
		
		if (node[v] == false && cost + tsp.weight[e] < tsp.BestCircuitValue)
			bfs(tsp, maxTime, v, cost + tsp.weight[e], node, t, circuit);
		
		// se achou o ciclo, verifica se a solucao melhora
		
		else if (v == circuit.front() && circuit.size() == tsp.NNodes) 
			if (cost + tsp.weight[e] <= tsp.BestCircuitValue)
				updateSolution(tsp, cost + tsp.weight[e], circuit);
	
	}

	// tira o vertice da solucao

	circuit.pop_back();	
	node[u] = false;
	
	// retorna verdadeiro se encontrou uma solucao otima
	
	if (circuit.empty() && tsp.BestCircuit.size() == tsp.NNodes)
		return true;
	return false;

}

void updateSolution (TSP_Data &tsp, double cost, vector<Node> circuit) {

	// acha o menor vertice na nova solucao

	int k = 0;

	for (int i = 1; i < circuit.size(); i++)
		if (strcmp(tsp.vname[circuit[k]],tsp.vname[circuit[i]]) > 0)
			k = i;
			
	// verifica se tem melhor solucao

	if (cost == tsp.BestCircuitValue) {	
		for (int i = 0; i < tsp.NNodes; i++) {
			int aux = strcmp(tsp.vname[tsp.BestCircuit[i]],tsp.vname[circuit[(k+i)%tsp.NNodes]]);
			if (aux > 0)
				break;
			if (aux < 0)				
				return;
		}
	}

	// atualiza a solucao

	for (int i = 0; i < tsp.NNodes; i++)
		tsp.BestCircuit[i] = circuit[(k+i)%tsp.NNodes];
	tsp.BestCircuitValue = cost;
								
}

//------------------------------------------------------------------------------

bool bnb(TSP_Data &tsp,  int maxTime) {

	

	return false;

}
//------------------------------------------------------------------------------
