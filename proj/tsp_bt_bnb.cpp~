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
#include <map>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "tsp_bt_bnb.h"

bool bt_dfs (TSP_Data &tsp, int maxTime, Node u, double cost, NodeBoolMap &node, clock_t t, vector<Node> circuit);
void updateSolution (TSP_Data &tsp, double cost, vector<Node> circuit);
bool bnb_dfs(TSP_Data &tsp, int maxTime, Node u, double cost, NodeBoolMap &node, clock_t t, vector<Node> circuit);
double getLowerBound (TSP_Data &tsp, NodeBoolMap &node);
bool mySort (const pair<Edge,double> &a, const pair<Edge,double> &b);

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
		return bt_dfs(tsp, maxTime, n, 0, node, t, circuit);
	return false;
	
}

bool bt_dfs (TSP_Data &tsp, int maxTime, Node u, double cost, NodeBoolMap &node, clock_t t, vector<Node> circuit) {

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
			bt_dfs(tsp, maxTime, v, cost + tsp.weight[e], node, t, circuit);
		
		// se achou o ciclo, verifica se a solucao melhora e atualiza
		
		else if (v == circuit.front() && circuit.size() == (unsigned) tsp.NNodes) 
			if (cost + tsp.weight[e] <= tsp.BestCircuitValue)
				;//updateSolution(tsp, cost + tsp.weight[e], circuit);
	
	}

	// tira o vertice da solucao

	circuit.pop_back();	
	node[u] = false;
	
	// retorna verdadeiro se encontrou uma solucao otima
	
	if (maxTime < (clock() - t) / CLOCKS_PER_SEC)
		return false;
	if (circuit.empty() && tsp.BestCircuit.size() == (unsigned) tsp.NNodes)
		return true;
	return false;

}

void updateSolution (TSP_Data &tsp, double cost, vector<Node> circuit) {

	// primeira solucao encontrada

	if (tsp.BestCircuit.empty()) {
		tsp.BestCircuit = circuit;
		tsp.BestCircuitValue = cost;
		return;
	}

	// acha o menor vertice na nova solucao

	int k = 0;
	
	for (unsigned int i = 1; i < circuit.size(); i++)
		if (tsp.vname[circuit[k]] > tsp.vname[circuit[i]])
			k = i;
			
	// verifica se tem menor ordem lexicografica

	if (cost == tsp.BestCircuitValue) {	
		for (int i = 0; i < tsp.NNodes; i++) {
			if (tsp.vname[tsp.BestCircuit[i]] > tsp.vname[circuit[(k+i)%tsp.NNodes]])
				break;
			if (tsp.vname[tsp.BestCircuit[i]] < tsp.vname[circuit[(k+i)%tsp.NNodes]])				
				return;
		}
	}

	// atualiza a solucao

	for (int i = 0; i < tsp.NNodes; i++)
		tsp.BestCircuit[i] = circuit[(k+i)%tsp.NNodes];
	tsp.BestCircuitValue = cost;
								
}

//------------------------------------------------------------------------------

bool bnb(TSP_Data &tsp, int maxTime) {

	clock_t t = clock();
	vector<Node> circuit;
	circuit.clear();
	tsp.BestCircuit.clear();
	NodeBoolMap node(tsp.g);
	EdgeBoolMap edge(tsp.g);

	for (ListGraph::NodeIt n(tsp.g); n != INVALID; ++n)
		node[n] = false;


	for (ListGraph::NodeIt n(tsp.g); n != INVALID; ++n)
		return bnb_dfs(tsp, maxTime, n, 0, node, t, circuit);
	return false;

}

bool bnb_dfs(TSP_Data &tsp, int maxTime, Node u, double cost, NodeBoolMap &node, clock_t t, vector<Node> circuit) {

	// verifica o tempo de execucao

	if (maxTime < (clock() - t) / CLOCKS_PER_SEC)
		return false;

	// poe o vertice na solucao 

	node[u] = true;
	circuit.push_back(u);

	double bound = getLowerBound(tsp, node);
	map <Edge,double> edges;

	// ordena as arestas que saem do vertice

	for (ListGraph::IncEdgeIt e(tsp.g, u); e != INVALID; ++e) 
		edges.insert(pair<Edge,double>(e,tsp.weight[e]));

	vector<pair<Edge,double>> ve(edges.begin(), edges.end());
	sort(ve.begin(), ve.end(), mySort);

	// passa pelos vizinhos do vertice

	for (unsigned int i = 0; i < ve.size(); i++) {

		Edge e = ve[i].first;
		Node v = tsp.g.u(e);
		if (v == u)
			v = tsp.g.v(e);
		
		// se existe uma potencial solucao, continua a busca
		
		if (node[v] == false && cost + bound < tsp.BestCircuitValue)
			bt_dfs(tsp, maxTime, v, cost + tsp.weight[e], node, t, circuit);
		
		// se achou o ciclo, verifica se a solucao melhora
		
		else if (v == circuit.front() && circuit.size() == (unsigned) tsp.NNodes) 
			if (cost + tsp.weight[e] <= tsp.BestCircuitValue)
				updateSolution(tsp, cost + tsp.weight[e], circuit);
	
	}
				
	// tira o vertice da solucao

	circuit.pop_back();	
	node[u] = false;
	
	// retorna verdadeiro se encontrou uma solucao otima
	
	if (maxTime < (clock() - t) / CLOCKS_PER_SEC)
		return false;
	if (circuit.empty() && tsp.BestCircuit.size() == (unsigned) tsp.NNodes)
		return true;
	return false;		
		
}

double getLowerBound (TSP_Data &tsp, NodeBoolMap &node) {

	double bound = 0;
	map<Node,list<Edge>> nemap;
	
	// acha as duas menores arestas que saem de cada vertice
	// considera apenas vertices que nao estao na solucao 
	
	for (ListGraph::NodeIt u(tsp.g); u != INVALID; ++u) {

		if (node[u])
			continue;

		for (ListGraph::IncEdgeIt e(tsp.g, u); e != INVALID; ++e) {

			if (nemap[u].size() == 0) 
				nemap[u].push_back(e);

			else if (nemap[u].size() == 1) {
				if (tsp.weight[nemap[u].front()] < tsp.weight[e])
					nemap[u].push_back(e);
				else
					nemap[u].push_front(e);
			}

			else if (tsp.weight[nemap[u].back()] > tsp.weight[e]) {
				nemap[u].pop_back();
				if (tsp.weight[nemap[u].front()] > tsp.weight[e])
					nemap[u].push_front(e);
				else
					nemap[u].push_back(e);				
			}

		}
	
	}
	
	// encontra o limitante inferior sa solucao otima
	
	for (auto i = nemap.begin(); i != nemap.end(); ++i)
		bound += tsp.weight[i->second.front()] + tsp.weight[i->second.back()];
	bound /= 2;
	
	return bound;
	
}

bool mySort (const pair<Edge,double> &a, const pair<Edge,double> &b) {

	return a.first < b.first;

}

//------------------------------------------------------------------------------
