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
bool bt(TSP_Data &tsp, int maxTime)
/*******************************************************************************
 * SUBSTITUIA O CONTEÚDO DESTE MÉTODO POR SUA IMPLEMENTAÇÃO DE BACKTRACKING.
 * ENTRETANTO, NÃO ALTERE A ASSINATURA DO MÉTODO.
 ******************************************************************************/
{

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
		
/*		cerr << "-------------------" << endl;		
		cerr << "BFS : " << visit << endl;
		cerr << "Node u : " << tsp.g.id(u) << endl;					
		cerr << "Node v : " << tsp.g.id(v) << endl;
		cerr << "Arc u->v : " << tsp.g.id(e) << endl;
		cerr << "Weight : " << tsp.weight[e] << endl;
		cerr << "Map v : " << node[v] << endl;		
*/				
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
bool bnb(TSP_Data &tsp,  int maxTime)
/*******************************************************************************
 * SUBSTITUIA O CONTEÚDO DESTE MÉTODO POR SUA IMPLEMENTAÇÃO DE BRANCH AND BOUND.
 * ENTRETANTO, NÃO ALTERE A ASSINATURA DO MÉTODO.
 ******************************************************************************/
{
	// Algoritmo guloso para o tsp
	
	EdgeBoolMap x(tsp.g);  // Maps a boolean x to each edge of the graph g
	for(ListGraph::EdgeIt e(tsp.g); e!=INVALID; ++e){  // We set every edge out of the the solution
		x[e] = false;
	}
	
	list<int> tour;
	double cost = 0.0;
	int nedges = 0;
	
	NodeBoolMap y(tsp.g);  // Maps a boolean y to each vertex of the graph g
	for(NodeIt o(tsp.g); o!=INVALID; ++o){
		y[o] = false;
		cerr << tsp.g.id(o) << ":" << y[o] << "  ";
	}
	cerr << endl;

	NodeIt nit(tsp.g);
	Node n = nit;
	Node f = nit;
	
	while(nedges != tsp.NNodes){
		y[n] = true;  // Put the vertex in the solution
		tour.push_back(tsp.g.id(n));
		
		cerr << "n: " << tsp.g.id(n) << "  y[n]: " << y[n] << endl;
		
		double wmin = DBL_MAX;  // min weight
		IncEdgeIt emin = INVALID;  // min inc edge of n
		Node nmin = INVALID;
		
		IncEdgeIt e(tsp.g, n);
		Node op = INVALID;
		
		cerr << "wmin: " << wmin << endl;
		
		for(; e != INVALID; ++e){

			
			op = tsp.g.v(e);
			if( op == n ){
				op = tsp.g.u(e);
			}
			
			cerr << "   (" << tsp.g.id(tsp.g.u(e)) << ", " << tsp.g.id(tsp.g.v(e)) << ")  x: " << x[e] << "  c: " << tsp.weight[e] << " op: " << tsp.g.id(op) << " y[op] " << y[op] << endl;
			
			if( ! y[ op ] ){        // The expression in [] returns the "destin" vertex of edge e
				if( tsp.weight[e] < wmin ){
					wmin = tsp.weight[e];
					emin = e;
					nmin = op;
				}
			}
		}
		
		if( wmin < DBL_MAX ){  // If got some edge
			cerr << "wmin: " << wmin << endl;
			x[emin] = true;  // Puts the edge e in the solution: this data will be visible outside this function
			nedges++;
			cost += wmin;
			n = nmin;
			cerr << "new n: " << tsp.g.id(n) << endl;
		}
		else{
			cerr << "Error: could not found a minimum weight value." << endl;
			exit(1);
		}
		
		if( nedges == tsp.NNodes - 1 ){
			y[f] = false;
		}
		
		cerr << "nedges: " << nedges << endl;
		cerr << endl;
	}

	
	if( nedges > 0 ){
		tsp.BestCircuitValue = cost;
	}

	return false;
}
//------------------------------------------------------------------------------
