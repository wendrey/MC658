#ifndef READEDGEGRAPH_DEFINE
#define READEDGEGRAPH_DEFINE

#include<lemon/list_graph.h>
#include<lemon/math.h>
#include <lemon/concepts/digraph.h>
using namespace std;
using namespace lemon;

/* le um grafo a partir de um arquivo texto, formato arestas */
#define MAXLABELNAME 200
#define MAXLINE 1000

// read a list digraph. If go_and_back is true, for a line [u,v,cost] the 
// routine insert the arc (u,v) and (v,u), both with cost custo. Otherwise,
// it insert only the arc (u,v).
bool ReadListDigraph(ListDigraph &g,
		     ListDigraph::NodeMap<string>& nodename,
		     ListDigraph::ArcMap<double>& custo,
		     const bool go_and_back,
		     char *filename);

// Read a geometric graph (points in the euclidean plane) or a list graph
// If the graph is geometric, the positions (posx and posy) are the given points. 
// If the graph is a list graph, the positions are computed by the neato program.
bool ReadGraph(ListGraph &g,
	       ListGraph::NodeMap<string>& nodename,
	       ListGraph::EdgeMap<double>& custo,
	       ListGraph::NodeMap<double>& posx,
	       ListGraph::NodeMap<double>& posy,
	       char *filename);

#endif
