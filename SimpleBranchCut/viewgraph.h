#ifndef VIEWGRAPH_DEFINE
#define VIEWGRAPH_DEFINE

#include<lemon/list_graph.h>
#include <string>
#include "color.h"
using namespace std;
using namespace lemon;

#define VIEW_DOT 0
#define VIEW_NEATO 1

int ViewEuclideanGraph(ListGraph &g,
		       ListGraph::NodeMap<string> &vname, // nome dos vertices
		       ListGraph::NodeMap<double> &posx, // coord. x dos vertices
		       ListGraph::NodeMap<double> &posy, // coord. y dos vertices
		       ListGraph::NodeMap<int> &vcolor,  // cor dos vertices
		       ListGraph::EdgeMap<int> &acolor);  // cor das arestas


int ViewGraph(ListGraph &g,
	      int DOT_or_NEATO,
	      ListGraph::NodeMap<string> &vname, // nome dos vertices
	      ListGraph::EdgeMap<string> &ename,  // nome das arestas (e.g. peso dela)
	      ListGraph::NodeMap<double>& posx,
	      ListGraph::NodeMap<double>& posy,
	      ListGraph::NodeMap<int> &vcolor,   // cor dos vertices
	      ListGraph::EdgeMap<int> &acolor);    // cor das arestas

int ViewDigraph(ListDigraph &g,
	      int DOT_or_NEATO,
	      ListDigraph::NodeMap<string> &vname, // nome dos vertices
	      ListDigraph::ArcMap<string> &ename,  // nome das arestas (e.g. peso dela)
	      ListDigraph::NodeMap<int> &vcolor,   // cor dos vertices
	      ListDigraph::ArcMap<int> &acolor);   // cor das arestas

#endif
