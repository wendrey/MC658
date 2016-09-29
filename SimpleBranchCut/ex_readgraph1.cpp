
//----------------------------------------------------------------------
// Example of program that uses the LEMON Library
//
// This program reads a graph and then show a pdf of the graph.
// The pdf of the graph is done by the graph neato, available in the package
// Graphviz. Of course, it uses a pdf reader.
//
// To compile this program, read the README file
//
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------

#include <iomanip>
#include <sstream>
#include <iostream>
#include "readgraph.h"
#include "viewgraph.h"
#include<lemon/list_graph.h>
#include <lemon/concepts/digraph.h>
#include <lemon/preflow.h>

int main()
{

  typedef ListDigraph::ArcMap<double> CapMap;
  typedef ListDigraph::NodeMap<string> nodenametype;
  typedef ListDigraph::Node Node;
  typedef ListDigraph::NodeIt NodeIt;
  typedef ListDigraph::Arc Arc;
  ListDigraph g;
  CapMap weight(g);

  ListDigraph::NodeMap<int> vcolor(g);
  ListDigraph::NodeMap<string> vname(g);
  ListDigraph::ArcMap<int> acolor(g);
  ListDigraph::ArcMap<string> ename(g);

  // Read the graph  gr_50. The parameter "true" makes the graph g with 
  // two arcs (forward and backward) for each edge in the file
  ReadListDigraph(g,vname,weight,true,"gr_50");

  // Color all vertices with color MAGENTA
  for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) vcolor[v] = MAGENTA;

  // color all edges with color BLUE and the label of each edge is its weight
  for (ListDigraph::ArcIt a(g); a!=INVALID; ++a) {
    std::ostringstream oss;
    oss << setprecision(2) << weight[a];
    acolor[a] = BLUE;
    ename[a]=oss.str(); 
  }
  
  // view the graph using the neato program.
  ViewDigraph(g,VIEW_NEATO,vname,ename,vcolor,acolor);
  return 0;
}
