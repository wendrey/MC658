
//----------------------------------------------------------------------
// Example of program that uses the LEMON Library
//
// This program find a minimum cut that separates two verticex. 
// See definition (in portuguese) of st-min cut in page 43 of the slides
// in http://www.ic.unicamp.br/~fkm/lectures/proglin.pdf
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

typedef ListDigraph::ArcMap<double> CapMap;
typedef ListDigraph::NodeMap<string> nodenametype;
typedef ListDigraph::NodeMap<bool> CutMap;
typedef Preflow<ListDigraph, CapMap> PType;
typedef ListDigraph::Node Node;
typedef ListDigraph::NodeIt NodeIt;
typedef ListDigraph::Arc Arc;


// Let's call the mincut routine with a more intuitive name and parameteres
double MinCut(ListDigraph &g, CapMap &weight, Node &s,Node &t, CutMap &vcut)
{
  PType preflow_test(g, weight, s, t); 
  preflow_test.run();
  preflow_test.minCutMap(vcut);
  return (preflow_test.flowValue());
}


int main()
{

  // Declare an oriented graph
  ListDigraph g;
  CutMap cut(g); // cut is a boolean value for each vertex (see def. of CutMap above)
  CapMap weight(g); // weight is a double value for each arc (see def. above)
  nodenametype nodename(g); // nodename is a string for each vertex
  ListDigraph::ArcMap<string> ename(g); // each arc has also a name (for visualization purposes)
  Node s,t; // used to define the s-t cut.

  // To visualize the graph, each vertex/arc will have a color and label
  ListDigraph::NodeMap<int> vcolor(g);
  ListDigraph::NodeMap<string> vname(g);
  ListDigraph::ArcMap<int> acolor(g);
  ListDigraph::ArcMap<string> aname(g);

  // this program read the graph gr_50. The "true" parameter means that
  // for each edge of the file, the graph will have two oriented arcs, 
  // the backward and forward arc, both with the same weight.
  ReadListDigraph(g,nodename,weight,true,"gr_50");

  // set the label of each arc a (ename[a]) as the weight of the arc
  for (ListDigraph::ArcIt a(g); a!=INVALID; ++a) {
    std::ostringstream oss;
    oss << setprecision(2) << weight[a];
    acolor[a] = BLACK;
    ename[a]=oss.str(); 
  }

  // Well,... I know that the vertices "2" and "20" leads to a more interesting
  // cut. So, the two loops below looks for these vertices and set them as the
  // variables s and t.
  for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) 
    if (nodename[v]=="2") { s = v; break; };
  for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) 
    if (nodename[v]=="20") { t = v; break; };

  // Obtain a min cut separating vertices s and t
  double c = MinCut(g,weight,s,t,cut);

  // Color the arcs in the cut with color RED
  for (ListDigraph::ArcIt a(g); a!=INVALID; ++a) 
    if (cut[g.source(a)] && !cut[g.target(a)]) acolor[a]=RED;

  // color the nodes in one side of the cut with MAGENTA and the other 
  // set with color CYAN.
  for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) if (cut[v])  vcolor[v] = MAGENTA;
  for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) if (!cut[v]) vcolor[v] = CYAN;

  // print the cost of the cut
  cout << "Cut has value: " << c << endl;

  // visualize the graph using NEATO
  ViewDigraph(g,VIEW_NEATO,nodename,ename,vcolor,acolor);
  return 0;
}
