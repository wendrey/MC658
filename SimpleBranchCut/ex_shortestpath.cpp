
//----------------------------------------------------------------------
// Example of program that uses the LEMON Library, Linear Programming
//
// To use this program, you must install the LEMON Library and a linear programming
// package used by LEMON (like glpk, clp, cplex,...)
//
// This program find a shortest path in a oriented graph. Arc costs must 
// be non-negative. This is by solving a linear program. The formulation used
// is given in the page 48 of the slides in the link (in portuguese)
//
// http://www.ic.unicamp.br/~fkm/lectures/proglin.pdf
//
// To compile this program, read the README file
//
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------

#include <lemon/lp.h>
#include "color.h"
#include <iomanip>
#include <sstream>
#include <iostream>
#include "readgraph.h"
#include "viewgraph.h"
#include<lemon/list_graph.h>
#include <lemon/concepts/digraph.h>
#include <lemon/preflow.h>

#define EPS 0.00001
string DoubleToString(double x)
{
  std::ostringstream oss;
  oss << x;
  return(oss.str());
}


int main()
{
  ListDigraph g;
  ListDigraph::ArcMap<double> weight(g);
  ListDigraph::NodeMap<int> vcolor(g);
  ListDigraph::NodeMap<string> vname(g);
  ListDigraph::ArcMap<int> ecolor(g);
  ListDigraph::ArcMap<string> ename(g);
  ListDigraph::ArcMap<Lp::Col> x(g);
  Lp lp;
  Lp::Expr e;
  char filename[]="gr_caminho";
  ReadListDigraph(g,vname,weight,false,filename);
  for (ListDigraph::ArcIt a(g); a!=INVALID; ++a)  {
    x[a] = lp.addCol();
    e += weight[a]*x[a];
    ecolor[a] = BLUE;
    lp.colLowerBound( x[a], 0 );    lp.colUpperBound( x[a], 1 );
  }
  lp.obj(e);    lp.min();

  for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) {
    vcolor[v] = MAGENTA;
    Lp::Expr e;
    for (ListDigraph::OutArcIt a(g,v); a != INVALID; ++a)  e += x[a];
    for (ListDigraph::InArcIt  a(g,v); a != INVALID; ++a)  e -= x[a];
    if (vname[v] == "s")      lp.addRow(e == 1);
    else if (vname[v] == "t") lp.addRow(e == -1);
    else                      lp.addRow(e == 0);
  }
  lp.solve();

  if (lp.primalType() == Lp::OPTIMAL) {
    for (ListDigraph::ArcIt a(g); a!=INVALID; ++a) 
      if (lp.primal(x[a]) > 1.0-EPS) {
	ecolor[a] = RED;
	cout<<vname[g.source(a)]<<" , "<< vname[g.target(a)]<<endl;
      }
    // pinta as arestas de azul e coloca o peso da aresta como label
    for (ListDigraph::ArcIt a(g); a!=INVALID; ++a) {
      std::ostringstream oss;
      oss << setprecision(2) << weight[a];
      ename[a]=oss.str(); 
    }
    ViewDigraph(g,VIEW_NEATO,vname,ename,vcolor,ecolor);
  }else
    cout << "Erro na resolucao do Programa linear.\n" << endl;
  return(0);
}
