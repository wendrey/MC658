//----------------------------------------------------------------------
// Example of program that uses the LEMON Library, Linear Programming
//
// To use this program, you must install the LEMON Library and a linear programming
// package used by LEMON (like glpk, clp, cplex,...)
//
// This program find a minimum Traveling Salesman Tour (TSP) via a branch
// and cut approach. This program is for oriented graphs.
//
// This program is slower than "ex_tsp_undirected.cpp" (for undirected graphs)
// as separation is done testing all pair of vertices (we do not use gomory-hu tree
// in this case).
//
// To compile this program, read the README file
//
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------

#include "branch_cut.h"
#include "readgraph.h"
#include "viewgraph.h"
#include <lemon/list_graph.h>
#include <lemon/concepts/digraph.h>


#include <lemon/preflow.h>
//#include "branch_cut.cpp"
//#include "readgraph.cpp"
//#include "viewgraph.cpp"

typedef ListDigraph::Node Node;
typedef ListDigraph::NodeIt NodeIt;
typedef ListDigraph::Arc Arc;
typedef ListDigraph::ArcMap<double> ArcWeight;
typedef ListDigraph::NodeMap<string> NodeName;
typedef ListDigraph::ArcMap<string> ArcName;
typedef ListDigraph::NodeMap<int> NodeColor;
typedef ListDigraph::ArcMap<int> ArcColor;
typedef ListDigraph::NodeMap<bool> CutMap;
typedef Preflow<ListDigraph, ArcWeight> PType;
typedef ListDigraph::ArcMap<Lp::Col> EdgeLp;



typedef struct {
  ListDigraph *g;
  EdgeLp *x;
  NodeName *vname;
  ArcName *ename;
  NodeColor *vcolor;
  ArcColor *ecolor;
} BC_TransferType;


void ViewFracDigraph(ListDigraph &g,
	      ListDigraph::NodeMap<string> &vname,  // label of vertices
	      ListDigraph::ArcMap<string> &ename,   // label of arcs
	      ListDigraph::ArcMap<double> &weight,  // weight of arcs
	      ListDigraph::NodeMap<int> &vcolor,    // color of vertices
	      ListDigraph::ArcMap<int> &ecolor)     // color of arcs
{
    for (ListDigraph::ArcIt a(g); a!=INVALID; ++a) {
      if (fabs(weight[a])<0.001) ecolor[a]=WHITE;
      else if (BC_IsFrac(weight[a])) ecolor[a]=RED;
      else ecolor[a]=BLUE;
    }
    ViewDigraph(g,VIEW_NEATO,vname,ename,vcolor,ecolor);
}



string DoubleToString2(double x)
{
  std::ostringstream oss;
  oss << x;
  return(oss.str());
}

double MinCut(ListDigraph &g, ArcWeight &weight, Node &s,Node &t, CutMap &vcut)
{
  PType preflow_test(g, weight, s, t); 
  preflow_test.run();
  preflow_test.minCutMap(vcut);
  return (preflow_test.flowValue());
}



bool TSP_InsertGlobalCuts(BCTree &T)
{
  BC_TransferType *transfer;
  ListDigraph *g;
  EdgeLp *x;
  double vcut,vlastcut;
  bool foundcut;
  transfer = (BC_TransferType *) BC_GetProblemData(T);
 
  if (transfer==NULL) {
    cout << "TSP_InsertGlobalCuts: Error - NULL problem Data address." << endl;
    exit(0);
  }
  g = transfer->g;  // obtain the pointer for the graph
  x = transfer->x;  // obtain the pointer for the variables of the lp
  ArcWeight lpedgevalue(*g);
  CutMap cut(*g);
  for (ListDigraph::ArcIt a(*g); a!=INVALID; ++a) 
    // obtain the value of the lp variable in the current LP
    lpedgevalue[a] = BC_GetLpVariableValue(T,(*x)[a]);

  // For each pair of vertices, obtain a min cut and see if the cut is violated
  // It is not optimized, as you obtain the same cut many times
  foundcut = false;
  vlastcut = BC_INF;
  for (ListDigraph::NodeIt u(*g); u!=INVALID; ++u) {
    for (ListDigraph::NodeIt v(*g);v!=INVALID; ++v) {
      if (u==v) continue;
      vcut = MinCut(*g,lpedgevalue, u , v, cut);
      if ((vcut < 1.0-0.0001)&&(vcut<vlastcut-BC_EPS)) {
	Lp::Expr e;
	vlastcut = vcut;
	foundcut = true;
	for (ListDigraph::ArcIt a(*g); a!=INVALID; ++a) 
	  if (cut[(*g).source(a)] && !cut[(*g).target(a)]) e += (*x)[a];
	BC_AddConstraint(T, e >= 1 );
      }
    }
  }
  return(foundcut);
}

int main()
{

  ListDigraph g;
  ArcWeight weight(g);
  ArcWeight lpvariablevalue(g);
  NodeName vname(g);
  ArcName ename(g);
  ArcName enametoview(g);
  ArcColor ecolor(g);
  NodeColor vcolor(g);

  int nvar,i;
  BC_TransferType Transfer;
  
  ListDigraph::ArcMap<Lp::Col> x(g);
  
  ReadListDigraph(g,vname,weight,true,"gr_50"); 
  // gr_50 tem custo de 20.2833 in 41s
  // gr_100 tem custo de 2078 in 2m42s  no não orientado leva 42s.

  nvar = countArcs(g); // obtain number of arcs.
  BCTree T(nvar,BC_MINIMIZE);
  BC_SetPrintIterations(T,true);
  Transfer.g = &g;
  Transfer.x = &x;
  Transfer.vname = &vname;
  Transfer.ename = &ename;
  Transfer.vcolor = &vcolor;
  Transfer.ecolor = &ecolor;
  
  BC_SetProblemData(T,(void *) &Transfer);
  BC_SetNodeSeparationAlgorithm(T,TSP_InsertGlobalCuts);

  // Define one variable for each edge
  for (ListDigraph::ArcIt a(g); a!=INVALID; ++a,i++) {
    ecolor[a] = BLUE;
    x[a] = BC_GetNewVariable(T);
    ename[a] = DoubleToString2(weight[a]);
    BC_SetVarBounds(T , x[a] , 0.0 , 1.0); 
  }
  
  {// Define objective function
    Lp::Expr e;
    for (ListDigraph::ArcIt a(g); a!=INVALID; ++a) e += weight[a]*x[a];
    BC_ObjectiveFunction(T, e ); 
  }

  // Inserindo as restrições
  for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) {
    Lp::Expr e,f;
    vcolor[v] = GREY;
    // Constraint: for each vertex, exactly one arc must leave
    for (ListDigraph::OutArcIt a(g, v); a != INVALID; ++a) e += x[a];
    BC_AddConstraint(T, e == 1 );

    // Constraint: for each vertex, exactly one arc must enter 
    for (ListDigraph::InArcIt a(g, v); a != INVALID; ++a) f += x[a];
    BC_AddConstraint(T, f == 1 );
  }

  // Solving the system
  if (BC_Solve(T)==BC_OPTIMUM) {
    cout << "Value of the objective function: " << BC_GetSolutionValue(T) << endl;
    for (ListDigraph::ArcIt a(g); a!=INVALID; ++a) {
      lpvariablevalue[a] = BC_GetSolutionVariableValue(T,x[a]);
      enametoview[a]=ename[a]+"|"+DoubleToString2(lpvariablevalue[a]);
    }
    ViewFracDigraph(g,vname,enametoview,lpvariablevalue,vcolor,ecolor);
  } else {cout << "Could not obtain optimum solution." << endl;}
  return 0;
}
