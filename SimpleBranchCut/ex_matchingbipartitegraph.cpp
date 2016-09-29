//----------------------------------------------------------------------
// Example of program that uses the LEMON Library, Linear Programming
//
// To use this program, you must install the LEMON Library and a linear programming
// package used by LEMON (like glpk, clp, cplex,...)
//
// Find a maximum matching in a bipartite graph
//
// Obs.: As the graph is bipartite, the matrix in the linear programming formulation
// is totally unimodular. So, the vertices of the defined polyhedron are integer
// and therefore, the solution of the linear program formulation is integer.
//
// The general formulation for this problem is given in page 53 in the slides
// http://www.ic.unicamp.br/~fkm/lectures/proglin.pdf
//
// To compile this program, read the README file
//
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------

#include<lemon/list_graph.h>
#include <lemon/preflow.h>
#include <lemon/lp.h>

using namespace std;
using namespace lemon;
typedef ListGraph Graph;
typedef Graph::EdgeMap<double> EdgeCost;
typedef Graph::Node Node;
typedef Graph::NodeIt NodeIt;
typedef Graph::EdgeIt EdgeIt;
typedef Graph::IncEdgeIt IncEdgeIt;
typedef Graph::Edge Edge;
int main()
{
  // declare a graph g
  Graph g;

  // declare weights for each edge of the graph
  EdgeCost weight(g);

  Lp lp; // Declare an LP variable (initially empty)

  // declare a linear programming variable for each edge of the graph g
  Graph::EdgeMap<Lp::Col> z(g);

  // The bipartite graph defined below is very particular but after the graph 
  // is defined, the formulation is valid for general bipartite graphs.
  // In this case, the graph has two sets of vertices: X, Y. 
  // The set X = {x1,...,x4} and the set Y={y1,...,y4}
  Node x1=g.addNode(), x2=g.addNode(), x3=g.addNode(), x4=g.addNode();
  Node y1=g.addNode(), y2=g.addNode(), y3=g.addNode(), y4=g.addNode();

  // Add the edges of the graph g and the weight for each edge
  Edge a;
  a=g.addEdge(x1,y1); weight[a]=10;
  a=g.addEdge(x1,y3); weight[a]=40;
  a=g.addEdge(x1,y2); weight[a]=30;
  a=g.addEdge(x2,y2); weight[a]=22;
  a=g.addEdge(x2,y3); weight[a]=23;
  a=g.addEdge(x2,y4); weight[a]=24;
  a=g.addEdge(x3,y4); weight[a]=30;
  a=g.addEdge(x3,y4); weight[a]=31;
  a=g.addEdge(x3,y4); weight[a]=32;
  a=g.addEdge(x4,y1); weight[a]=3;
  a=g.addEdge(x4,y2); weight[a]=2;
  a=g.addEdge(x4,y3); weight[a]=4;

  // Add a linear programming variable for each edge
  lp.addColSet(z);

  // add the constraint: Number of chosen edges in each vertex is at most 1
  for (NodeIt n(g); n != INVALID; ++n) {
    Lp::Expr e;
    for (IncEdgeIt a(g, n); a != INVALID; ++a) e += z[a];
    lp.addRow(e <= 1 );
  }

  // print the edges of the graph g
  cout << "Edges of the graph\n";
  for (EdgeIt a(g); a != INVALID; ++a) 
    cout << "(" << g.id(g.u(a)) << "," << g.id(g.v(a)) << "," << weight[a] << ")\n";

  // Build the objective function and set the lower and upper bound for each variable
  Lp::Expr o;
  for (EdgeIt a(g); a != INVALID; ++a) {
    o += weight[a]*z[a];
    lp.colLowerBound( z[a], 0.0 ); lp.colUpperBound( z[a], 1.0 ); 
  }
  lp.obj(o);

  // it is a maximization problem
  lp.max();

  // Solve the LP problem
  lp.solve();
  
  // Print the variables/solution
  if (lp.primalType() == Lp::OPTIMAL) {
    cout << "Value of the matching: " << lp.primal() << endl;
    cout << "Edges\n";
    for (EdgeIt a(g); a != INVALID; ++a) {
      if (lp.primal(z[a])>0.999)
	cout << "(" << g.id(g.u(a)) << "," << g.id(g.v(a)) << "," << weight[a] << ")\n";
    }


  }else if (lp.primalType() == Lp::INFEASIBLE) {
    cout << "The linear program is infeasible" << endl;
  }else if (lp.primalType() == Lp::UNBOUNDED) {
    cout << "The linear program is unbounded" << endl;
  }else {cout << "No optimum solution found." << endl;}
  return 0;
}
