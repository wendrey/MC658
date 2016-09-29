//----------------------------------------------------------------------
// Example of program that uses the LEMON Library, Linear Programming
//
// To use this program, you must install the LEMON Library and a linear programming
// package used by LEMON (like glpk, clp, cplex,...)
//
// This program find a minimum Traveling Salesman Tour (TSP) via a branch
// and cut approach. The formulation used is given in the page 100 of the
// slides in the link (in portuguese)
//
// http://www.ic.unicamp.br/~fkm/lectures/proglin.pdf
//
// The obtain the cuts that are violated, it use the gomory-hu subroutine, 
// available from LEMON package. For a short explanation why it is interesting
// to use Gomory-Hu tree, see page 142 of the above slides (in portuguese).
//
// To compile this program, read the README file
//
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------
#include <float.h>
#include "branch_cut.h"
#include "readgraph.h"
#include "viewgraph.h"
#include "adjacencymatrix.h"
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/gomory_hu.h>

#include "tempo.h"
#include <fstream>

//#define TSP_INF 100000000000

typedef ListGraph::Node Node;
typedef ListGraph::NodeIt NodeIt;
typedef ListGraph::Edge Edge;
typedef ListGraph::EdgeMap<double> EdgeWeight;
typedef ListGraph::EdgeMap<int> EdgeIndex;
typedef ListGraph::NodeMap<string> NodeName;
typedef ListGraph::EdgeMap<string> EdgeName;
typedef ListGraph::NodeMap<bool> CutMap;
typedef Preflow<ListGraph, EdgeWeight> PFType;
typedef ListGraph::EdgeMap<Lp::Col> EdgeLp;
typedef ListGraph::NodeMap<int> NodeColor;
typedef ListGraph::EdgeMap<int> EdgeColor;

void FreeNotNull(void *p){  if (p) free(p);  }


// This is the type used to obtain the pointer to the problem data. This pointer
// is stored in the branch and cut tree. And when we define separation routines,
// we can recover the pointer and access the problem data again.
class TSP_Data {
public:
  TSP_Data(ListGraph *graph);
  ListGraph *g;
  int NNodes,NEdges;
  EdgeLp x;
  NodeName vname;
  EdgeName ename;
  NodeColor vcolor;
  EdgeColor ecolor;
  EdgeWeight weight;
  AdjacencyMatrix AdjMat; // adjacency matrix
  EdgeIndex edge2index;    // for each edge, gives the index of lp variable
  vector<Edge> index2edge; // for each lp variable, gives the edge
  vector<Node> BestCircuit; // vector containing the best circuit found
  double BestCircuitValue;
};

// int NumberOfNodes(ListGraph &g)
// {
//   int n=0;
//   for (ListGraph::NodeIt u(g); u!=INVALID; ++u) n++;
//   return (n);
// }
// int NumberOfEdges(ListGraph &g)
// {
//   int n=0;
//   for (ListGraph::EdgeIt e(g); e!= INVALID; ++e) n++;
//   return (n);
// }

TSP_Data::TSP_Data(ListGraph *graph):
  x(*graph),
  vname(*graph),
  ename(*graph),
  vcolor(*graph),
  ecolor(*graph),
  weight(*graph),
  AdjMat(*graph,weight,-1.0),
  edge2index(*graph),
  index2edge(countNodes(*graph)),
  BestCircuit(countEdges(*graph)) 
  //index2edge(this->NNodes),
  //BestCircuit(this->NEdges) 
{
  g=graph;
  NNodes=countNodes(*graph);
  NEdges=countEdges(*graph);
  BestCircuitValue = DBL_MAX;
}



bool coin(double x)
{
  if (((double) rand()/(double) RAND_MAX) < x) return(true);
  return(false);
}

bool Heuristic_2_OPT(AdjacencyMatrix &A,vector<Node> &Circuit,int &NNodesCircuit)
{
  double CurrentWeight=0.0,Remove,Insert;
  bool globalimproved,improved;
  int i;
  for (int i=0;i<NNodesCircuit-1;i++) 
    CurrentWeight += A.Cost(Circuit[i],Circuit[i+1]);
  CurrentWeight += A.Cost(Circuit[NNodesCircuit-1],Circuit[0]);

  i=0;
  globalimproved = false;
  do {
    int i1,i2,j1,j2;
    improved=false;
    for (i1=i;i1!=i;i1=(i1+1)%NNodesCircuit) {
      i2 = (i1+1)%NNodesCircuit;
      for (j1=(i2+1)%NNodesCircuit, j2 = (j1+1)%NNodesCircuit; 
	   j2!=i1; j1=(j1+1)%NNodesCircuit, j2 = (j1+1)%NNodesCircuit) {
	Remove = A.Cost(Circuit[i1],Circuit[i2])+A.Cost(Circuit[j1],Circuit[j2]);
	Insert = A.Cost(Circuit[i1],Circuit[j1])+A.Cost(Circuit[i2],Circuit[j2]);
	if (Remove-Insert > 0) {
	  improved = true; globalimproved=true;
	}
      }
    }
    i = (i+1)%NNodesCircuit;
  }while (improved);
  return(globalimproved);
}

// Convert a circuit given by a binary vector of edges (edges belongs to
// the circuit if and only if the vector entry has value 1.0) to sequence 
// of vertices. So, a circuit with k vertices is represented by a sequence
// Circuit=(v0,v1,...,v_{k-1}), where k=NNodesCircuit.
bool Convert_Edge_Sequence(ListGraph &g,EdgeWeight &x,
		       vector<Node> &Circuit,int &NNodesCircuit)
{
  // given a binary vector x indexed on the edges, this routine produces 
  // vectors adj1 and adj2, indexed on the vertices, such that
  // adj1[v] and adj2[v] are the adjacent vertices of v.
  ListGraph::NodeMap<Node> Adj1(g);
  ListGraph::NodeMap<Node> Adj2(g);
  Node FirstNode=INVALID,u,ant;
  int i;

  NNodesCircuit=0;
  
  for (NodeIt v(g); v != INVALID; ++v) { Adj1[v]=INVALID; Adj2[v]=INVALID; }
  for (EdgeIt e(g); e != INVALID; ++e) { 
    if (x[e] > 1.0-BC_EPS) {
      Node u,v;
      NNodesCircuit++;
      u = (g.u(e)); v = (g.v(e));

      if (Adj1[u]==INVALID) Adj1[u] = v;
      else if (Adj2[u]==INVALID) Adj2[u] = v;
      else { cout << "Vector x in Convert_Edge_Sequence is not a circuit.1\n"; exit(0); }

      if (Adj1[v]==INVALID) Adj1[v] = u;
      else if (Adj2[v]==INVALID) Adj2[v] = u;
      else { cout << "Vector x in Convert_Edge_Sequence is not a circuit.2\n"; exit(0); }
    }
  }
  if (NNodesCircuit==0) return(false);
  for (NodeIt v(g); v != INVALID; ++v) 
    if (Adj1[v]!=INVALID) { FirstNode = v; break; }

  i = 1;
  Circuit[0] = FirstNode;
  ant = FirstNode;
  u = Adj1[FirstNode];
  while (u!=FirstNode) {
    Circuit[i] = u;
    // obtain the next vertex u
    if (Adj1[u]==ant) {ant = u; u=Adj2[u]; }
    else if (Adj2[u]==ant) { ant = u; u=Adj1[u]; }
    else { cout << "Vector x in Convert_Edge_Sequence is not a circuit.3\n"; exit(0); }
    if (u==INVALID) {cout << "Vector x in Convert_Edge_Sequence is not a circuit.4\n"; exit(0); }
    i++;
  }

  if (i!=NNodesCircuit) {cout << "Vector x in Convert_Edge_Sequence is not a circuit.5\n"; exit(0); }
  return(true);
}


// This routine must be called when the vector x (indexed on the edges) is integer.
// The contained circuit is transformed into a circuit represented by a sequence of nodes.
bool Update_Circuit(BCTree &T)
{
  TSP_Data *data;
  data = (TSP_Data *) BC_GetProblemData(T);
  ListGraph::NodeMap<Node> Adj1(*data->g);
  ListGraph::NodeMap<Node> Adj2(*data->g);
  double CircuitValue;
  Node FirstNode=INVALID,u,ant;
  int n,i,NNodesCircuit;

  NNodesCircuit=0; // used to count nodes, only to verify correctness
  // Adj1[v] and Adj2[v] have the two adjacent nodes of v in the circuit
  for (NodeIt v(*data->g); v != INVALID; ++v) { Adj1[v]=INVALID; Adj2[v]=INVALID; }
  CircuitValue = 0.0;
  for (EdgeIt e((*data->g)); e != INVALID; ++e) { 

    if (T.lp.primal(T.lpvar[data->edge2index[e]]) > 1.0-BC_EPS) {
      Node u,v;
      NNodesCircuit++;
      u = ((*data->g).u(e)); v = ((*data->g).v(e));
      CircuitValue += data->weight[e];

      if (Adj1[u]==INVALID) Adj1[u] = v;
      else if (Adj2[u]==INVALID) Adj2[u] = v;
      else { cout << "Vector x in Update_Circuit is not a circuit.1\n"; exit(0); }

      if (Adj1[v]==INVALID) Adj1[v] = u;
      else if (Adj2[v]==INVALID) Adj2[v] = u;
      else { cout << "Vector x in Update_Circuit is not a circuit.2\n"; exit(0); }
    }
  }
  
  if (CircuitValue > data->BestCircuitValue-BC_EPS) return(false);

  // First verify if the circuit is valid
  if (NNodesCircuit==0) return(false); // otherwise, the circuit is non null
  // find one node of the circuit
  for (NodeIt v(*data->g); v != INVALID; ++v) 
    if (Adj1[v]!=INVALID) {FirstNode = v; break;} // find one node of the circuit
  n = 1;
  ant = Adj1[FirstNode];
  u = Adj2[FirstNode];
  // count the number of nodes in the circuit
  while (u!=FirstNode) {
    if (u==INVALID) {cout << "Vector x in Update_Circuit is not a circuit.\n"; exit(0); }
    if (Adj1[u]==ant) {ant = u; u=Adj2[u]; }
    else if (Adj2[u]==ant) { ant = u; u=Adj1[u]; }
    else { cout << "Vector x in Update_Circuit is not a circuit.3\n"; exit(0); }
    n++;
  }
  if (n!= data->NNodes) return(false); // circuit does not have correct number of nodes
  
  // the circuit is ok. So, update the best circuit
  data->BestCircuitValue = CircuitValue; 
  i = 1;
  data->BestCircuit[0] = FirstNode;
  u = Adj2[FirstNode];
  ant = FirstNode;
  while (u!=FirstNode) {
    data->BestCircuit[i] = u;      i++;
    if (Adj1[u]==ant) { ant = u; u=Adj2[u]; }
    else              { ant = u; u=Adj1[u]; }
  }
  return(true);
}



bool TSP_InsertNodeCuts(BCTree &T);

bool TSP_DivingHeuristic1(BCTree &T)
{
  TSP_Data *problemdata;
  ListGraph *g;
  EdgeLp *x;
  double original_lb_var[T.num_var],original_ub_var[T.num_var];
  int imaxfrac,nfracvar,fracvar[T.num_var];
  double threshold,maxfrac,valfracvar[T.num_var];
  int maxit=7;

  problemdata = (TSP_Data *) BC_GetProblemData(T);
  
  for (int i=0;i<T.num_var;i++) {
    original_lb_var[i] = T.global_lb_var[i];
    original_ub_var[i] = T.global_ub_var[i];
  }

  // obtain the graph g, the linear programming variables x
  g = problemdata->g;
  x = &problemdata->x;

  threshold = 1.0;

  for (int k=0;k<maxit;k++) { // try 5 times to find integer solutions by diving
    threshold = threshold * 0.98;
    
    cout << "[" << T.best_ub << "] Rounding+Diving [" << k+1 << " of " << maxit << "] / rounding above " << threshold << ": ";

    while (1) { // repeats while there is frac. var. and is feasible/promising
	
      // while can insert cutting planes, solve the lp
      while (TSP_InsertNodeCuts(T)) T.lp.solve();
      if (T.lp.primalType()!=Lp::OPTIMAL) break; // lp is infeasible
      if (T.lp.primal() > T.best_ub-BC_EPS) break; // lp is feasible but not promising
      // Round 
      EdgeWeight cap(*g);
      nfracvar = 0;
      maxfrac = -10000000;
      imaxfrac = -1;
      for (int i=0;i<T.num_var;i++) 
	if (T.lp.primal(T.lpvar[i]) > threshold)  
	  if (coin(T.lp.primal(T.lpvar[i]))) T.lp.colLowerBound( T.lpvar[i], 1.0 );
      T.lp.solve();
      if (T.lp.primalType()!=Lp::OPTIMAL) break; // lp is infeasible
      if (T.lp.primal() > T.best_ub-BC_EPS) break; // lp is feasible but not promisin

      for (int i=0;i<T.num_var;i++) {
	// obtain the value of the variable in the current LP
	double f = T.lp.primal(T.lpvar[i]);
	if (BC_IsFrac(f)) {
	  valfracvar[nfracvar] = f;
	  fracvar[nfracvar] = i;
	  if (f>maxfrac) { maxfrac=f; imaxfrac=nfracvar; }
	  nfracvar++;
	}
      }
      cout << ".";
      fflush(stdout);
      if (nfracvar==0) { // obtained integer solution
	for (int i=0;i<T.num_var;i++) T.best_x[i] = T.lp.primal(T.lpvar[i]);
	T.best_ub = T.lp.primal();
	//cout << "\nNew solution: " << T.best_ub;
	Update_Circuit(T);
	break;
      }else{ // there is fractional variable
	int r=rand()%nfracvar; // choose one of the frac. var. randomly
	// and then round it to 0 or to 1
	if (maxfrac > 0.599) { r = imaxfrac;}
	
	if (valfracvar[r]>=0.5) T.lp.colLowerBound( T.lpvar[fracvar[r]], 1.0 ); 
	else T.lp.colUpperBound( T.lpvar[fracvar[r]], 0.0 ); 
	T.lp.solve();
      }
    }
    cout << endl;
    for (int i=0;i<T.num_var;i++) {
      T.lp.colLowerBound(T.lpvar[i],original_lb_var[i]);
      T.lp.colUpperBound(T.lpvar[i],original_ub_var[i]);
    }
    T.lp.solve();
  }
  return true;
}


bool TSP_DivingHeuristic2(BCTree &T)
{
  TSP_Data *problemdata;
  ListGraph *g;
  EdgeLp *x;
  double original_lb_var[T.num_var],original_ub_var[T.num_var];
  int nfracvar,fracvar[T.num_var];
  double totfrac,valfracvar[T.num_var];
  int maxit = 200;

  problemdata = (TSP_Data *) BC_GetProblemData(T);
  
  for (int i=0;i<T.num_var;i++) {
    original_lb_var[i] = T.global_lb_var[i];
    original_ub_var[i] = T.global_ub_var[i];
  }

  // obtain the graph g, the linear programming variables x
  g = problemdata->g;
  x = &problemdata->x;

  for (int k=0;k<maxit;k++) { // try 20 times to find integer solutions by diving
    cout << "[" << T.best_ub << "] Diving [" << k+1 << " of " << maxit << "]: ";

    while (1) { // repeats while there is frac. var. and is feasible/promising
	
      // while can insert cutting planes, solve the lp
      while (TSP_InsertNodeCuts(T)) T.lp.solve();
      if (T.lp.primalType()!=Lp::OPTIMAL) break; // lp is infeasible
      if (T.lp.primal() > T.best_ub-BC_EPS) break; // lp is feasible but not promisin

      // for (int i=0;i<T.num_var;i++) {
      // 	// obtain the value of the variable in the current LP
      // 	if (T.lp.primal(T.lpvar[i]) > 0.85)
      // 	  T.lp.colLowerBound(T.lpvar[i],1.0);
      // }
      // while (TSP_InsertNodeCuts(T)) T.lp.solve();
      // if (T.lp.primalType()!=Lp::OPTIMAL) break; // lp is infeasible
      // if (T.lp.primal() > T.best_ub-BC_EPS) break; // lp is feasible but not promisin
      
      
      // Now, let see if there is fractional variables
      // obtain the value of the linear programming variables, for each edge
      EdgeWeight cap(*g);
      nfracvar = 0;
      totfrac = 0.0;
      for (int i=0;i<T.num_var;i++) {
	// obtain the value of the variable in the current LP
	double f = T.lp.primal(T.lpvar[i]);
	if (BC_IsFrac(f)) {
	  valfracvar[nfracvar] = f;
	  fracvar[nfracvar] = i;
	  nfracvar++;
	  totfrac += (f*10)*(f*10);
	}
      }
      // cout << ".";
      cout << T.lp.primal() << "[" << nfracvar << "] ";
      fflush(stdout);
      if (nfracvar==0) { // obtained integer solution
	for (int i=0;i<T.num_var;i++) T.best_x[i] = T.lp.primal(T.lpvar[i]);
	T.best_ub = T.lp.primal();
	// cout << "\nNew solution: " << T.best_ub;
	Update_Circuit(T);
	break;
      }else{ // there is fractional variable
	int tmax = 2;
	for (int t=0;t<tmax;t++) {
	  int r=-1;
	  double bound;
	  double s=0.0, tr=totfrac*((double) rand())/((double) RAND_MAX);
	  for (int i=0;i<nfracvar;i++) {
	    s += (valfracvar[i]*10)*(valfracvar[i]*10);
	    if (s < tr+BC_EPS) r=i;
	    else break;
	  }
	  if (coin(valfracvar[r])) bound = 1.0;
	  else bound = 0.0;
	  T.lp.colLowerBound( T.lpvar[fracvar[r]], bound );
	  T.lp.colUpperBound( T.lpvar[fracvar[r]], bound ); 
	}
	T.lp.solve();
	if ((T.lp.primalType()!=Lp::OPTIMAL) || // lp is infeasible
	    (T.lp.primal() > T.best_ub-BC_EPS)) break;
      }
    }
    cout << endl;
    for (int i=0;i<T.num_var;i++) {
      T.lp.colLowerBound(T.lpvar[i],original_lb_var[i]);
      T.lp.colUpperBound(T.lpvar[i],original_ub_var[i]);
    }
    T.lp.solve();
  }
  return true;
}


bool TSP_RootHeuristics(BCTree &T)
{
  TSP_Data *problemdata;
  ListGraph *g;
  vector<Node> Circuit( problemdata->NNodes );
  int NNodesCircuit;
  problemdata = (TSP_Data *) BC_GetProblemData(T);
  g = problemdata->g;
  //Convert_Edge_Sequence(*g,T.best_x,Circuit,NNodesCircuit);

  TSP_DivingHeuristic1(T);
  TSP_DivingHeuristic2(T);
  Heuristic_2_OPT(problemdata->AdjMat,Circuit,NNodesCircuit);
  return(true);
}


// Routine to visualize a fractional graph.
void ViewFracGraph(ListGraph &g,
	      ListGraph::NodeMap<string> &vname,  // nome dos vertices
	      ListGraph::EdgeMap<double> &weight,  // peso das arestas
	      ListGraph::EdgeMap<string> &ename,   // nome das arestas
	      ListGraph::NodeMap<double>& posx,
	      ListGraph::NodeMap<double>& posy,
	      ListGraph::NodeMap<int> &vcolor,    // cor dos vertices
	      ListGraph::EdgeMap<int> &ecolor)     // cor das arestas
{
    for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
      if (fabs(weight[a])<0.001) ecolor[a]=WHITE;
      else if (BC_IsFrac(weight[a])) ecolor[a]=RED;
      else ecolor[a]=BLUE;
    }
    ViewGraph(g,VIEW_NEATO,vname,ename,posx,posy,vcolor,ecolor);
}


// Convert a double value to string
string DoubleToString2(double x)
{
  std::ostringstream oss;
  oss << x;
  return(oss.str());
}

// obtain a mincut separating vertices 's' and 't'
// The input is given by the graph 'g'. Each edge of g has a "weight".
// The returned cut is given by the vector of nodes 'cut' (boolean
// vector, nodes in one side have value false and nodes in the other side
// have value true).
double MinCut(ListGraph &g, EdgeWeight &weight, Node &s,Node &t, CutMap &cut)
{
  PFType pf(g, weight, s, t); 
  pf.runMinCut();
  pf.minCutMap(cut);
  return (pf.flowValue());
}

// Routine that insert global cutting planes, for each branch and cut node.
bool TSP_InsertNodeCuts(BCTree &T)
{
  TSP_Data *problemdata;
  ListGraph *g;
  EdgeLp *x;
  double vcut;
  bool foundcut;
  problemdata = (TSP_Data *) BC_GetProblemData(T);
 
  if (problemdata==NULL) {
    cout << "TSP_InsertNodeCuts: Error - NULL problem Data address." << endl;
    exit(0);
  }
  // obtain the graph g, the linear programming variables x
  g = problemdata->g;
  x = &problemdata->x;

  // obtain the value of the linear programming variables, for each edge
  EdgeWeight cap(*g);
  for (ListGraph::EdgeIt a(*g); a!=INVALID; ++a) 
    // obtain the value of the variable in the current LP
    cap[a] = BC_GetLpVariableValue(T,(*x)[a]);

  // Obtain the Gomory-Hu tree using the linear programming variables for each
  // edge as capacity.
  GomoryHu<ListGraph,EdgeWeight> ght((*g),cap);  
  ght.run();

  // The Gomory-Hu tree is given as a rooted directed tree. Each node has
  // an arc that points to its father. The root node has father -1.
  // Remember that each arc in this tree represents a cut and the value of
  // the arc is the weight of the corresponding cut. So, if an arc has weight
  // less than 2 than we found a violated cut and in this case, we insert the
  // corresponding constraint.
  foundcut = false;
  for (NodeIt u(*g); u != INVALID; ++u) {
    Lp::Expr e; 	// Insert constraint

    if ((*g).id(ght.predNode(u))==-1) continue;
    vcut =  ght.predValue(u);
    if (vcut > 2.0 - BC_EPS) continue;
    foundcut = true;
    for(GomoryHu<ListGraph,EdgeWeight>::MinCutEdgeIt a(ght,u,ght.predNode(u));a!=INVALID;++a) 
      e += (*x)[a];
    BC_AddConstraint(T, e >= 2.0 );
  }
  return foundcut;
}

// This is a routine that looks for a violated cut testing the mincut for //
// all pair of vertices. Clearly, the above version is much more efficient.
bool TSP_InsertNodeCutsDummyVersion(BCTree &T)
{
  TSP_Data *problemdata;
  ListGraph *g;
  EdgeLp *x;
  double vcut,vmincut;
  bool foundcut;
  problemdata = (TSP_Data *) BC_GetProblemData(T);
 
  if (problemdata==NULL) {
    cout << "TSP_InsertNodeCuts: Error - NULL problem Data address." << endl;
    exit(0);
  }
  g = problemdata->g;
  x = &problemdata->x;
  EdgeWeight lpedgevalue(*g);
  CutMap cut(*g);
  for (ListGraph::EdgeIt a(*g); a!=INVALID; ++a) 
    // obtain the value of the variable in the current LP
    lpedgevalue[a] = BC_GetLpVariableValue(T,(*x)[a]);

  foundcut = false;
  vmincut = 2.0;
  for (ListGraph::NodeIt u(*g); u!=INVALID; ++u) {
    for (ListGraph::NodeIt v(*g);v!=INVALID; ++v) {
      if (u==v) continue;
      vcut = MinCut(*g,lpedgevalue, u , v, cut);
      if ((vcut < 2.0-0.0001)&&(vcut<vmincut-0.0001)) { // Found violated cut
	vmincut = vcut; // next cut must have value smaller than the previous cut
	foundcut = true; // so, we guarantee that the new cut is different
	Lp::Expr e; 	// Insert constraint
	for (ListGraph::EdgeIt a(*g); a!=INVALID; ++a) 
	  if (cut[(*g).u(a)] != cut[(*g).v(a)]) e += (*x)[a];
	BC_AddConstraint(T, e >= 2 );
      }
    }
  }
  return(foundcut);
}

// Only to see if a file exists. It (tries to) open and close the file.
bool FileExists(const char *filename)
{
  bool r;
  ifstream ifile(filename);
  r = ifile;
  ifile.close();
  return r;
}


int main(int argc, char *argv[]) 
{

  // All the problem information in only one variable. This will help when
  // transfering information for the Branch and Cut tree

  ListGraph g;
  EdgeWeight weight(g);
  EdgeWeight lpvariablevalue(g);
  NodeName vname(g);
  EdgeName ename(g);
  EdgeName enametoview(g);
  EdgeColor ecolor(g);
  NodeColor vcolor(g);
    
  ListGraph::NodeMap<double> posx(g),posy(g);

  int nvar,Nnodes;
  
  ListGraph::EdgeMap<Lp::Col> x(g);

  srand(1);
 
  if (argc!=2) {
    cout << "Usage: " << argv[0] << " <graph_filename>" << endl;
    exit(0);
  } else if (!FileExists(argv[1])) {
    cout << "File " << argv[1] << " does not exist." << endl;
    exit(0);
  }
  if (!ReadGraph(g,vname,weight,posx,posy,argv[1])) {
    cout << "Error reading graph file " << argv[1] << "." << endl;
    exit(0);
  }
  Nnodes = countNodes(g);

  
  vector<Node> BestCircuit(Nnodes);
  double BestCircuitValue;

  nvar=0;// round down the edge weight to integer, and count the number of edges
  for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
    weight[a] = round(weight[a]);
    nvar++;
  }
  EdgeIndex edge2index(g);
  vector<Edge> index2edge(nvar);

  // obtain an adjacency matrix of A
  AdjacencyMatrix A(g,weight,-1.0);

  // Initialize the branch and cut tree with nvar variables and set the
  //  objective sense to minimization
  BCTree T(nvar,BC_MINIMIZE);

  BC_SetPrintIterations(T,true);

  // Set the node separation routine as the routine TSP_InsertNodeCuts
  BC_SetNodeSeparationAlgorithm(T,TSP_InsertNodeCuts);

  // Simple dive heuristics. For large graphs these heuristics are too slow
  // and does not give good results
  BC_SetRootHeuristic(T,TSP_RootHeuristics);

  // Obtain a linear programming variable for each edge and define the lower and
  // upper bound for each variable.
  for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
    ecolor[a] = BLUE;
    x[a] = BC_GetNewVariable(T);
    ename[a] = DoubleToString2(weight[a]);
    BC_SetVarBounds(T , x[a] , 0.0 , 1.0); 
    edge2index[a]= BC_GetLpVarIndex(T,x[a]);
    index2edge[BC_GetLpVarIndex(T,x[a])] = a;
  }

  // ProblemData is a variable that has all information of the problem. Note that 
  // it has only pointers for the original data.
  // The ProblemData variable will be used by separation routines, that will obtain
  // the problem data (this will be a pointer to ProblemData) from the branch and cut
  // tree. For example, see the separation routines TSP_InsertNodeCuts.
  TSP_Data ProblemData(&g);
  //ProblemData.g = &g;
  ProblemData.x = x;
  ProblemData.NNodes =   countNodes(g);
  ProblemData.NEdges =   countEdges(g);
  ProblemData.vname = &vname;
  ProblemData.ename = &ename;
  ProblemData.vcolor = &vcolor;
  ProblemData.ecolor = &ecolor;
  ProblemData.AdjMat = &A;
  ProblemData.index2edge = &index2edge;
  ProblemData.edge2index = &edge2index;
  ProblemData.BestCircuit = &BestCircuit;
  ProblemData.BestCircuitValue = &BestCircuitValue;
  BC_SetProblemData(T,(void *) &ProblemData);

  
  {// Define the objective function
    Lp::Expr e;
    for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) e += weight[a]*x[a];
    BC_ObjectiveFunction(T, e ); 
  }

  // Insert a constraint for each node v
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    Lp::Expr e;
    vcolor[v] = GREY;

    // Insert constrant: In the sum of edge varibles incident to a vertex 
    // is equal to 2
    for (ListGraph::IncEdgeIt a(g, v); a != INVALID; ++a) e += x[a];
    BC_AddConstraint(T, e == 2.0 );
  }

  // Solve the linear programming with branch and cut
  if (BC_Solve(T)==BC_OPTIMUM) {
    double soma=0.0;
    cout << "Value of the objective function: " << BC_GetSolutionValue(T) << endl;

    // Paint the edges of the solution with RED. 
    // Edges not in the solution are not displayed
    for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
      lpvariablevalue[a] = BC_GetSolutionVariableValue(T,x[a]);
      if (lpvariablevalue[a]>0.9999) {
	ecolor[a] = RED;
	soma += A.Cost(a);
      } else ecolor[a] = WHITE; // white edges are not displayed
    }

    { // The code below is only to check if the solution is correct.
      vector<Node> Circuit(A.Nnodes);
      int NNodesCircuit;
      soma=0.0;
      Convert_Edge_Sequence(g,lpvariablevalue,Circuit,NNodesCircuit);
      for (int i=0;i<NNodesCircuit;i++) 
	soma += A.Cost(Circuit[i],Circuit[(i+1)%NNodesCircuit]);
    }

    ViewFracGraph(g,vname,lpvariablevalue,enametoview,posx,posy,vcolor,ecolor);
    cout << "Conference of solution cost = "<< soma << endl;
    cout << "Resolution time: ";
    printtime(BC_GetSolverTime(T));
    cout << endl;
    //ViewEuclideanGraph(g,vname,posx,posy,vcolor,ecolor);
  } else {cout << "No optimum solution was obtained." << endl;}
  return 0;
 }
 
