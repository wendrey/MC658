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
#include <assert.h>
#include <math.h>
#include <sstream>
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
typedef ListGraph::NodeMap<double> NodePos;
typedef ListGraph::EdgeMap<string> EdgeName;
typedef ListGraph::NodeMap<bool> CutMap;
typedef Preflow<ListGraph, EdgeWeight> PFType;
typedef ListGraph::EdgeMap<Lp::Col> EdgeLp;
typedef ListGraph::NodeMap<int> NodeColor;
typedef ListGraph::EdgeMap<int> EdgeColor;

// This is the type used to obtain the pointer to the problem data. This pointer
// is stored in the branch and cut tree. And when we define separation routines,
// we can recover the pointer and access the problem data again.
class TSP_Data {
public:
  TSP_Data(ListGraph &graph,NodeName &nodename,NodePos &posicaox,
	   NodePos &posy,EdgeWeight &eweight,EdgeLp &xx);
  ListGraph &g;
  int NNodes,NEdges;
  int max_perturb2opt_it; // maximum number of iterations for heuristic Perturb2OPT
  int max_diving_it; // maximum number of iterations for heuristic Diving
  EdgeLp &x;
  NodeName &vname;
  EdgeName ename;
  NodeColor vcolor;
  EdgeColor ecolor;
  EdgeWeight &weight;
  NodePos &posx;
  NodePos &posy;
  AdjacencyMatrix AdjMat; // adjacency matrix
  EdgeIndex edge2index;    // for each edge, gives the index of lp variable
  vector<Edge> index2edge; // for each lp variable, gives the edge
  vector<Node> BestCircuit; // vector containing the best circuit found
  double BestCircuitValue;
};

TSP_Data::TSP_Data(ListGraph &graph,NodeName &nodename,NodePos &posicaox,
		   NodePos &posicaoy,EdgeWeight &eweight,EdgeLp &xx):
  g(graph),
  x(xx),
  vname(nodename),
  ename(graph),
  vcolor(graph),
  ecolor(graph),
  weight(eweight),
  posx(posicaox),
  posy(posicaoy),
  AdjMat(graph,eweight,100000000.0),
  edge2index(graph),
  index2edge(countEdges(graph)),
  BestCircuit(countEdges(graph)) 
{
  NNodes=countNodes(this->g);
  NEdges=countEdges(this->g);
  BestCircuitValue = DBL_MAX;
  max_perturb2opt_it = 3000; // default value
  max_diving_it = 5; // default value
}



bool coin(double x)
{
  if (((double) rand()/(double) RAND_MAX) < x) return(true);
  return(false);
}

bool Heuristic_2_OPT(AdjacencyMatrix &A,vector<Node> &Circuit,double &BestCircuitValue, int &NNodesCircuit)
{
  double CurrentWeight=0.0,Remove,Insert;
  bool globalimproved,improved;
  vector<Node> CircuitAux(NNodesCircuit);
  int i,j,k,l;
  for (int i=0;i<NNodesCircuit-1;i++) 
    CurrentWeight += A.Cost(Circuit[i],Circuit[i+1]);
  CurrentWeight += A.Cost(Circuit[NNodesCircuit-1],Circuit[0]);

  i=0;
  globalimproved = false;
  do {
    int i1,i2,j1,j2;
    improved=false;
    i1 = i;
    do {
      // one edge is (i1 , i1+1)
      i2 = (i1+1)%NNodesCircuit;
      for (j=2;j<NNodesCircuit-1;j++) {
	j1 = (i1+j)%NNodesCircuit;
	j2 = (j1+1)%NNodesCircuit;
	Remove = A.Cost(Circuit[i1],Circuit[i2])+A.Cost(Circuit[j1],Circuit[j2]);
	Insert = A.Cost(Circuit[i1],Circuit[j1])+A.Cost(Circuit[i2],Circuit[j2]);
	if (Remove-Insert > 0) {
	  k = 0;
	  CircuitAux[k++] = Circuit[i1];
	  for (l=j1;l!=i1;l=(l-1+NNodesCircuit)%NNodesCircuit)
	    CircuitAux[k++] = Circuit[l];
	  for (l=j2;l!=i1;l=(l+1)%NNodesCircuit)
	    CircuitAux[k++] = Circuit[l];
	  //assert(k!=NNodesCircuit);

	  for (k=0;k<NNodesCircuit;k++) Circuit[k] = CircuitAux[k];
	  CurrentWeight = CurrentWeight - Remove + Insert;
	  if (CurrentWeight < BestCircuitValue-BC_EPS) {
	    BestCircuitValue = CurrentWeight;
	    globalimproved = true;
	  }
	  improved = true; 
	}
      }
      i1=(i1+1)%NNodesCircuit;
    } while (i1!=i);
    i = (i+1)%NNodesCircuit;
  }while (improved);
  if (globalimproved) cout << "[Heuristic: 2OPT] New solution of value " << BestCircuitValue << "\n";
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

      assert((Adj1[u]==INVALID) || (Adj2[u]==INVALID));
      if (Adj1[u]==INVALID) Adj1[u] = v; else Adj2[u] = v;

      assert((Adj1[v]==INVALID) || (Adj2[v]==INVALID));
      if (Adj1[v]==INVALID) Adj1[v] = u; else Adj2[v] = u;
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
    assert((Adj1[u]==ant) || (Adj2[u]==ant));
    if (Adj1[u]==ant) {ant = u; u=Adj2[u]; } else { ant = u; u=Adj1[u]; }
    assert(u!=INVALID);
    i++;
  }
  assert(i==NNodesCircuit);
  return(true);
}


// This routine must be called when the vector x (indexed on the edges) is integer.
// The contained circuit is transformed into a circuit represented by a sequence of nodes.
bool Update_Circuit(BCTree &T)
{
  TSP_Data *tsp = (TSP_Data *) BC_GetProblemData(T);
  ListGraph::NodeMap<Node> Adj1(tsp->g), Adj2(tsp->g);
  double CircuitValue,f;
  Node FirstNode=INVALID,u,ant;
  int n,i,NNodesCircuit;

  NNodesCircuit=0; // used to count nodes, only to verify correctness
  // Adj1[v] and Adj2[v] have the two adjacent nodes of v in the circuit
  for (NodeIt v(tsp->g); v != INVALID; ++v) { Adj1[v]=INVALID; Adj2[v]=INVALID; }
  CircuitValue = 0.0;
  for (EdgeIt e(tsp->g); e != INVALID; ++e) { 
    f = T.lp.primal(T.lpvar[tsp->edge2index[e]]);
    if (BC_IsFrac(f)) return(false);
    if (f > 1.0-BC_EPS) {
      Node u,v;
      NNodesCircuit++;
      u = (tsp->g.u(e)); v = (tsp->g.v(e));
      CircuitValue += tsp->weight[e];
      assert((Adj1[u]==INVALID)||(Adj2[u]==INVALID));
      assert((Adj1[v]==INVALID)||(Adj2[v]==INVALID));
      if (Adj1[u]==INVALID) Adj1[u] = v; else Adj2[u] = v;
      if (Adj1[v]==INVALID) Adj1[v] = u; else Adj2[v] = u;
    }
  }
  if (CircuitValue > tsp->BestCircuitValue-BC_EPS) return(false);

  // First verify if the circuit is valid
  if (NNodesCircuit==0) return(false); // otherwise, the circuit is non null
  for (NodeIt v(tsp->g); v != INVALID; ++v) 
    assert((Adj1[v]!=INVALID) && (Adj2[v]!=INVALID));

  // find one node of the circuit
  for (NodeIt v(tsp->g); v != INVALID; ++v) if (Adj1[v]!=INVALID) {FirstNode = v; break;} 
  n = 1;     ant = FirstNode;     u = Adj1[FirstNode];
  // verify if each node have two adjacent nodes
  while (u!=FirstNode) {
    assert(u!=INVALID);
    if (Adj1[u]==ant) {ant = u; u=Adj2[u];} else {assert(Adj2[u]==ant); ant=u; u=Adj1[u]; }
    n++;
  }
  assert(n==tsp->NNodes);

  // the circuit is ok. So, update the best circuit
  tsp->BestCircuitValue = CircuitValue; 
  i = 1;
  tsp->BestCircuit[0] = FirstNode;
  u = Adj2[FirstNode];
  ant = FirstNode;
  while (u!=FirstNode) {
    tsp->BestCircuit[i] = u;      i++;
    if (Adj1[u]==ant) {ant=u; u=Adj2[u];} else {ant=u; u=Adj1[u];}
  }
  Heuristic_2_OPT(tsp->AdjMat,tsp->BestCircuit,tsp->BestCircuitValue,tsp->NNodes);
  assert(tsp->BestCircuitValue < T.best_ub);
  T.best_ub = tsp->BestCircuitValue;
  return(true);
}


bool TSP_InsertNodeCuts(BCTree &T);

bool TSP_DivingHeuristic(BCTree &T)
{
  TSP_Data *tsp;
  ListGraph *g;
  EdgeLp *x;
  double original_lb_var[T.num_var],original_ub_var[T.num_var],
    totfrac,valfracvar[T.num_var];
  int nfracvar,fracvar[T.num_var];
  tsp = (TSP_Data *) BC_GetProblemData(T);
  if (tsp->max_diving_it <1) return(false);
  g = &tsp->g;  x = &tsp->x;// obtain the graph g and the LP variables x
  
  for (int i=0;i<T.num_var;i++) {
    original_lb_var[i] = T.global_lb_var[i];
    original_ub_var[i] = T.global_ub_var[i];}

  for (int k=0;k<tsp->max_diving_it;k++) { // try max_diviing_it times to find integer solutions by diving
    cout << "[" << T.best_ub << "] Diving [" << k+1 << " of " << tsp->max_diving_it << "] ";

    while (1) { // repeats while there is frac. var. and is feasible/promising
      // while can insert cutting planes, solve the lp
      while (TSP_InsertNodeCuts(T)) T.lp.solve();
      if (T.lp.primalType()!=Lp::OPTIMAL) break; // lp is infeasible
      if (T.lp.primal() > T.best_ub-BC_EPS) break; // lp is feasible but not promisin

      // Choose a fractional variable
      EdgeWeight cap(*g);
      nfracvar = 0;
      totfrac = 0.0;
      for (int i=0;i<T.num_var;i++) {
	double f=T.lp.primal(T.lpvar[i]); // obtain the variable LP value 
	if (BC_IsFrac(f)) { // if it is frac, store in the vector of frac. vars.
	  valfracvar[nfracvar] = f;
	  fracvar[nfracvar] = i;
	  nfracvar++;
	  totfrac += (f*10)*(f*10);
	}
      }
      //cout << T.lp.primal() << "[" << nfracvar << "] ";
      fflush(stdout);
      if (nfracvar==0) { // obtained integer solution
	for (int i=0;i<T.num_var;i++) T.best_x[i] = T.lp.primal(T.lpvar[i]);
	T.best_ub = T.lp.primal();
	cout<<"[Heuristic: Diving] New solution of value "<<T.best_ub<<endl;fflush(stdout);
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

void ChangeNode(Node &a,Node &b)
{ Node c;  c = a; a = b; b = c; }
  

// This routine starts with some solution (current best or a random generated
// solution) and iteratively perturb changing some nodes and applying 2OPT.
bool TSP_Perturb2OPT(BCTree &T)
{
  TSP_Data *tsp;
  ListGraph *g;
  int nchanges,i,j;
  tsp = (TSP_Data *) BC_GetProblemData(T);
  g = &tsp->g;
  vector<Node> Circuit(tsp->NNodes), BestCircuit(tsp->NNodes);
  double BestCircuitValue;
  //nchanges = 2 (tests with instance a280 indicate that nchanges=1 is better
  nchanges = 1;

  // Start with a initial solution (if there is no solution, generate any sequence)
  if (tsp->BestCircuitValue < DBL_MAX) {
    BestCircuitValue = tsp->BestCircuitValue;
    for (int i=0; i<tsp->NNodes; i++) BestCircuit[i] = tsp->BestCircuit[i];
  } else {
    int i=0;
    BestCircuitValue = 0.0;
    for (ListGraph::NodeIt v(*g); v!=INVALID; ++v) BestCircuit[i++]=v;
    for (i=0;i<tsp->NNodes;i++) 
      BestCircuitValue += 
	tsp->AdjMat.Cost(BestCircuit[i] , BestCircuit[(i+1)%tsp->NNodes]);
  }

  for (int it=0;it<tsp->max_perturb2opt_it;it++) {
    if (!(it%100)) printf("[Heuristic: Perturbation+2OPT] it = %d (of %d)\n",it+1,tsp->max_perturb2opt_it);
    for (int k=0;k<tsp->NNodes;k++) Circuit[k] = BestCircuit[k];
    for (int nc=0;nc<nchanges;nc++) {
      i = (int) (drand48()*tsp->NNodes);
      j = (int) (drand48()*tsp->NNodes);
      if (i!=j) ChangeNode(Circuit[i],Circuit[j]);
    }
    Heuristic_2_OPT(tsp->AdjMat,Circuit,BestCircuitValue,tsp->NNodes);
    if (BestCircuitValue < T.best_ub) {
      tsp->BestCircuitValue = BestCircuitValue;
      T.best_ub = BestCircuitValue;
      for (int i=0;i<tsp->NNodes;i++) {
	BestCircuit[i] = Circuit[i];
	tsp->BestCircuit[i] = Circuit[i];
      }
    }
  }
  return(true);
}

bool TSP_RootHeuristics(BCTree &T)
{
  TSP_DivingHeuristic(T);
  TSP_Perturb2OPT(T);
  return(true);
}

// Convert a double value to string
string DoubleToString2(double x)
{ std::ostringstream oss; oss << x; return(oss.str()); }


void ViewTspCircuit(TSP_Data &tsp)
{
  ListGraph h;
  ListGraph::NodeMap<string> h_vname(h);  // node names
  ListGraph::NodeMap<Node> h_g2h(tsp.g);  // maps a node of g to a node of h
  ListGraph::NodeMap<double> h_posx(h);
  ListGraph::NodeMap<double> h_posy(h);
  ListGraph::NodeMap<int> vcolor(h);   // color of the vertices
  ListGraph::EdgeMap<int> acolor(h);  // color of edges
  ListGraph::EdgeMap<string> aname(h);  // name of edges
  for (ListGraph::NodeIt v(tsp.g); v!=INVALID; ++v) {
    Node hv;
    hv = h.addNode();
    h_g2h[v] = hv;
    h_posx[hv] = tsp.posx[v];
    h_posy[hv] = tsp.posy[v];
    h_vname[hv] = tsp.vname[v];
    vcolor[hv] = BLUE;
  }
  for (int i=0;i<tsp.NNodes;i++) {
    ListGraph::Node u,v;
    ListGraph::Edge a;
    u = tsp.BestCircuit[i]; 
    v = tsp.BestCircuit[(i+1) % tsp.NNodes]; 
    a = h.addEdge(h_g2h[u] , h_g2h[v]);
    aname[a] = "";
    acolor[a] = BLUE;
  }
  ViewGraph(h,VIEW_NEATO,h_vname,aname,h_posx,h_posy,vcolor,acolor);
}


// Routine that insert global cutting planes, for each branch and cut node.
bool TSP_InsertNodeCuts(BCTree &T)
{
  TSP_Data *tsp;
  ListGraph *g;
  EdgeLp *x;
  double vcut;
  bool foundcut;
  tsp = (TSP_Data *) BC_GetProblemData(T);
   if (tsp==NULL) {
    cout << "TSP_InsertNodeCuts: Error - NULL problem data address." << endl;
    exit(0); }
  // obtain the graph g, the linear programming variables x
  g = &tsp->g; x = &tsp->x;

  // obtain the value of the linear programming variables, for each edge
  EdgeWeight cap(*g);
  for (ListGraph::EdgeIt a(*g); a!=INVALID; ++a) 
    cap[a] = BC_GetLpVariableValue(T,(*x)[a]); // get the capacity of the edge

  // Obtain the Gomory-Hu tree using the LP variables for each edge as capacity.
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
  if (!foundcut) if (Update_Circuit(T)) cout<<"New solution of value "<<T.best_ub;
  return foundcut;
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
  int nvar;
  // Declare the graph variable and its data
  ListGraph g;
  EdgeLp x(g);
  EdgeWeight weight(g);
  EdgeWeight lpvar(g);
  NodeName vname(g);
  ListGraph::NodeMap<double> posx(g),posy(g);
  srand(1);

  if (argc!=2) {cout<<"Usage: "<< argv[0]<<" <graph_filename>"<<endl; exit(0);} 
  else if (!FileExists(argv[1])) {cout<<"File "<<argv[1]<<" does not exist."<<endl; exit(0);}

  // Read the graph
  if (!ReadGraph(g,vname,weight,posx,posy,argv[1])) 
    {cout<<"Error reading graph file "<<argv[1]<<"."<<endl;exit(0);}

  // round down the edge weight to integer
  // for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) weight[a] = round(weight[a]);

  // All the problem information in only one variable. The graph must have
  // already valid data. This will help when transfering information for the
  // Branch and Cut tree Tsp is a variable that has all information of
  // the problem. 
  TSP_Data tsp(g,vname,posx,posy,weight,x); 

  // Initialize the branch and cut tree with number of variables equal to number
  // of edges, and set the objective sense to minimization
  nvar = countEdges(g);
  BCTree T(nvar,BC_MINIMIZE);
  BC_SetProblemData(T,(void *) &tsp);
  
  tsp.max_perturb2opt_it=10000; // increase/reduce for +/- iterations of Perturb2OPT heuristic
  tsp.max_diving_it=0;   // increase/reduce this value for +/- iterations of Diving heuristic

  // variables that maps each edge of the graph to a linear programming variable
  // and vice-versa. Also, set the bounds of each variable to the interval [0,1]
  for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
    x[a] = BC_GetNewVariable(T);
    BC_SetVarBounds(T , x[a] , 0.0 , 1.0); 
    tsp.edge2index[a]= BC_GetLpVarIndex(T,x[a]);
    tsp.index2edge[BC_GetLpVarIndex(T,x[a])] = a;
  }

  BC_SetPrintIterations(T,true); // print some information in each iteration
  BC_SetNodeSeparationAlgorithm(T,TSP_InsertNodeCuts); // node separation routine (subtour elimination)
  BC_SetRootHeuristic(T,TSP_RootHeuristics); // TSP_RootHeuristics call diving and local search heuristics

  {// Define the objective function
    Lp::Expr e;
    for (ListGraph::EdgeIt a(g); a!=INVALID; ++a)  e += weight[a] * x[a];
    BC_ObjectiveFunction(T, e ); 
  }

  // Insert a constraint for each node v
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    Lp::Expr e;
    for (ListGraph::IncEdgeIt a(g, v); a != INVALID; ++a) e += x[a];
    BC_AddConstraint(T, e == 2.0 ); // Degree constraints
  }

  // Solve the linear program with branch and cut
  if (BC_Solve(T)==BC_OPTIMUM) {
    cout << "Value of the objective function: " << BC_GetSolutionValue(T) << endl;
    Update_Circuit(T); // Stores the LP integer solution in the circuit form (sequence of nodes)
    ViewTspCircuit(tsp); // display the solution (uses neato/graphviz and a pdf viewer) 
    cout << "Resolution time: ";  printtime(BC_GetSolverTime(T)); cout << endl;
  } else {cout << "No optimum solution was obtained." << endl;}
  return 0;
 }
 
