
//----------------------------------------------------------------------
// Example of program that uses the LEMON Library and Linear Programming
// 
// Implementation of a simple Branch & Cut for didactic purpouses.
//
// Obs.: A high level description of this branch and cut (in portuguese) is
// available in http://www.ic.unicamp.br/~fkm/lectures/proglin.pdf
//
// This implementation does not use: tailing off, pool of inequalities,
// elimination of variables by reduced cost, or other optimization strategies.
// Any cut (linear constraint) that is inserted is global. 
// The branching process is made on variables. All variables must be integer.
// It works as a maximization or minimization problem, but it is implemented
// as a minimization problem. When it is used for a maximization problem the
// objective funciton is multiplied by (-1). The branch and cut tree is 
// represented by a heap (data structure available in the standard library of c++).
//
// The idea is to maintain the code very simple for didactic purpouses.
//
// To use this program, must install the LEMON Library and a linear programming
// package used by LEMON (like glpk, clp, cplex,...)
//
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------


#include <cstdio>
#include "branch_cut.h"
#include <time.h>
#include "tempo.h"

// Function used in the priority queue.
// Nodes with small lower bound are removed first.
bool CompareBCNode::operator()(BCNode &a, BCNode &b)
{  return(a.lb_value > b.lb_value);  }
//{  return(a.nvarfrac > b.nvarfrac);  }

// Node constructor
BCNode::BCNode(int maxvar)
{
  lb_var.resize(maxvar);
  ub_var.resize(maxvar);
}

// Branch and Cut Tree Constructor
BCTree::BCTree(int maxvar,int sense)
{
  printiterations=true;
  iterations=0;
  num_constraints=0;
  num_bbnodes=0;
  lpvar.resize(maxvar);
  global_lb_var.resize(maxvar);
  global_ub_var.resize(maxvar);
  best_x.resize(maxvar);
  best_lb = BC_MINUS_INF;
  best_ub = BC_INF;
  node_separation_algorithm = NULL; // initially no separation algorithm for nodes (including the root node)
  root_separation_algorithm = NULL; // initially no separation algorithm on the root node
  node_heuristic = NULL; // initially no heuristic in the nodes (including the root node)
  root_heuristic = NULL; // initially no heuristic in the root
  problem_data = NULL; // initially no problem data associated
  max_var = maxvar; // maximum number of variables
  num_var = 0; // number of current variables

  if ((sense!=BC_MINIMIZE) && (sense!=BC_MAXIMIZE)){
    printf("Objective function must be BC_MINIMIZE (for minimization) "
	   "or BC_MAXIMIZE (for maximization).\n");
    exit(0);
  }
  obj_sense = sense;
  for (int i=0;i<max_var;i++) {
    global_lb_var[i] = BC_MINUS_INF;
    global_ub_var[i] = BC_INF;
    lpvar[i] = lp.addCol();  // para cada i, associa a variável no lp
    if (lp.id(lpvar[i])!=i) {
      printf("BCTree constructor: Error - found situation "
	     "where new variables have non-sequential lp.id().\n");
      exit(0);
    }
  }
  initial_time = 0;
  ending_time = 0;
}


// Set the lower bound and the upper bound of a variable
void BC_SetVarBounds(BCTree &T,const Lp::Col var,double lb,double ub)
{
  T.lp.colLowerBound( var, lb ); 
  T.lp.colUpperBound( var, ub ); 
  T.global_lb_var[T.lp.id(var)] = lb;
  T.global_ub_var[T.lp.id(var)] = ub;
}

// Add a global linear constraint to the formulation
void BC_AddConstraint(BCTree &T,const Lp::Constr &C)
{  T.num_constraints++; T.lp.addRow( C );  }


// Set the objective function of the problem. It can be BC_MINIMIZE or BC_MAXIMIZE.
// Note that all branch and cut program is for minimization. So, if the user problem
// is of maximization, the objective function is multiplied by (-1).
void BC_ObjectiveFunction(BCTree &T,const Lp::Expr &E)
{ 
  if (T.obj_sense==BC_MINIMIZE) { // Minimization problem
     T.lp.obj( E );
  } else { // Maximization problem is transformed into a minimization problem
    Lp::Expr E2;
    E2 = -1.0*(E);
    T.lp.obj( E2 );
  }
  T.lp.min();
}

// Return true if variable x is fractional (within a certain small error).
bool BC_IsFrac(double x)
{
  double f;
  f = ceil(x)-x;
  if (f<BC_EPS) return(false);
  if (f>1.0-BC_EPS) return(false);
  return(true);
}

// obtain the fractional part of a number and returns the distance from the integer
double BC_FracLevel(double x) 
{
  double f;
  f = x-floor(x); // fractional part
  if ((f<BC_EPS) || (f>1.0-BC_EPS)) return(0.0);
  return(0.5-fabs(0.5 - f)); 
}

// return the index of the 'most' fractional variable. This routine is not
// used. It was replaced by the function BC_FracVar below. Probably, this routine
// can substitute the other when you have a good initial solutions (say obtained
// by heuristics) and the task is mainly to prove optimality.
int BC_FracVar2(BCTree &T)
{
  int imax=-1;
  double vmax=BC_EPS,f;
  for (int i=0;i<T.num_var;i++){
    f = BC_FracLevel(T.lp.primal(T.lpvar[i]));
    if (f > vmax) {   imax = i;    vmax = f;    }
  }
  return(imax);
}


// return the index of a fractional variable
// in 1/10 of the cases (see below) it returns the most fractional variable
// in 9/10 of the cases it returns the less fractional variable (not integer).
// See also the comments of the above routine.
int BC_FracVar(BCTree &T,int &nvarfrac)
{
  int imax=-1,imin=-1,nfrac=0;
  double vmax=BC_EPS,vmin=1.0,f;
  for (int i=0;i<T.num_var;i++){
    f = BC_FracLevel(T.lp.primal(T.lpvar[i]));
    if (f > vmax) {   imax = i;    vmax = f;    }
    if ((f > 0) && (f<vmin))  {   imin = i;    vmin = f;    }
    if (BC_IsFrac(f)) nfrac++;
  }
  nvarfrac=nfrac;
  if (rand()%10 == 0) return(imax);
  return(imin);
  // testing "rand()%10" in "ex_tsp_undirected.e gr_200" solved in 2min11seg. 
  // testing "rand()%5" in "ex_tsp_undirected.e gr_200" solved in 6min. 
  // testing "rand()%2" in "ex_tsp_undirected.e gr_200" solved in 7min. 
  // you must test the best for your problem.
}


void BC_PrintDouble(double x)
{
  if (x==BC_MINUS_INF) printf("-INF  ");
  else if (x==BC_INF) printf("INF  ");
  else printf("%6.2lf",x);
}


// To solve a node, set the lower and upper bound of the variables of 
// the linear program and then solve the linear program. If all variables are
// integer, it obtained a feasible solution and it will be compared with the 
// best upper bound found before.
bool BC_SolveNode(BCNode &No)
{ int nvar;
  BCTree *T;
  T = No.T;
  nvar=T->num_var;
  for (int i=0;i<nvar;i++) {
    T->lp.colLowerBound(T->lpvar[i],No.lb_var[i]);
    T->lp.colUpperBound(T->lpvar[i],No.ub_var[i]);
  }
  T->lp.solve();
  // If there is cutting plane algorithm, call it until cannot separate
  if (T->node_separation_algorithm) 
    while ((*(T->node_separation_algorithm))(*T)) T->lp.solve();
  // run the node heuristic, if it exists
  if (T->node_heuristic) (*(T->node_heuristic))(*T);
  if (T->lp.primalType() != Lp::OPTIMAL) {
    if (T->lp.primalType() == Lp::INFEASIBLE)     No.lp_status = BC_INFEASIBLE;
    else if (T->lp.primalType() == Lp::UNBOUNDED) No.lp_status = BC_UNBOUNDED;
    else if (T->lp.primalType() == Lp::UNDEFINED) No.lp_status = BC_UNDEFINED;
    else if (T->lp.primalType() == Lp::FEASIBLE)  No.lp_status = BC_UNDEFINED;
    return(false);
  }
  No.lb_value = T->lp.primal();
  if (No.lb_value >T->best_ub-BC_EPS) return(false); // the node will be pruned

  No.ind_frac = BC_FracVar(*T,No.nvarfrac);
  if (No.ind_frac==-1) {// no fractional variable (found integer solution)
    No.lp_status = BC_INTEGER;
    if (No.lb_value < T->best_ub) { // compare with the previous best integer sol.
      T->best_ub = No.lb_value; // the best integer sol. is updated.
      for (int i=0;i<nvar;i++) T->best_x[i] = T->lp.primal(T->lpvar[i]);
    }
    return(false); // the node will be pruned
  }
  No.lp_status = BC_FEASIBLE;
  No.val_frac = T->lp.primal(T->lpvar[No.ind_frac]);
  return(true);
}


// print some informations in each iteration
void BC_PrintIterations(BCTree &T)
{
  printf("It[%2lu] ",T.iterations);
  if (T.obj_sense==BC_MINIMIZE) {
    printf("Bnds["); 
    BC_PrintDouble(T.best_lb);        printf(",");  BC_PrintDouble(T.best_ub);
    printf("]");
  }else{
    printf("Bnds:["); 
    BC_PrintDouble((-1.0)*T.best_ub); printf(",");  BC_PrintDouble((-1.0)*T.best_lb);
    printf(" ]");
  }
  printf(" Act[%2lu] TotNod[%2lu] Constr[%2lu] ",T.BCHeap.size(),T.num_bbnodes,T.num_constraints);
  printf("Time[");
  shortprinttime(time70()-T.initial_time);
  printf("]\n");
}

// Boolean function to determine if a stop condition is satisfied.
bool BC_SatisfyStopCondition(BCTree &T)
{
  if (T.BCHeap.empty()) 
    return(true); // no more branch and bound nodes in the heap
  if (T.best_ub - T.best_lb < BC_EPS) 
    return(true); // solution obtained is already very close to opt.
                  // if you know that all solutions have integer values, you
                  // may also use this information to stop when the difference 
                  // from the upper and lower bounds is less than 1.
  
  // put other conditions here: e.g., by time

  return(false);
}


// This is the main routine of the branch and cut process.
// It starts the initial information and already solve the root node.
// If the root node is infeasible, unbounded or undefined, all the 
// integer formulation is also infeasible, unbounded or undefined.
// Then the routine inserts the root node in the heap (set of active nodes).
// This heap only receives nodes that have fractional variables.
// The main loop of BC_Solve routine remove a node, looks for a fractional variable,
// say x_i with value F. Then, it splits the node in two nodes changing the bounds
// of x_i. In one node it insert the bound x_i<=floor(F) and in the other node it
// insert the bound x_i >= ceil(F).
//
BC_Status BC_Solve(BCTree &T)
{
  BCNode No(T.max_var);
  double f,lb_var,ub_var,No_lower_bound,No_upper_bound;
  int ifrac;
  bool EmptyTree;
  srand(1);
  T.initial_time = time70();
  T.ending_time = time70();

  //  Define root node
  No.T = &T;
  No.ub_value = T.best_ub;  No.lb_value = T.best_lb;
  for (int i=0;i<T.num_var;i++) {
    No.lb_var[i] = T.global_lb_var[i]; No.ub_var[i] = T.global_ub_var[i]; }

  BC_SolveNode(No);
  if (No.lp_status==BC_INFEASIBLE) return(BC_INFEASIBLE);
  if (No.lp_status==BC_UNBOUNDED)  return(BC_UNBOUNDED);
  if (No.lp_status==BC_UNDEFINED) return(BC_UNDEFINED);

  // run the root separation heurist, if it exists
  if (T.root_separation_algorithm) 
    while ((*(T.root_separation_algorithm))(T)) T.lp.solve();

  // run the root heuristic, if it exists
  if (T.root_heuristic) (*(T.root_heuristic))(T);

  T.best_lb = T.lp.primal();
  if (No.lp_status==BC_INTEGER) { // Root node already integer
    No.lp_status = BC_OPTIMUM;
    return(No.lp_status); 
  }

  T.BCHeap.push(No); // Insert the root node in the heap (set of active nodes).

  T.num_bbnodes=1;
  while (!BC_SatisfyStopCondition(T)) {
    T.iterations++;  if (T.printiterations) BC_PrintIterations(T);
    
    // obtain (and remove) the next node from the heap
    No =  T.BCHeap.top();   T.BCHeap.pop();
    T.best_lb = No.lb_value; // Update global lower bound
    No_lower_bound = max(No.lb_value,T.best_lb); // Update node bounds
    No_upper_bound = min(No.ub_value,T.best_ub);
    ifrac = No.ind_frac;  f = No.val_frac;

    // Adding child node 1
    No.ub_value = No_upper_bound; No.lb_value = No_lower_bound;
    lb_var = No.lb_var[ifrac];  ub_var = No.ub_var[ifrac];
    No.ub_var[ifrac] = floor(f);  
    if (BC_SolveNode(No)) { T.BCHeap.push(No); T.num_bbnodes++; }
    
    // Adding child node 2
    No.lb_value = max(No_lower_bound,T.best_lb);
    No.ub_value = min(No_upper_bound,T.best_ub);
    No.ub_var[ifrac] = ub_var;
    No.lb_var[ifrac] = ceil(f);  
    if (BC_SolveNode(No)) { T.BCHeap.push(No);  T.num_bbnodes++; }    
  }
  EmptyTree = T.BCHeap.empty();

  T.ending_time = time70(); 
  while (!T.BCHeap.empty()) T.BCHeap.pop();  // Remove remaining nodes in the Heap
  if (T.printiterations) BC_PrintIterations(T);
  if (EmptyTree) { // performed complete enumeration
    if (T.best_ub < BC_INF-BC_EPS) return(BC_OPTIMUM);
    return(BC_INFEASIBLE);
  } else { // could not perform the complete enumeration
    if (T.best_ub - T.best_lb < BC_EPS) return(BC_OPTIMUM); //sol. considered optimum
    else if (T.best_ub < BC_INF-BC_EPS) return(BC_FEASIBLE);
    return(BC_UNDEFINED);
  }
}

// Return a new linear programming variable. Note that the maximum number of
// variables is defined when you declare the branch and cut tree T.
Lp::Col BC_GetNewVariable(BCTree &T)
{
  if (T.num_var < T.max_var) {
    T.num_var++;
    return(T.lpvar[T.num_var-1]);
  }else{
    printf("BC_GetNewVariable: "
	   "Error: maximum number of variables (of %d) attained.\n",T.max_var);
    exit(0);
  }
}

// Return the elapsed time since the starting time of the BC_Solve routine.
long BC_GetSolverTime(BCTree &T)
{  return(T.ending_time - T.initial_time); }

// Return the value of the linear programming variable.
double BC_GetSolutionVariableValue(BCTree &T,Lp::Col v)
{  return(T.best_x[T.lp.id(v)]);}

// return the value of the objective function. Note that for maximization problem
// the function was multiplied by (-1). So, the returned value must also be
// multiplied by (-1).
double BC_GetSolutionValue(BCTree &T)
{ if (T.obj_sense==BC_MINIMIZE) return(T.best_ub); else return((-1.0)*T.best_ub); }

// return the global lower bound of the problem. Note that for maximization problem
// the objective function was transformed into a minimization function. So, the 
// current lower bound is in fact an upper bound.
void BC_SetGlobalLowerBound(BCTree &T,double value)
{ if (T.obj_sense==BC_MINIMIZE) T.best_lb = value; else T.best_ub = value; }

// return the global upper bound of the problem. See note for the global 
// lower bound above.
void BC_SetGlobalUpperBound(BCTree &T,double value)
{ if (T.obj_sense==BC_MINIMIZE) T.best_ub = value; else T.best_lb = value; }

// Place the pointer of the problem in the branch and cut structure
// Possibly, if you use some heuristic based on the linear programming variables,
// you may need to access the original data of the problem.
void BC_SetProblemData(BCTree &T,void *ProblemData)
{ T.problem_data = ProblemData; }

// Get the pointer of the problem from the branch and cut structure
void *BC_GetProblemData(BCTree &T)
{ return(T.problem_data); }

// Set the routine that insert global cutting planes (linear constraints) in 
// the linear program, when solving a node (even if the node is the root node).
void BC_SetNodeSeparationAlgorithm(BCTree &T, bool (*SeparationAlgorithm)(BCTree &T))
{ T.node_separation_algorithm = SeparationAlgorithm; }

// Set the routine that insert global cutting planes only in the root node.
// In some problems, it is interesting to insert cutting planes only in 
// the root node. 
void BC_SetRootSeparationAlgorithm(BCTree &T, bool (*SeparationAlgorithm)(BCTree &T))
{ T.root_separation_algorithm = SeparationAlgorithm; }

// Set the heuristic that tries to obtain primal/feasible solution, 
// considering the fractional variables in the node
void BC_SetNodeHeuristic(BCTree &T, bool (*NodeHeuristic)(BCTree &T))
{ T.node_heuristic = NodeHeuristic; }

// Set the heuristic that tries to obtain primal/feasible solution,
// considering the fractional variables in the root node
void BC_SetRootHeuristic(BCTree &T, bool (*RootHeuristic)(BCTree &T))
{ T.root_heuristic = RootHeuristic; }

// returns the value of linear programming variable of the last solved LP.
double BC_GetLpVariableValue(BCTree &T,Lp::Col v)
{  return(T.lp.primal(v)); }

// print some information in each iteration if PrintIterations is true
void BC_SetPrintIterations(BCTree &T,bool PrintIterations)
{ T.printiterations = PrintIterations; }

// return the current number of branch and cut nodes.
unsigned long BC_GetNumberBBNodes(BCTree &T)
{ return(T.num_bbnodes); } // possibly used in cutting plane subroutines

// return the current number of iterations  
unsigned long BC_GetNumberIterations(BCTree &T)
{ return(T.iterations); }

// return the index of the lp variable
unsigned long BC_GetLpVarIndex(BCTree &T, Lp::Col v)
{ return(T.lp.id(v)); }

  
