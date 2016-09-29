
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


#ifndef BC_DEFINE
#define BC_DEFINE
#include <queue>
#include <lemon/lp.h>
#include <time.h>
using namespace std;
using namespace lemon;

// Este programa foi feito apenas para efeitos didáticos e não foi projetado
// para ter bom desempenho. 
#define BC_EPS 0.001
#define BC_MINIMIZE 0
#define BC_MAXIMIZE 1
// Teste com números muito grandes faz LP do LEMON+GLPK dar erro para alguns problemas
//#define BC_INF (numeric_limits<double>::max())
//#define BC_MINUS_INF (-numeric_limits<double>::max())
#define BC_INF 1000000000
#define BC_MINUS_INF -1000000000

typedef enum 
  {BC_OPTIMUM,    // Integer optimum solution
   BC_INTEGER,    // Obtained an integer solution, but possibly not optimum
   BC_FEASIBLE,   // Could not obtain any integer solution
   BC_UNBOUNDED,  // Integer program is unbounded
   BC_INFEASIBLE, // Integer program is infeasible
   BC_UNDEFINED   // Could not obtain anything
  } BC_Status;
  
class BCNode {
public:
  BCNode(int maxvar);
  class BCTree *T;  // para obter informacoes da arvore
  double lb_value,ub_value;
  vector <double> lb_var;
  vector <double> ub_var;
  int ind_frac;
  double val_frac;
  int nvarfrac;
  BC_Status lp_status; 
};

// class used in the priority queue to compare nodes.
class CompareBCNode {
public:
  bool operator()(BCNode &a, BCNode &b);
};

// This branch_and_cut code only consider global cuts. 
// That is, cuts that are valid for all nodes.
// To simplify the code, all variables must be integer.
// It is an exercise to change this code to solve mixed integer programs
class BCTree {
 public:
  BCTree(int maxvar,int sense); // constructor with obj sense and max. num. var.
  int obj_sense; // says if the problem is of minimization (BC_MINIMIZE) or maximization (BC_MAXIMIZE)
  double best_lb,best_ub; // stores the global lower and upper bounds
  int max_var, // maximum number of variables that are possible to use
    num_var; // maximum number of variables being used in the linear program
  // Branch and Cut tree represented by a priority queue
  priority_queue <BCNode,vector<BCNode>,CompareBCNode> BCHeap; 
  Lp lp;
  vector <Lp::Col> lpvar;    // variables of the linear program
  vector <double> global_lb_var;  // original lower bound of the variables
  vector <double> global_ub_var;  // original upper bound of the variables
  vector <double> best_x;  // best solution found so far
  bool (*node_separation_algorithm)(BCTree &); // initially no separation algorithm
  bool (*root_separation_algorithm)(BCTree &); // initially no separation algorithm
  bool (*node_heuristic)(BCTree &); // initially no node heuristic
  bool (*root_heuristic)(BCTree &); // initially no root heuristic
  void *problem_data;  // pointer to the data of the user problem
  bool printiterations; // flag if print informations for each iteration
  unsigned long iterations, // count the number of iterations
    num_bbnodes, // count the number of branch and bound/cut nodes
    num_constraints; // count the total number of inserted constraints
  time_t initial_time,ending_time;
};

// Return the elapsed time since the starting time of the BC_Solve routine.
long BC_GetSolverTime(BCTree &T);

// Return true if variable x is fractional (within a certain small error).
bool BC_IsFrac(double x);

// Set the lower bound and the upper bound of a variable
void BC_SetVarBounds(BCTree &T,const Lp::Col var, double lb,double ub);

// Add a global linear constraint to the formulation
void BC_AddConstraint(BCTree &T,const Lp::Constr &C);

// Set the objective function of the problem. It can be BC_MINIMIZE or BC_MAXIMIZE.
// Note that all branch and cut program is for minimization. So, if the user problem
// is of maximization, the objective function is multiplied by (-1).
void BC_ObjectiveFunction(BCTree &T,const Lp::Expr &E);

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
BC_Status BC_Solve(BCTree &T);

// Return a new linear programming variable. Note that the maximum number of
// variables is defined when you declare the branch and cut tree T.
Lp::Col BC_GetNewVariable(BCTree &T);

// Return the value of the linear programming variable.
double BC_GetSolutionVariableValue(BCTree &T,const Lp::Col var);

// return the value of the objective function. Note that for maximization problem
// the function was multiplied by (-1). So, the returned value must also be
// multiplied by (-1).
double BC_GetSolutionValue(BCTree &T);

// return the global lower bound of the problem. Note that for maximization problem
// the objective function was transformed into a minimization function. So, the 
// current lower bound is in fact an upper bound.
void BC_SetGlobalLowerBound(BCTree &T,double value);

// return the global upper bound of the problem. See note for the global 
// lower bound above.
void BC_SetGlobalUpperBound(BCTree &T,double value);

// Place the pointer of the problem in the branch and cut structure
// Possibly, if you use some heuristic based on the linear programming variables,
// you may need to access the original data of the problem.
void BC_SetProblemData(BCTree &T,void *ProblemData);

// Get the pointer of the problem from the branch and cut structure
void *BC_GetProblemData(BCTree &T);

// Set the routine that insert global cutting planes (linear constraints) in 
// the linear program, when solving a node (even if the node is the root node).
void BC_SetNodeSeparationAlgorithm(BCTree &T,bool (*SeparationAlgorithm)(BCTree &T));

// Set the routine that insert global cutting planes only in the root node.
// In some problems, it is interesting to insert cutting planes only in 
// the root node. 
void BC_SetRootSeparationAlgorithm(BCTree &T,bool (*SeparationAlgorithm)(BCTree &T));

// Set the heuristic that tries to obtain primal/feasible solution, 
// considering the fractional variables in the node
void BC_SetNodeHeuristic(BCTree &T,bool (*NodeHeuristic)(BCTree &T));

// Set the heuristic that tries to obtain primal/feasible solution,
// considering the fractional variables in the root node
void BC_SetRootHeuristic(BCTree &T,bool (*RootHeuristic)(BCTree &T));

// returns the value of linear programming variable of the last solved LP.
double BC_GetLpVariableValue(BCTree &T,Lp::Col v);

// print some information in each iteration if PrintIterations is true
void BC_SetPrintIterations(BCTree &T,bool PrintIterations);

// return the current number of branch and cut nodes.
unsigned long BC_GetNumberBBNodes(BCTree &T);

// return the current number of iterations  
unsigned long BC_GetNumberIterations(BCTree &T);

// return the index of the lp variable
unsigned long BC_GetLpVarIndex(BCTree &T,Lp::Col v);
#endif
