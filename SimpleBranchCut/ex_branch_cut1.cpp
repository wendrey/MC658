//----------------------------------------------------------------------
// Example of program that uses the LEMON Library, Linear Programming
// and the branch and cut code in this directory to solve a integer linear
// program.
//
// To use this program, you must install the LEMON Library and a linear programming
// package used by LEMON (like glpk, clp, cplex,...)
// Read the README file for more information.
//
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------



#include "branch_cut.h"

int main()
{
  int nvar=2;
  Lp::Col x1,x2;
  
  // Declare a branch and cut tree with nvar variables and minimization obj. function
  BCTree T(nvar,BC_MINIMIZE);

  // Obtain the two variables
  x1 = BC_GetNewVariable(T);
  x2 = BC_GetNewVariable(T);

  // Print informations in each iteration
  BC_SetPrintIterations(T,true);

  // Define the lower and upper bound of each variable
  BC_SetVarBounds(T , x1 , 0.0 , 4.0);    // var. x1 has lower and upper bounds [0,4]
  BC_SetVarBounds(T , x2 , 0.0 , BC_INF); // var. x2 has lower and upper bounds [0,inf]
  
  // Add two constraints
  BC_AddConstraint(T,   x1 +    x2 <= 5.9 );
  BC_AddConstraint(T, 2*x1 + 3*x2  <= 12.9);

  // Define the objective function
  BC_ObjectiveFunction(T, -2*x1 - x2 ); 

  // Solve the problem and verify if it obtained an optimum solution
  if (BC_Solve(T)==BC_OPTIMUM) {
    cout << "Objective function value: " << BC_GetSolutionValue(T) << endl;
    cout << "x1 = " << BC_GetSolutionVariableValue(T,x1) << endl;
    cout << "x2 = " << BC_GetSolutionVariableValue(T,x2) << endl;
  } else {cout << "No optimum solution found." << endl;}
  return 0;
}
