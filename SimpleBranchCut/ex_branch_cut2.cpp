//----------------------------------------------------------------------
// Example of program that uses the LEMON Library, Linear Programming
// and the branch and cut code in this directory to solve a integer linear
// program.
//
// To use this program, you must install the LEMON Library and a linear programming
// package used by LEMON (like glpk, clp, cplex,...). 
// Read the README file for more information.
//
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------

#include "branch_cut.h"
//#include "branch_cut.cpp"
int main()
{
  Lp::Col x1,x2;
  
  // Declare a branch and cut tree with nvar variables and maximization obj. function
  BCTree T(2,BC_MAXIMIZE);

  // Obtain the two variables
  x1 = BC_GetNewVariable(T);
  x2 = BC_GetNewVariable(T);

  // Print no informations in each iteration
  BC_SetPrintIterations(T,false);

  // Define the lower and upper bound of each variable
  BC_SetVarBounds(T , x1 , 0.0 , BC_INF); 
  BC_SetVarBounds(T , x2 , 0.0 , BC_INF); 
  
  // Define the objective function
  BC_ObjectiveFunction(T, 3*x1 + 4*x2 ); 

  // Add three constraints
  BC_AddConstraint(T, -3*x1 + 2*x2  <= 2 );
  BC_AddConstraint(T,    x1 + 3*x2  <= 11);
  BC_AddConstraint(T,    x1 +   x2  <= 6);

  // Solve the problem and verify if it obtained an optimum solution
  if (BC_Solve(T)==BC_OPTIMUM) {
    cout << "Objective function value: " << BC_GetSolutionValue(T) << endl;
    cout << "x1 = " << BC_GetSolutionVariableValue(T,x1) << endl;
    cout << "x2 = " << BC_GetSolutionVariableValue(T,x2) << endl;
  } else {cout << "No optimum solution found." << endl;}
  return 0;
}
