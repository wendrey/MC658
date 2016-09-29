//----------------------------------------------------------------------
// Example of program that uses the LEMON Library, Linear Programming
//
// To use this program, you must install the LEMON Library and a linear programming
// package used by LEMON (like glpk, clp, cplex,...)
//
// This program solve the diet problem using the linear programming formulation
// given in the page 12 of the slides in the link (in portuguese)
//
// http://www.ic.unicamp.br/~fkm/lectures/proglin.pdf
//
// To compile this program, read the README file
//
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------

#include <lemon/lp.h>
using namespace lemon;
using namespace std;
int main()
{ 

  // Declaring the LP variable
  Lp lp;
  typedef Lp::Col LPvar;
  
  // Add three variables
  LPvar x1 = lp.addCol();  
  LPvar x2 = lp.addCol();  
  LPvar x3 = lp.addCol();

  // Set the objective function
  lp.obj( 4.78*x1 + 1.48*x2 + 0.85*x3 );   

  // set the problem as a minimization problem
  lp.min();

  // Add two constraints
  lp.addRow( 800 <= 110*x1 + 170*x2 + 1150*x3 <= 1200 );
  lp.addRow(  10 <=  29*x1 +  15*x2           <= 18 );

  // Add lower bounds for each variable
  lp.colLowerBound( x1, 0 );   
  lp.colLowerBound( x2, 0 );
  lp.colLowerBound( x3, 0 );

  // Solve the linear program
  lp.solve();

  // If found optimum solution, print the values of the variables
  if (lp.primalType() == Lp::OPTIMAL) {
    cout << "Value of the objective function: " << lp.primal() << endl;
    cout << "x1 = " << lp.primal(x1) << endl;
    cout << "x2 = " << lp.primal(x2) << endl;
    cout << "x3 = " << lp.primal(x3) << endl;
  } else {cout << "No solution found." << endl;}
  return 0;
}
