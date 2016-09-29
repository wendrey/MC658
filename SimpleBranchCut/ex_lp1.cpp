//----------------------------------------------------------------------
// Example of program that uses the LEMON Library, Linear Programming
//
// To use this program, you must install the LEMON Library and a linear programming
// package used by LEMON (like glpk, clp, cplex,...)
//
// Solving a small linear program
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
  // Declaring/generating an LP (initially empty)
  Lp lp;

  // Adding two variables (columns of the LP)
  Lp::Col x1 = lp.addCol();
  Lp::Col x2 = lp.addCol();

  // Adding constraints
  lp.addRow(   x1 +    x2 <= 5.9 );
  lp.addRow( 2*x1 + 3*x2  <= 12.9);

  // Setting lower and upper bounds for each variable
  lp.colLowerBound( x1, 0 );   lp.colUpperBound( x1, 4 );
  lp.colLowerBound( x2, 0 );

  // Defining the objective function and its objective sense
  lp.obj( -2*x1 - x2 );   lp.min();

  // Solving the linear program
  lp.solve();

  // Printing the value of the variables
  if (lp.primalType() == Lp::OPTIMAL) {
    cout << "Value of the objective function: " << lp.primal() << endl;
    cout << "x" << lp.id(x1) << " = " << lp.primal(x1) << endl;
    cout << "x" << lp.id(x2) << " = " << lp.primal(x2) << endl;
  } else {cout << "No optimum solution found." << endl;}
  return 0;
}
