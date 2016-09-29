//----------------------------------------------------------------------
// Example of program that uses the LEMON Library, Linear Programming
//
// To use this program, you must install the LEMON Library and a linear programming
// package used by LEMON (like glpk, clp, cplex,...)
//
// This program is the formulation of the Traveling Salesman Problem 
// in pages 144 to 148 of the slides in the link (in portuguese)
//
// http://www.ic.unicamp.br/~fkm/lectures/proglin.pdf
//
// To compile this program, read the README file
//
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------
//
#include <lemon/lp.h>
using namespace lemon;
using namespace std;
int main()
{
  // Declare a linear programming variable (initially empty)
  Lp lp;

  // Add linear programming variables: xa,...,xl
  Lp::Col xa = lp.addCol();
  Lp::Col xb = lp.addCol();
  Lp::Col xc = lp.addCol();
  Lp::Col xd = lp.addCol();
  Lp::Col xe = lp.addCol();
  Lp::Col xf = lp.addCol();
  Lp::Col xg = lp.addCol();
  Lp::Col xh = lp.addCol();
  Lp::Col xi = lp.addCol();
  Lp::Col xj = lp.addCol();
  Lp::Col xk = lp.addCol();
  Lp::Col xl = lp.addCol();

  // Each variable has bounds in [0,1]
  lp.colLowerBound( xa, 0 );   lp.colUpperBound( xa, 1 );
  lp.colLowerBound( xb, 0 );   lp.colUpperBound( xb, 1 );
  lp.colLowerBound( xc, 0 );   lp.colUpperBound( xc, 1 );
  lp.colLowerBound( xd, 0 );   lp.colUpperBound( xd, 1 );
  lp.colLowerBound( xe, 0 );   lp.colUpperBound( xe, 1 );
  lp.colLowerBound( xf, 0 );   lp.colUpperBound( xf, 1 );
  lp.colLowerBound( xg, 0 );   lp.colUpperBound( xg, 1 );
  lp.colLowerBound( xh, 0 );   lp.colUpperBound( xh, 1 );
  lp.colLowerBound( xi, 0 );   lp.colUpperBound( xi, 1 );
  lp.colLowerBound( xj, 0 );   lp.colUpperBound( xj, 1 );
  lp.colLowerBound( xk, 0 );   lp.colUpperBound( xk, 1 );
  lp.colLowerBound( xl, 0 );   lp.colUpperBound( xl, 1 );

  // Define the objective function and set as minimization 
  lp.min();
  lp.obj( 7*xa + 5*xb + 8*xc + 2*xd + 3*xe + 2*xf + 
	  9*xg + 4*xh + 3*xi + 8*xj + 2*xk + 6*xl);   

  // Add constraints (see the sequence of slides). In fact, to obtain the
  // solution of each slide, I have to run the solve lp.solve() after inserting 
  // each constraint. 
  lp.addRow( xa +   xb == 2 );
  lp.addRow( xa + xc + xd + xe == 2 );
  lp.addRow( xb + xe + xf + xg == 2 );
  lp.addRow( xd + xf + xh + xi == 2 );
  lp.addRow( xc + xh + xj + xk == 2 );
  lp.addRow( xg + xi + xj + xl == 2 );
  lp.addRow( xk + xl == 2 );
  lp.solve();
  lp.addRow( xc + xd + xf + xg >= 2 );
  lp.addRow( xc + xg + xh + xi >= 2 );
  lp.solve();
  
  // Print the value of the variables
  if (lp.primalType() == Lp::OPTIMAL) {
    cout << "Value of the objective function: " << lp.primal() << endl;
    cout << "x_a = " << lp.primal(xa) << endl;
    cout << "x_b = " << lp.primal(xb) << endl;
    cout << "x_c = " << lp.primal(xc) << endl;
    cout << "x_d = " << lp.primal(xd) << endl;
    cout << "x_e = " << lp.primal(xe) << endl;
    cout << "x_f = " << lp.primal(xf) << endl;
    cout << "x_g = " << lp.primal(xg) << endl;
    cout << "x_h = " << lp.primal(xh) << endl;
    cout << "x_i = " << lp.primal(xi) << endl;
    cout << "x_j = " << lp.primal(xj) << endl;
    cout << "x_k = " << lp.primal(xk) << endl;
    cout << "x_l = " << lp.primal(xl) << endl;
  } else {cout << "No optimum solution found." << endl;}
  return 0;
}
