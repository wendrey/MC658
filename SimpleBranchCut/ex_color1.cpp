
// A simple example to use the definitions in color.h
#include <iostream>
#include "color.h"
#include <stdio.h>
using namespace std;
int main()
{
  Color c;
  c = BLUE;
  cout << ColorName[c] << endl;
}

