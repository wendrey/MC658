#include <iostream>
#include <string>
#include <stdio.h>
#include <vector>
#include <queue>
#include <iomanip>

using namespace std;

typedef struct {
  int year,month,day;
} Date;

bool comparedate (const Date &a,const Date &b) { 
  if (a.year < b.year) return(true);
  if (a.year > b.year) return(false);
  if (a.month < b.month) return(true);
  if (a.month > b.month) return(false);
  if (a.day < b.day) return(true);
  if (a.day > b.day) return(false);
  return(false);
}

void printdate(const Date &d) 
{
  cout << setfill('0');
  cout << setw(4) << d.year << '/' 
       << setw(2) << d.month << '/' 
       << setw(2) << d.day;
  cout << setfill(' ');
}

  
void generaterandomdatevector(vector<Date> &d,int n)
{
  int i;
  for (i=0;i<n;i++) {
    d[i].year = 1900 + rand()%100;
    d[i].month = 1 + rand()%12;
    d[i].day = 1 + rand()%28;
  }
}


void printdatevector(vector<Date> &d,int n)
{
  int i;
  for (i=0;i<n;i++) {
    cout << setw(2) << i+1 << "  "; 
    printdate(d[i]);  
    cout << endl;
  }
}

#define N 20
int main()
{
  vector<Date> d(N);
  Date x;
  int i,n;
  n = N;
  generaterandomdatevector(d,n);
  printdatevector(d,n);
  cout << "Numero de elementos: " << d.size() << endl;
  cout << "Vetor original" << endl;
  printdatevector(d,n);
  sort(d.begin(),d.begin()+n,comparedate);
  cout << "Vetor ordenado" << endl;
  printdatevector(d,n);
 
}
