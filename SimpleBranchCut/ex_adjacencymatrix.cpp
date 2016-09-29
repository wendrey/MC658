#include "readgraph.h"
#include "viewgraph.h"
#include "adjacencymatrix.h"
#include <fstream>
#include <iostream>
#include <iomanip>


// Only to see if a file exists. It (try to) open and close the file.
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

  ListGraph g;
  EdgeWeight weight(g);
  EdgeWeight lpvariablevalue(g);
  NodeName vname(g);
  EdgeName ename(g);
  EdgeName enametoview(g);
  EdgeColor ecolor(g);
  NodeColor vcolor(g);
  ListGraph::NodeMap<double> posx(g),posy(g);

  
  srand(1);
 
  if (argc!=2) {
    cout << "Usage: " << argv[0] << " <graph_filename>" << endl;
    exit(0);
  } else if (!FileExists(argv[1])) {
    cout << "File " << argv[1] << " does not exist." << endl;
    exit(0);
  }
  if (!ReadGraph(g,vname,weight,posx,posy,argv[1])) {
    cout << "Error reading graph file " << argv[1] << "." << endl;
    exit(0);
  }

  cout << "number of nodes = " << countNodes(g) << endl;
  cout << "number of edges  = " << countEdges(g) << endl;
  
  for (EdgeIt a(g); a!=INVALID; ++a) {
    cout << "w(" << g.id(g.u(a)) << "," << g.id(g.v(a)) << ")=" << weight[a] << endl;
  }

  cout << "---------------------------" << endl;

  
  AdjacencyMatrix Adj(g,weight,-1.0);
  cout << "    ";
  for (ListGraph::NodeIt u(g); u!=INVALID; ++u) cout << setw(5) << g.id(u) << "  ";
  cout << endl;

  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    cout << g.id(v) << "  ";
    for (ListGraph::NodeIt u(g); u!=INVALID; ++u) 
      cout << setw(5) << Adj.Cost(u,v) << "  ";
    cout << endl;
  }
  return 0;
 }
 
