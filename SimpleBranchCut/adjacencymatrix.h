#ifndef ADJMAT_DEFINE
#define ADJMAT_DEFINE

#include<lemon/list_graph.h>
#include<lemon/concepts/graph.h>
#include<lemon/concepts/maps.h>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<string>
#include <fstream> 

using namespace lemon;
using namespace std;


typedef ListGraph::Node Node;
typedef ListGraph::NodeIt NodeIt;
typedef ListGraph::EdgeIt EdgeIt;
typedef ListGraph::Edge Edge;
typedef ListGraph::EdgeMap<double> EdgeWeight;
typedef ListGraph::NodeMap<string> NodeName;
typedef ListGraph::EdgeMap<string> EdgeName;
typedef ListGraph::NodeMap<bool> CutMap;
typedef ListGraph::NodeMap<int> NodeColor;
typedef ListGraph::EdgeMap<int> EdgeColor;
typedef ListGraph::EdgeMap<int> EdgeIndex;
typedef ListGraph::NodeMap<int> NodeIndex;

class AdjacencyMatrix {
public:
  AdjacencyMatrix(ListGraph &graph,EdgeWeight &graphweight,double NonEdgeValue);
  ~AdjacencyMatrix();
  double *AdjMatrix;
  ListGraph *g;
  EdgeWeight *weight;
  int Nnodes,Nedges,Nmatrix;
  double NonEdgeValue;
  Node *Index2Node;
  Edge *Index2Edge;
  double Cost(Node,Node);
  double Cost(Edge);
  NodeIndex Node2Index;
  EdgeIndex Edge2Index;
};


#endif

  
