#include <iostream> 
#include <fstream> 
#include "adjacencymatrix.h"

int globalcounter=0;

void ADJMAT_FreeNotNull(void *p){  if (p) free(p);  }

// Define an adjacency matrix, so as we have a fast access for the edges of a graph,
// given a pair of vertices. This is mainly used in the subroutine 2opt, for 
// dense graphs.
// The adjacency matrix is stored in a strictly triangular inferior matrix.

AdjacencyMatrix::AdjacencyMatrix(ListGraph &graph,EdgeWeight &graphweight,double NonEdgValue):
  Node2Index(graph),Edge2Index(graph)
{
  int i;
  g = &graph;
  NonEdgeValue = NonEdgValue;
  weight = &graphweight;
  Nnodes = countNodes(graph); // number of nodes in the input graph
  Nedges = countEdges(graph); // number of edges in the input graph
  Nmatrix = (Nnodes*(Nnodes-1))/2; // no. of edges/elem. in strict. triang. inf. matrix

  AdjMatrix = (double *) malloc(sizeof(double)*Nmatrix);
  Index2Node = (Node *) malloc(sizeof(Node)*Nnodes);
  Index2Edge = (Edge *) malloc(sizeof(Edge)*Nedges);

  if ((AdjMatrix==NULL)||(Index2Node==NULL)||(Index2Edge==NULL)) { 
    cout << "Out of memory in constructor of AdjacencyMatrix\n"; 
    ADJMAT_FreeNotNull(AdjMatrix); ADJMAT_FreeNotNull(Index2Node); ADJMAT_FreeNotNull(Index2Edge);
    exit(0);}

  i = 0;
  for (NodeIt v(*g); v != INVALID; ++v) {
    Index2Node[i] = v;
    AdjacencyMatrix::Node2Index[v]=i;
    i++;
  }

  // Initially all edges have infinity weight
  for (int i=0;i<Nmatrix;i++) AdjMatrix[i] = NonEdgeValue;
  // Then, update the existing edges with the correct weight
  for (EdgeIt e(graph); e != INVALID; ++e) {
    Node u,v;    int i_u,i_v;
    u = graph.u(e);  v = graph.v(e);  // obtain the extremities of e
    i_u = Node2Index[u];
    i_v = Node2Index[v];
    if (i_u > i_v) AdjMatrix[i_u*(i_u-1)/2+i_v] = graphweight[e];
    else if (i_u < i_v) AdjMatrix[i_v*(i_v-1)/2+i_u] = graphweight[e];
    else {
      cout << "Out of memory in constructor of AdjacencyMatrix\n"; 
      exit(0);}
  }
}

double AdjacencyMatrix::Cost(Node u,Node v)
{
  int i_u,i_v;
  i_u = Node2Index[u];
  i_v = Node2Index[v];
  globalcounter++;
  //cout << "globalcounter1 = " << globalcounter << endl;
  try{
    if (i_u > i_v) return(AdjMatrix[i_u*(i_u-1)/2+i_v]);
    else if (i_u < i_v) return(AdjMatrix[i_v*(i_v-1)/2+i_u]);
    else return(NonEdgeValue);
  }catch (...) {
    cout << "LANCOU UMA EXCECAO: " << globalcounter << endl;
    exit(1);
  }
}

double AdjacencyMatrix::Cost(Edge e)
{
  int i_u,i_v;
  Node u,v;
  globalcounter++;
  //cout << "globalcounter2 = " << globalcounter << endl;
  u = (*g).u(e);  v = (*g).v(e);
  i_u = Node2Index[u];
  i_v = Node2Index[v];
  if (i_u > i_v) return(AdjMatrix[i_u*(i_u-1)/2+i_v]);
  else return(AdjMatrix[i_v*(i_v-1)/2+i_u]);
}


AdjacencyMatrix::~AdjacencyMatrix()
{
  ADJMAT_FreeNotNull(AdjMatrix); ADJMAT_FreeNotNull(Index2Node); ADJMAT_FreeNotNull(Index2Edge);
}

