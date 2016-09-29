/* 

     Flávio K. Miyazawa
     Envie erros/comentários para fkm@icunicamp.br.
*/

#include <iomanip>
#include <sstream>
#include <iostream>
#include "readgraph.h"
#include "viewgraph.h"
#include<lemon/list_graph.h>
#include <lemon/concepts/digraph.h>
#include <lemon/preflow.h>

using namespace std;

int main()
{

  ListGraph g;
  ListGraph::EdgeMap<double> weight(g);
  ListGraph::NodeMap<string> vname(g);
  ListGraph::NodeMap<double> posx(g);
  ListGraph::NodeMap<double> posy(g);
  ListGraph::NodeMap<int> vcolor(g);
  ListGraph::EdgeMap<int> acolor(g);

  // le o grafo do arquivo gr_50. O true faz cada aresta ser ida e volta
  ReadEuclideanGraph(g,vname,weight,posx,posy,"gr_berlin52");

  // pinta os vertices de azul
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) vcolor[v] = BLUE;

  // pinta as arestas de azul e coloca o peso da aresta como label
  for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
    std::ostringstream oss;
    oss << setprecision(2) << weight[a];
    acolor[a] = RED;
  }
  
  ViewEuclideanGraph(g,vname,posx,posy,vcolor,acolor);
  return 0;
}
