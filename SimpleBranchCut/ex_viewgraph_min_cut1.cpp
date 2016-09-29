
#include "viewgraph.h"
#include <lemon/preflow.h>

double cutValue (const ListDigraph& g,
              const ListDigraph::NodeMap<bool>& cut,
              const ListDigraph::ArcMap<double>& cap) {

  double c=0.0;
  for(ListDigraph::ArcIt e(g); e!=INVALID; ++e) {
    if (cut[g.source(e)] && !cut[g.target(e)]) c+=cap[e];
  }
  return c;
}

string DoubleToString(double x)
{
  std::ostringstream oss;
  oss << x;
  return(oss.str());
}

int main()
{
  typedef ListDigraph::ArcMap<double> CapMap;
  typedef ListDigraph::NodeMap<bool> CutMap;
  ListDigraph g;
  CutMap cut(g);
  CapMap weight(g);
  typedef ListDigraph::Node Node;
  typedef ListDigraph::NodeIt NodeIt;
  typedef ListDigraph::Arc Arc;

  typedef Preflow<ListDigraph, CapMap> PType;

  ListDigraph::NodeMap<int> vcolor(g);
  ListDigraph::NodeMap<string> vname(g);
  ListDigraph::ArcMap<int> acolor(g);
  ListDigraph::ArcMap<string> aname(g);


  // Cria cada vértice, seu nome e sua cor
  Node n1=g.addNode(); vname[n1]="1"; vcolor[n1]=BLUE;
  Node n2=g.addNode(); vname[n2]="2"; vcolor[n2]=BLUE;
  Node n3=g.addNode(); vname[n3]="3"; vcolor[n3]=BLUE;
  Node n4=g.addNode(); vname[n4]="4"; vcolor[n4]=BLUE;
  Node n5=g.addNode(); vname[n5]="5"; vcolor[n5]=BLUE;

  Arc a;
  // Cria cada aresta, define sua cor inicial, seu peso, e seu label (igual ao peso)
  a=g.addArc(n1,n2); acolor[a]=BLUE; weight[a]=1.1; aname[a]=DoubleToString(weight[a]); 
  a=g.addArc(n2,n3); acolor[a]=BLUE; weight[a]=1.2; aname[a]=DoubleToString(weight[a]); 
  a=g.addArc(n3,n5); acolor[a]=BLUE; weight[a]=3.3; aname[a]=DoubleToString(weight[a]); 
  a=g.addArc(n5,n4); acolor[a]=BLUE; weight[a]=1.3; aname[a]=DoubleToString(weight[a]); 
  a=g.addArc(n4,n1); acolor[a]=BLUE; weight[a]=1.4; aname[a]=DoubleToString(weight[a]); 
  a=g.addArc(n2,n4); acolor[a]=BLUE; weight[a]=2.2; aname[a]=DoubleToString(weight[a]); 
  a=g.addArc(n3,n4); acolor[a]=BLUE; weight[a]=1.5; aname[a]=DoubleToString(weight[a]); 
  
  // Encontrar o corte separando os vertices n2 e n4
  PType preflow_test(g, weight, n2, n4);
  preflow_test.run();

  CutMap min_cut(g);
  preflow_test.minCutMap(cut);
  double min_cut_value=cutValue(g,cut,weight);
  
  for (ListDigraph::ArcIt a(g); a!=INVALID; ++a) 
    if (cut[g.source(a)] && !cut[g.target(a)]) acolor[a]=RED;
  
  cout << "Nos do lado A:" << endl;
  for (ListDigraph::NodeIt v(g); v!=INVALID; ++v)
    {
      if (cut[v])  {
	cout << g.id(v)+1 << ", " << endl;
	vcolor[v] = BLUE;
      }
    }
  cout << "\nNos do lado B:" << endl;
  for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) {
    if (!cut[v]) {
      cout << g.id(v)+1 << ", " << endl;
      vcolor[v] = RED;
    }
  }
  
  cout << "Corte tem valor: " << min_cut_value << endl;

  ViewDigraph(g,VIEW_DOT,vname,aname,vcolor,acolor);
  return 0;
}
