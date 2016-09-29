#include "readgraph.h"
#include <math.h>
#include "color.h"
#include <iostream> 
#include <fstream> 
#include <cstdlib>
#include <cstring>

typedef struct {
  ListDigraph::Node v;
  string vert_arquivo;
} tipomapadigraph;

typedef struct {
  ListGraph::Node v;
  string vert_arquivo;
} tipomapagraph;

string IntToString(int number)
{ stringstream s;  s << number;  return s.str(); }


int comparamapadigraph (const void *pa,const void *pb) { 
  tipomapadigraph *a,*b;
  a = (tipomapadigraph *) pa;
  b = (tipomapadigraph *) pb;
  if (a->vert_arquivo > b->vert_arquivo)return(1);
  if (a->vert_arquivo < b->vert_arquivo)return(-1);
  return(0);
}

int comparamapagraph (const void *pa,const void *pb) { 
  tipomapagraph *a,*b;
  a = (tipomapagraph *) pa;
  b = (tipomapagraph *) pb;
  if (a->vert_arquivo > b->vert_arquivo)return(1);
  if (a->vert_arquivo < b->vert_arquivo)return(-1);
  return(0);
}

void PulaBrancoComentario(ifstream &ifile)
{
  char c;
  string line;
  while (!ifile.eof()) {
    c = ifile.get();
    while ((c==' ') && (!ifile.eof())) c = ifile.get();
    if (!ifile.eof()) {
      if (c=='#')  getline(ifile,line);
      else {ifile.unget(); break; }
    }
  }
}


bool WriteGraphGraphviz(ListGraph &g,
		   ListGraph::NodeMap<string> &vname, // vertex names
		   ListGraph::EdgeMap<string> &ename,  // edge names
		   ListGraph::NodeMap<int> &vcolor,   // vertex colors
		   ListGraph::EdgeMap<int> &acolor,    // edge colors
		   string filename)
{
  ofstream out;
  string linha;
  
  out.open(filename.c_str());
  if (out.is_open()) return(false);

  out << "graph g {\n";
  out << "\tsize = \"8, 11\";\n";
  out << "\tnode [shape = \"circle\"];\n";
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    linha = "\t";    linha += vname[v].c_str();   linha += " [color=";
    linha += ColorName[vcolor[v]].c_str();        linha += "];\n";
    out << linha;
  }
  for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
    if (acolor[a]!=WHITE) {
      linha = "\t";
      linha += vname[g.u(a)].c_str() ;
      linha += "  -- ";
      linha += vname[g.v(a)].c_str();
      linha += " [label = \"";  linha += ename[a].c_str();  linha += "\"";
      linha += ", color=\""; linha += ColorName[acolor[a]].c_str(); linha += "\" ];\n";
      out << linha;
    }
  }
  out << "}\n";
  out.close();
  return(true);
}



bool ReadListDigraph(ListDigraph &g,
		 ListDigraph::NodeMap<string>& vname,
		 ListDigraph::ArcMap<double>& custo,
		 const bool dupla,char *filename)
{
  ifstream ifile;
  int i,n,m;
  ListDigraph::Node nu,nv;
  double peso;
  ListDigraph::Arc a;
  char nomeu[100],nomev[100];
  
  ifile.open(filename);  if (!ifile) return(false);
  PulaBrancoComentario(ifile);
  ifile >> n;    ifile >> m;
  if ((m<0)||(n<0)) {
    cout << "File " << filename << " is not a digraph given by edges.\n";
    exit(0);
  }
  tipomapadigraph mapa[n],k,*q;
  for (i=0;i<n;i++) {
    if (ifile.eof()) {
      cout << "Reached unexpected end of file " << filename << ".\n";
      exit(0);
    }
    ifile >> nomev;
    mapa[i].vert_arquivo = nomev;
    mapa[i].v = g.addNode();
    vname[mapa[i].v] = nomev;
  }
  qsort(mapa,n,sizeof(tipomapagraph),comparamapagraph);
  for (i=0;i<m;i++) {
    if (ifile.eof()) {
      cout << "Reached unexpected end of file " << filename << ".\n";
      exit(0);
    }
    ifile >> nomeu;  ifile >> nomev; ifile >> peso;
    // cout << "nomeu = " << nomeu << " nomev = " << peso << endl;
    
    k.vert_arquivo = nomeu; 
    q = (tipomapadigraph *) bsearch(&k,mapa,n,sizeof(tipomapagraph),comparamapagraph);
      if (q) nu = q->v;
      else { cout<<"ERRO: Vertice "<<q->vert_arquivo<<" nao encontrado.\n"<<endl;exit(0);}
      k.vert_arquivo = nomev; 
      q = (tipomapadigraph *) bsearch(&k,mapa,n,sizeof(tipomapagraph),comparamapagraph);
      if (q) nv = q->v;
      else { cout<<"ERRO: Vertice "<<q->vert_arquivo<<" nao encontrado.\n"<<endl;exit(0);}
      a = g.addArc(nu,nv);
      custo[a] = peso;
      if (dupla) {
	a = g.addArc(nv,nu);
	custo[a] = peso;
      }
  }
  ifile.close();
  return(true);
}

bool ReadGraph2(ListGraph &g,
	       ListGraph::NodeMap<string>& nodename,
	       ListGraph::EdgeMap<double>& custo,
	       char *filename)
{
  FILE *fp;
  int i,n,m;
  ListGraph::Node nu,nv;
  double peso;
  tipomapagraph mapa[5000],k,*q;
  ListGraph::Edge a;
  char nomeu[100],nomev[100];
  
  fp = fopen(filename,"r");
  if (fp==NULL) return(false);
  fscanf(fp,"%d %d",&n,&m);
  if (n>5000) {
    printf("ERRO: Numero de vertices do grafo maior 5000.\n");
    return(false);
  }
  for (i=0;i<n;i++) {
    fscanf(fp,"%s",nomev);
    mapa[i].vert_arquivo = nomev;
    mapa[i].v = g.addNode();
    nodename[mapa[i].v] = nomev;
    // cout << "node = " << mapa[i].vert_arquivo << endl;
  }
  qsort(mapa,n,sizeof(tipomapagraph),comparamapagraph);
  for (i=0;i<m;i++) {
    fscanf(fp,"%s %s %lf",nomeu,nomev,&peso);
    
    k.vert_arquivo = nomeu; 
    q = (tipomapagraph *) bsearch(&k,mapa,n,sizeof(tipomapagraph),comparamapagraph);
    if (q) nu = q->v;
    else { cout<<"ERRO: Vertice "<<q->vert_arquivo<<" nao encontrado.\n"<<endl;exit(0);}
    k.vert_arquivo = nomev; 
    q = (tipomapagraph *) bsearch(&k,mapa,n,sizeof(tipomapagraph),comparamapagraph);
    if (q) nv = q->v;
    else { cout<<"ERRO: Vertice "<<q->vert_arquivo<<" nao encontrado.\n"<<endl;exit(0);}
    a = g.addEdge(nu,nv);
    custo[a] = peso;
  }
  fclose(fp);
  return(true);
}

bool GenerateVertexPositions(ListGraph &g,
			     ListGraph::EdgeMap<double>& custo,
			     ListGraph::NodeMap<double>& posx,
			     ListGraph::NodeMap<double>& posy);

bool ReadListGraph(ListGraph &g,
		   ListGraph::NodeMap<string>& vname,
		   ListGraph::EdgeMap<double>& custo,
		   ListGraph::NodeMap<double>& posx,
		   ListGraph::NodeMap<double>& posy,
		   string filename)
{
  int i,n,m;
  ListGraph::Node nu,nv;
  double peso;
  tipomapagraph k,*q;
  ListGraph::Edge a;
  string nomeu,nomev;

  ifstream ifile;  
  ifile.open(filename.c_str());  
  if (!ifile) return(false);
  PulaBrancoComentario(ifile);
  ifile >> n;    ifile >> m;
  if (m<0) {
    cout << "File " << filename << " is not a graph given by edges.\n";
    exit(0);
  }
  tipomapagraph mapa[n];
  string STR;
  int nt;

  for (i=0;i<n;i++) {
    if (ifile.eof()) {
      cout << "Reached unexpected end of file " << filename << ".\n";
      exit(0);
    }
    getline(ifile,STR); 
    while (STR=="") getline(ifile,STR); 

    {
      string token;
      istringstream ins; // Declare an input string stream.
      ins.str(STR);        // Specify string to read.
      nt = 0;
      while( getline(ins, token, ' ') ) {
	if (nt==0) { 
	  nomev = token;
	  mapa[i].vert_arquivo = nomev;
	  mapa[i].v = g.addNode();
	  vname[mapa[i].v] = nomev;
	} else if (nt==1)  {
	  posx[mapa[i].v] = atof(token.c_str());
	} else if (nt==2)  {
	  posy[mapa[i].v] = atof(token.c_str());
	}
	nt++;
      }
    }
  }
  qsort(mapa,n,sizeof(tipomapagraph),comparamapagraph);
  for (i=0;i<m;i++) {
    if (ifile.eof()) {
      cout << "Reached unexpected end of file " << filename << ".\n";
      exit(0);
    }
    ifile >> nomeu;  ifile >> nomev; ifile >> peso;
    
    k.vert_arquivo = nomeu; 
    q = (tipomapagraph *) bsearch(&k,mapa,n,sizeof(tipomapagraph),comparamapagraph);
      if (q) nu = q->v;
      else { cout<<"ERRO: Vertice "<<q->vert_arquivo<<" nao encontrado.\n"<<endl;exit(0);}
      k.vert_arquivo = nomev; 
      q = (tipomapagraph *) bsearch(&k,mapa,n,sizeof(tipomapagraph),comparamapagraph);
      if (q) nv = q->v;
      else { cout<<"ERRO: Vertice "<<q->vert_arquivo<<" nao encontrado.\n"<<endl;exit(0);}
      a = g.addEdge(nu,nv);
      custo[a] = peso;
  }
  ifile.close();
  if (nt!=3) GenerateVertexPositions(g,custo,posx,posy);
  return(true);
}



bool ReadEuclideanGraph(ListGraph &g,
			ListGraph::NodeMap<string>& vname,
			ListGraph::EdgeMap<double>& custo,
			ListGraph::NodeMap<double>& posx,
			ListGraph::NodeMap<double>& posy,
			string filename)
{
  int i,n,m;
  ListGraph::Node nu,nv;
  ListGraph::Edge a;
  char nomev[100];
  ListGraph::Node v;
  double px,py;
  
  ifstream ifile;  ifile.open(filename.c_str());  if (!ifile) return(false);
  PulaBrancoComentario(ifile);
  ifile >> n;    ifile >> m;
  if (m!=-1) {
    printf("Wrong format in the euclidean graph of file %s.\n",filename.c_str());
    return(false);
  }
  for (i=0;i<n;i++) {
    ifile >> nomev;  ifile >> px; ifile >> py;
    v = g.addNode();
    vname[v] = nomev;
    posx[v] = px;
    posy[v] = py;
  }

  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    ListGraph::NodeIt u(g);
    u=v;
    for (++u; u!=INVALID; ++u) {
      a = g.addEdge(u,v);
      custo[a] = sqrt((posx[u]-posx[v])*(posx[u]-posx[v]) + 
		      (posy[u]-posy[v])*(posy[u]-posy[v]));
    }
  }
  ifile.close();
  return(true);
}

/* Dado um texto e um padrao, devolve a posicao no
   texto onde o padrao ocorre. Devolve (-1) se nao achar. */
int gr_busca_padrao(char *texto,char *busca)
{
  int i,j,tam_texto,tam_busca;
  char *p1,*p2;
  bool achou;
  tam_texto = strlen(texto);
  tam_busca = strlen(busca);
  for (i=0;i<tam_texto - tam_busca + 1;i++) {
      p1 = &texto[i];
      p2 = busca;
      for (j=0,achou = true;j<tam_busca && achou == true;j++)
	if (p1[j] != p2[j]) achou = false;
      if (achou) return(i);
  }
  return(-1);
}


// This routine use the neato program to generate positions.
bool GenerateVertexPositions(ListGraph &g,
			     ListGraph::EdgeMap<double>& custo,
			     ListGraph::NodeMap<double>& posx,
			     ListGraph::NodeMap<double>& posy)
{
  size_t t=0;
  double x,y;
  char tempname[1000],tempnamedot[1000],tempnameposdot[1000],cmd[1000];
  ofstream out;
  ifstream in;
  string linha,substring;

  // obtain a temporary file name
  tmpnam(tempname);
  strcpy(tempnamedot,tempname);   strcat(tempnamedot,".dot");
  strcpy(tempnameposdot,tempname);strcat(tempnameposdot,"_pos.dot");
  
  out.open(tempnamedot);
  if (!out.is_open()) return(false);

  out << "graph g {\n";
  out << "\tsize = \"11, 11\";\n";
  out << "\tnode [shape = \"circle\"];\n";
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    linha = "\t";    linha += IntToString(g.id(v));   linha += ";\n";
    out << linha;
  }
  for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
    linha = "\t";  linha += IntToString(g.id(g.u(a)));
    linha += "  -- ";  linha += IntToString(g.id(g.v(a))); linha += ";\n";
    out << linha;
  }
  out << "}\n";
  out.close();
  sprintf(cmd,"neato -Goverlap=false %s -o %s",tempnamedot,tempnameposdot); 
  system(cmd); // gera outro arquivo do neato, mas com posições

  in.open(tempnameposdot);
  if (!in.is_open()) return(false);
  while (!in.eof()) {
    getline(in,linha);
    t = linha.find("[bb=\"");
    if (t!=string::npos) break;
  }
  if (t==string::npos) {
    cout << "Temp Graphviz file is not appropriate for GenerateVertexPositions.\n";
    exit(0);
  }
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    getline(in,linha);
    t = linha.find("pos=\"");
    if (t!=string::npos) {
      stringstream s;
      substring = linha.substr(t+5);
      sscanf(substring.c_str(),"%lf,%lf",&x,&y);
      posx[v] = x;      posy[v] = y;
    } else {
      printf("GenerateVerexPositions: Error to obtain vertex coordinates.\n");
      return(false);
    }
  }
  return(true);
}


bool ReadGraph(ListGraph &g,
	       ListGraph::NodeMap<string>& nodename,
	       ListGraph::EdgeMap<double>& custo,
	       ListGraph::NodeMap<double>& posx,
	       ListGraph::NodeMap<double>& posy,
	       char *filename)
{
  int n,m;
  bool r;
  ifstream ifile;  
  ifile.open(filename);  if (!ifile) return(false);

  PulaBrancoComentario(ifile);
  ifile >> n;    ifile >> m;
  ifile.close();
  if ((n<=0)||(m<-1)) {
      cout << "Wrong number of vertices/edges in file " << filename << ".\n";
      exit(0);
  }
  if (m==-1) r = ReadEuclideanGraph(g,nodename,custo,posx,posy,filename);
  else       r = ReadListGraph(g,nodename,custo,posx,posy,filename);
  return(r);
}

  
