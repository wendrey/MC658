#include<lemon/math.h>
#include "thirdpartprograms.h"
#include "viewgraph.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

// To view a pdf file, you need to set the correct pdf reader.
// For example, in Mac OS (apple) system, we can open a pdf file useing
// open file.pdf
// In linux, we can use the "xpdf" program runing
// xpdf file.pdf
char pdfreader[]=PDF_READER; // pdfreader for MACOS
//char pdfreader[]="xpdf"; // pdfreader for LINUX

// View an undirected graph with DOT or NEATO program (both available in the 
// graphviz package. If you do not have the positions of the vertices of a graph,
// you can use the routine GenerateVertexPositions, available in the file 
// readgraph.cpp.
int ViewGraph(ListGraph &g,
	      int DOT_or_NEATO, // can be VIEW_DOT or VIEW_NEATO
	      ListGraph::NodeMap<string> &vname, // name of the vertices
	      ListGraph::EdgeMap<string> &ename,  // name of edges
	      ListGraph::NodeMap<double>& posx, // position x of the vertices
	      ListGraph::NodeMap<double>& posy, // position y of the vertices
	      ListGraph::NodeMap<int> &vcolor,   // color of the vertices
	      ListGraph::EdgeMap<int> &acolor)  // color of edges
{
  char tempname[1000],cmd[1000];
  FILE *fp;
  double minposx=10000000,minposy=10000000,maxposx=-10000000,maxposy=-10000000,
    delta,factor;

  tmpnam (tempname);
  fp = fopen(tempname,"w+");
  // fp = stdout;
  if (fp==NULL) {
    cout << "Erro ao abrir arquivo para visualizar o grafo.\n";
    return(0);
  }
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    if (posx[v] < minposx) minposx = posx[v];
    if (posx[v] > maxposx) maxposx = posx[v];
    if (posy[v] < minposy) minposy = posy[v];
    if (posy[v] > maxposy) maxposy = posy[v];
  }
  factor = 30;
  delta = fmax(maxposx - minposx,maxposy - minposy);

  fprintf(fp,"graph g {\n");
  // fprintf(fp,"\t\"start\" ;\n");
  fprintf(fp,"\tsize = \"10, 10\";\n");
  fprintf(fp,"\tnode [shape = \"circle\"];\n");
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    fprintf(fp,"\t%s [color=\"%s\", pos = \"%lf,%lf!\" ];\n",vname[v].c_str(),ColorName[vcolor[v]].c_str(),factor*(posx[v]-minposx)/delta,factor*(posy[v]-minposy)/delta);
  }
  for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
    if (acolor[a]!=WHITE)
      fprintf(fp,"\t%s  -- %s [label = \"%s\", color=\"%s\" ];\n",vname[g.u(a)].c_str(),vname[g.v(a)].c_str(),ename[a].c_str(),ColorName[acolor[a]].c_str());
  }
  fprintf(fp,"}\n");
  fclose(fp);

  if (DOT_or_NEATO==VIEW_NEATO) sprintf(cmd,"neato -Tpdf %s -o %s.pdf",tempname,tempname); 
  else sprintf(cmd,"dot -Tpdf %s -o %s.pdf",tempname,tempname);
  system(cmd);

  sprintf(cmd,"%s %s.pdf",pdfreader,tempname); system(cmd);
  return(1);
}
int ViewDigraph(ListDigraph &g,
	      int DOT_or_NEATO,
	      ListDigraph::NodeMap<string> &vname, // nome dos vertices
	      ListDigraph::ArcMap<string> &ename,  // nome das arestas (e.g. peso dela)
	      ListDigraph::NodeMap<int> &vcolor,   // cor dos vertices
	      ListDigraph::ArcMap<int> &acolor)    // cor das arestas
{
  char tempname[1000],cmd[1000];
  FILE *fp;

  tmpnam (tempname);
  fp = fopen(tempname,"w+");
  // fp = stdout;
  if (fp==NULL) {
    cout << "Erro ao abrir arquivo para visualizar o grafo.\n";
    return(0);
  }
  fprintf(fp,"digraph g {\n");
  // fprintf(fp,"\t\"start\" ;\n");
  // fprintf(fp,"\tgraph;\n");
  fprintf(fp,"\tedge [fontsize=9];\n");
  fprintf(fp,"\tnode [style=filled, fontsize=10];\n");
  fprintf(fp,"\tnodesep = .3;\n");
  fprintf(fp,"\tranksep = .5;\n");
  //fprintf(fp,"\tsize = \"expand\";\n");
  fprintf(fp,"\tnode [shape = circle];\n");
  for (ListDigraph::NodeIt v(g); v!=INVALID; ++v) {
    fprintf(fp,"\t%s [color=%s, width=.3, height=.3];\n",vname[v].c_str(),ColorName[vcolor[v]].c_str());
  }
  for (ListDigraph::ArcIt a(g); a!=INVALID; ++a) {
    if (acolor[a]!=WHITE)
      fprintf(fp,"\t%s  -> %s [label = \"%s\", color=%s ];\n",vname[g.source(a)].c_str(),vname[g.target(a)].c_str(),ename[a].c_str(),ColorName[acolor[a]].c_str());
  }
  fprintf(fp,"}\n");
  fclose(fp);

  if (DOT_or_NEATO==VIEW_NEATO) sprintf(cmd,"neato -Tpdf %s -o %s.pdf",tempname,tempname); 
  else sprintf(cmd,"dot -Tpdf %s -o %s.pdf",tempname,tempname);
  system(cmd);

  sprintf(cmd,"%s %s.pdf",pdfreader,tempname); system(cmd);
  return(1);
}


/*---------------------------------------------------------------------*/
#define MAXPOINTPOSITION 6000 /* cada coordenada estara' 
				  no intervalo 0..MAXPOINTPOSITION*/
void getepscolor(char *epscolorname,int cor)
{
  switch(cor) {
  case BLACK: strcpy(epscolorname,"col0"); return;
  case BLUE: strcpy(epscolorname,"col1"); return;
  case GREEN: strcpy(epscolorname,"col2"); return;
  case RED: strcpy(epscolorname,"col4"); return;
  case WHITE: strcpy(epscolorname,"col7"); return;
  }
  strcpy(epscolorname,"col0"); /* sem cor definida fica com preto */
}



int ViewEuclideanGraph(ListGraph &g,
		       ListGraph::NodeMap<string> &vname, // nome dos vertices
		       ListGraph::NodeMap<double> &posx, // coord. x dos vertices
		       ListGraph::NodeMap<double> &posy, // coord. y dos vertices
		       ListGraph::NodeMap<int> &vcolor,  // cor dos vertices
		       ListGraph::EdgeMap<int> &ecolor)  // cor das arestas
{
  char tempname[1000],cmd[1000];
  FILE *fp;
  double gap,maxx, maxy, minx, miny,
    telax,posxu,posxv,posyu,posyv;
  char epscolor[100];

  tmpnam (tempname);
  fp = fopen(tempname,"w+");
  if (fp==NULL) {
    cout << "Erro ao abrir arquivo para visualizar o grafo.\n";
    return(0);
  }

  maxx = -999999;
  maxy = -999999;
  minx =  999999;
  miny =  999999;
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    if (posx[v] > maxx) maxx =  posx[v];
    if (posx[v] < minx) minx =  posx[v];
    if (posy[v] > maxy) maxy =  posy[v];
    if (posy[v] < miny) miny =  posy[v];
  }
  telax = 500;
  fprintf(fp,"%%!PS-Adobe-2.0 EPSF-2.0\n");
  fprintf(fp,"%%%%Title: x.eps\n");
  fprintf(fp,"%%%%Creator: fig2dev Version 3.2 Patchlevel 3c\n");
  fprintf(fp,"%%%%CreationDate: Thu Sep 12 13:02:34 2002\n");
  fprintf(fp,"%%%%For: fkm@hobbes.dcc.unicamp.br ()\n");
  fprintf(fp,"%%%%BoundingBox: 0 0 %d %d\n",(int) telax,(int) telax);
  fprintf(fp,"%%%%Magnification: 1.0000\n");
  fprintf(fp,"%%%%EndComments\n");
  fprintf(fp,"/$F2psDict 200 dict def\n");
  fprintf(fp,"$F2psDict begin\n");
  fprintf(fp,"$F2psDict /mtrx matrix put\n");
  fprintf(fp,"/col-1 {0 setgray} bind def\n");
  fprintf(fp,"/col0 {0.000 0.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col1 {0.000 0.000 1.000 srgb} bind def\n");
  fprintf(fp,"/col2 {0.000 1.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col3 {0.000 1.000 1.000 srgb} bind def\n");
  fprintf(fp,"/col4 {1.000 0.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col5 {1.000 0.000 1.000 srgb} bind def\n");
  fprintf(fp,"/col6 {1.000 1.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col7 {1.000 1.000 1.000 srgb} bind def\n");
  fprintf(fp,"/col8 {0.000 0.000 0.560 srgb} bind def\n");
  fprintf(fp,"/col9 {0.000 0.000 0.690 srgb} bind def\n");
  fprintf(fp,"/col10 {0.000 0.000 0.820 srgb} bind def\n");
  fprintf(fp,"/col11 {0.530 0.810 1.000 srgb} bind def\n");
  fprintf(fp,"/col12 {0.000 0.560 0.000 srgb} bind def\n");
  fprintf(fp,"/col13 {0.000 0.690 0.000 srgb} bind def\n");
  fprintf(fp,"/col14 {0.000 0.820 0.000 srgb} bind def\n");
  fprintf(fp,"/col15 {0.000 0.560 0.560 srgb} bind def\n");
  fprintf(fp,"/col16 {0.000 0.690 0.690 srgb} bind def\n");
  fprintf(fp,"/col17 {0.000 0.820 0.820 srgb} bind def\n");
  fprintf(fp,"/col18 {0.560 0.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col19 {0.690 0.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col20 {0.820 0.000 0.000 srgb} bind def\n");
  fprintf(fp,"/col21 {0.560 0.000 0.560 srgb} bind def\n");
  fprintf(fp,"/col22 {0.690 0.000 0.690 srgb} bind def\n");
  fprintf(fp,"/col23 {0.820 0.000 0.820 srgb} bind def\n");
  fprintf(fp,"/col24 {0.500 0.190 0.000 srgb} bind def\n");
  fprintf(fp,"/col25 {0.630 0.250 0.000 srgb} bind def\n");
  fprintf(fp,"/col26 {0.750 0.380 0.000 srgb} bind def\n");
  fprintf(fp,"/col27 {1.000 0.500 0.500 srgb} bind def\n");
  fprintf(fp,"/col28 {1.000 0.630 0.630 srgb} bind def\n");
  fprintf(fp,"/col29 {1.000 0.750 0.750 srgb} bind def\n");
  fprintf(fp,"/col30 {1.000 0.880 0.880 srgb} bind def\n");
  fprintf(fp,"/col31 {1.000 0.840 0.000 srgb} bind def\n");
  fprintf(fp,"\n");
  fprintf(fp,"end\n");
  fprintf(fp,"save\n");
  fprintf(fp,"newpath 0 %d moveto 0 0 lineto %d 0 lineto %d %d "
	  "lineto closepath clip newpath\n",(int) telax,(int) telax,(int) telax,(int) telax);
  fprintf(fp,"%d %d translate\n",-10,(int) telax+10);
  fprintf(fp,"1 -1 scale\n");
  fprintf(fp,"\n");
  fprintf(fp,"/cp {closepath} bind def\n");
  fprintf(fp,"/ef {eofill} bind def\n");
  fprintf(fp,"/gr {grestore} bind def\n");
  fprintf(fp,"/gs {gsave} bind def\n");
  fprintf(fp,"/sa {save} bind def\n");
  fprintf(fp,"/rs {restore} bind def\n");
  fprintf(fp,"/l {lineto} bind def\n");
  fprintf(fp,"/m {moveto} bind def\n");
  fprintf(fp,"/rm {rmoveto} bind def\n");
  fprintf(fp,"/n {newpath} bind def\n");
  fprintf(fp,"/s {stroke} bind def\n");
  fprintf(fp,"/sh {show} bind def\n");
  fprintf(fp,"/slc {setlinecap} bind def\n");
  fprintf(fp,"/slj {setlinejoin} bind def\n");
  fprintf(fp,"/slw {setlinewidth} bind def\n");
  fprintf(fp,"/srgb {setrgbcolor} bind def\n");
  fprintf(fp,"/rot {rotate} bind def\n");
  fprintf(fp,"/sc {scale} bind def\n");
  fprintf(fp,"/sd {setdash} bind def\n");
  fprintf(fp,"/ff {findfont} bind def\n");
  fprintf(fp,"/sf {setfont} bind def\n");
  fprintf(fp,"/scf {scalefont} bind def\n");
  fprintf(fp,"/sw {stringwidth} bind def\n");
  fprintf(fp,"/tr {translate} bind def\n");
  fprintf(fp,"/tnt {dup dup currentrgbcolor\n");
  fprintf(fp,"  4 -2 roll dup 1 exch sub 3 -1 roll mul add\n");
  fprintf(fp,"  4 -2 roll dup 1 exch sub 3 -1 roll mul add\n");
  fprintf(fp,"  4 -2 roll dup 1 exch sub 3 -1 roll mul add srgb}\n");
  fprintf(fp,"  bind def\n");
  fprintf(fp,"/shd {dup dup currentrgbcolor 4 -2 roll mul 4 -2 roll mul\n");
  fprintf(fp,"  4 -2 roll mul srgb} bind def\n");
  fprintf(fp," /DrawEllipse {\n");
  fprintf(fp,"	/endangle exch def\n");
  fprintf(fp,"	/startangle exch def\n");
  fprintf(fp,"	/yrad exch def\n");
  fprintf(fp,"	/xrad exch def\n");
  fprintf(fp,"	/y exch def\n");
  fprintf(fp,"	/x exch def\n");
  fprintf(fp,"	/savematrix mtrx currentmatrix def\n");
  fprintf(fp,"	x y tr xrad yrad sc 0 0 1 startangle endangle arc\n");
  fprintf(fp,"	closepath\n");
  fprintf(fp,"	savematrix setmatrix\n");
  fprintf(fp,"	} def\n");
  fprintf(fp,"\n");
  fprintf(fp,"/$F2psBegin {$F2psDict begin "
	  "/$F2psEnteredState save def} def\n");
  fprintf(fp,"/$F2psEnd {$F2psEnteredState restore end} def\n");
  fprintf(fp,"\n");
  fprintf(fp,"$F2psBegin\n");
  fprintf(fp,"%%%%Page: 1 1\n");
  fprintf(fp,"10 setmiterlimit\n");
  fprintf(fp," %10.8lf %10.8lf sc\n",(double) MAXPOINTPOSITION/100000.0,(double) MAXPOINTPOSITION/100000.0);
  fprintf(fp,"%%\n");
  fprintf(fp,"%% Fig objects follow\n");
  fprintf(fp,"%%\n");
  fprintf(fp,"25.000 slw\n");
  gap = 300;
  for (ListGraph::EdgeIt a(g); a!=INVALID; ++a) {
    ListGraph::Node u,v;
    if (ecolor[a]==WHITE) continue;

    u = g.u(a);   v = g.v(a);
    posxu = (int) (MAXPOINTPOSITION*((double)(posx[u] -minx))/
		   ((double) (maxx-minx)))+gap;
    posyu = MAXPOINTPOSITION - (int) (MAXPOINTPOSITION*((double)(posy[u] -miny))/
				      ((double) (maxy-miny)))+gap;
    posxv = (int) (MAXPOINTPOSITION*((double)(posx[v] -minx))/
		   ((double) (maxx-minx)))+gap;
    posyv = MAXPOINTPOSITION - (int) (MAXPOINTPOSITION*((double)(posy[v] -miny))/
				      ((double) (maxy-miny)))+gap;
    getepscolor(epscolor,ecolor[a]);
    fprintf(fp,"n %d %d m\n %d %d l gs %s s gr \n",(int) posxu,(int) posyu,(int) posxv,(int) posyv,epscolor);
  }
  for (ListGraph::NodeIt v(g); v!=INVALID; ++v) {
    posxv = (int) (MAXPOINTPOSITION*((double)(posx[v] -minx))/
		  ((double) (maxx-minx)))+gap;
    posyv = MAXPOINTPOSITION - (int) (MAXPOINTPOSITION*((double)(posy[v] -miny))/
				     ((double) (maxy-miny)))+gap;
    getepscolor(epscolor,vcolor[v]);
    /* 45 eh o raio (ambos)*/
    fprintf(fp,"n %d %d 45 45 0 360 DrawEllipse gs %s 1.00 shd ef gr \n", 
	    (int) posxv,(int) posyv,epscolor);
  }
  fprintf(fp,"$F2psEnd\n");
  fprintf(fp,"rs\n");
  fclose(fp);

  sprintf(cmd,"mv %s %s.eps",tempname,tempname);           system(cmd);
  sprintf(cmd,"convert %s.eps %s.pdf",tempname,tempname);  system(cmd);
  sprintf(cmd,"%s %s.pdf",pdfreader,tempname);  system(cmd);
  return(true);
}
