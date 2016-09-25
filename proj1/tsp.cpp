/*******************************************************************************
 * MC658 - Projeto e Análise de Algoritmos III - 2s2016
 * Prof: Flavio Keidi Miyazawa
 * PED: Mauro Henrique Mulati
 * Usa ideias e código de Rafael Arakaki e Flávio Keidi Miyazawa 
 ******************************************************************************/

/*******************************************************************************
 * ATENÇÃO: NÃO ALTERE ESTE ARQUIVO
 ******************************************************************************/

#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "myutils.h"
#include "tsp.h"
#include "tsp_bt_bnb.h"

//------------------------------------------------------------------------------
TSP_Data::TSP_Data(ListGraph &graph,
                   NodeStringMap &nodename,
                   NodePosMap &posicaox,
                   NodePosMap &posicaoy,
                   EdgeValueMap &eweight):
                   g(graph),
                   vname(nodename),
                   ename(graph),
                   vcolor(graph),
                   ecolor(graph),
                   weight(eweight),
                   posx(posicaox),
                   posy(posicaoy),
                   AdjMat(graph,eweight,MY_INF),
                   BestCircuit(countEdges(graph))
{
	NNodes=countNodes(this->g);
	NEdges=countEdges(this->g);
	BestCircuitValue = DBL_MAX;
	// max_perturb2opt_it = 3000;  // default value
}
//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	// Variables to be obtained from parameters
	int    exec    = 0;      // 0: not set, 1: BT, 2: BNB
	int    maxTime = 0;      // 0: not set
   bool verbose   = false;
   string inputFile_name;   // Input graph file
   string outputFile_name;  // Output sol file
	
	// Reading program arguments
   for(int i = 1; i < argc; i++){
      const string arg(argv[i]);
      string next;
      if((i+1) < argc){
         next = string(argv[i+1]);
		}
      else{
         next = string("");
		}
		
      if( exec != 0 && (arg.find("-k") == 0 || arg.find("-a") == 0 ) ){
         cout << "Erro ao ler parametro \"" << arg << "\", somente pode haver um parametro de modo de execucao." << endl;
         showUsage();
         exit(1);
      }
      else if( arg.find("-k") == 0 ){
         exec = 1;
      }
      else if( arg.find("-a") == 0 ){
         exec = 2;
      }
      else if( arg.find("-v") == 0 ){
        verbose = true;
      }
      else if( arg.find("-t") == 0 && next.size() > 0){
         maxTime = atoi(next.c_str()); i++; continue;
      }
      else if( arg.find("-i") == 0 && next.size() > 0){
         inputFile_name = next; i++; continue;
      }
      else if( arg.find("-o") == 0 && next.size() > 0){
         outputFile_name = next; i++; continue;
      }
      else{
         cout << "Parametro invalido: \"" << arg << "\"" << " (ou faltando argumento)" << endl;
         showUsage();
         exit(1);
      }
   }

   // Required parameters
   if( exec == 0 ){
      cout << "Nenhum modo de execucao selecionado dentre: -k ou -a" << endl;
      showUsage(); 
		exit(1);
   }
   
   if( inputFile_name.size() < 1 ){
      cout << ((inputFile_name.size() < 1)? "nome do arq de grafo, ":"") 
			  << endl;
      showUsage(); 
		exit(1);
   }

   if( outputFile_name.size() < 1 ){
      cout << ((outputFile_name.size() < 1)? "nome do arq de saida.":"") 
			  << endl;
      showUsage(); 
		exit(1);
   }
   
   if( maxTime == 0 ){
		maxTime = 600;  // Default of 600s = 10m
      // cout << "Argumentos obrigatorios faltando: " 
		//      << ((maxTime == 0)?"-t <tempo_limite_segundos> ":"")
      //      << endl;
      // showUsage(); 
		// exit(1);
   }

	// int seed=1;     // mhmulati
	// srand48(seed);  // mhmulati
	
	// Variables that represent the weighted undirected graph of the asymmetric tsp
	ListGraph     g;
	EdgeValueMap  weight(g);
	NodeStringMap vname(g);
	NodePosMap    posx(g), 
	              posy(g);

	// Read the graph from the imput file
	if ( ! ReadListGraph(inputFile_name, g, vname, weight, posx, posy) ) {
	  cout << "Erro na leitura do arquivo de entrada " << inputFile_name << endl;
	  exit(1);
	}
	
	// Init the graph data structure
	TSP_Data tsp(g, vname, posx, posy, weight);
	
   double elapsed_time = DBL_MAX;
   clock_t before = clock();
   bool foundOptimalSolution = false;
	
		vector<Node> BestCircuit; // vector containing the best circuit found
	double BestCircuitValue;
	
	switch(exec){
		case 1:{
			foundOptimalSolution = bt(tsp, maxTime);
			break;
		}
		case 2:{
			foundOptimalSolution = bnb(tsp, maxTime);
			break;
		}
	}
   clock_t after = clock();
   elapsed_time = (double) (after-before) / CLOCKS_PER_SEC;

	// Imprimir a solucao em arquivo
   writeOutputFile(outputFile_name, inputFile_name, tsp, elapsed_time, maxTime, exec, foundOptimalSolution);
	
  // Imprime a solução na tela
   if(verbose){
      if( tsp.BestCircuitValue < DBL_MAX ){
			if(exec == 1){
				cout << "BACKTRACKING";
			}
			else{
				cout << "BRANCH AND BOUND";
			}
			
			cout << endl;
			cout << "Optimal?        : " <<  foundOptimalSolution << endl;
			cout << "Instance        : " << inputFile_name << endl;
         cout << "BestCircuitValue: " << tsp.BestCircuitValue << endl;
			cout << "ElapsedTime     : " << elapsed_time << endl;
			cout << "TimeLimit       : " << maxTime << endl;
			// ViewTspCircuit(tsp);
         cout << endl;
      }
      else{
         cout << "Nenhuma solucao viavel encontrada." << endl;
			return 1;
      }
   }

    return 0;
}
//------------------------------------------------------------------------------
void ViewTspCircuit(TSP_Data &tsp)
{
	ListGraph h;
	NodeStringMap h_vname(h);  // node names
	NodeNodeMap g2h(tsp.g);  // maps a node of g to a node of h
	NodePosMap h_posx(h);
	NodePosMap h_posy(h);
	NodeColorMap vcolor(h);   // color of the vertices
	EdgeColorMap acolor(h);  // color of edges
	EdgeStringMap aname(h);  // name of edges
	for (NodeIt v(tsp.g); v!=INVALID; ++v) {
		Node hv;
		hv = h.addNode();
		g2h[v] = hv;
		h_posx[hv] = tsp.posx[v];
		h_posy[hv] = tsp.posy[v];
		h_vname[hv] = tsp.vname[v];
		vcolor[hv] = BLUE;
	}
	for (int i=0;i<tsp.NNodes;i++) {
		Node u,v;
		Edge a;
		u = tsp.BestCircuit[i]; 
		v = tsp.BestCircuit[(i+1) % tsp.NNodes]; 
		a = h.addEdge(g2h[u] , g2h[v]);
		aname[a] = "";
		acolor[a] = BLUE;
	}
	ViewListGraph(h,h_vname,aname,h_posx,h_posy,vcolor,acolor,"TSP Circuit with cost "+DoubleToString(tsp.BestCircuitValue));
}
//------------------------------------------------------------------------------
void showUsage()
// Usage information
{
	cout << "Usage:"
	     << "./tsp <modo_operacao> (um dentre: -k backtracking, -a branch_and_bound) -t <tempo_max_em_segundos> {-v: mostra solução na tela}"
	     << "-i <nome_arquivo_entrada> -o <nome_arquivo_saida>"
	     << endl;
}
//------------------------------------------------------------------------------
void writeOutputFile(string outputfile, string graphname, TSP_Data &tsp, double elapsed_time, int max_time, int exec, bool opt)
{
   ofstream myfile;
   myfile.open(outputfile);
	if(exec == 1){
		myfile << "BACKTRACKING" << '\t'<< opt << '\t' << graphname << '\t';
		if( tsp.BestCircuitValue < DBL_MAX ){
			myfile << tsp.BestCircuitValue << '\t' << elapsed_time << '\t' << max_time << endl;
		}
		else{
			myfile << "NAO ENCONTROU SOLUCAO VIAVEL" << endl;
		}
	}
	else{
		myfile << "BRANCH AND BOUND" << '\t'<< opt  << '\t' << graphname << '\t';
		if( tsp.BestCircuitValue < DBL_MAX ){
			myfile << tsp.BestCircuitValue << '\t' << elapsed_time << '\t' << max_time << endl;
		}
		else{
			myfile << "NAO ENCONTROU SOLUCAO VIAVEL" << endl;
		}
	}
   
   myfile.close();
}
//------------------------------------------------------------------------------
