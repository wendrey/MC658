
CC = g++

#CFLAGS = -Wall -g  
CFLAGS = -Wall -O3  

# Using the CLP (Coin-OR Linear Programming) is installed
LIBDIR = -L/Users/fkm/tiziu/academico/otimizacao/coin-or/clp/lib
LDFLAGS = -lemon -lClp -lCoinUtils

LIBSOURCES = branch_cut.cpp readgraph.cpp viewgraph.cpp tempo.cpp adjacencymatrix.cpp
OBJLIB = $(LIBSOURCES:.cpp=.o)

EX = ex_shortestpath.cpp ex_tsp1_slides.cpp ex_tsp_directed.cpp ex_tsp_undirected.cpp ex_readgraph1.cpp ex_min_cut1.cpp ex_lp1.cpp ex_matchingbipartitegraph.cpp ex_color1.cpp ex_dieta1.cpp ex_viewgraph_min_cut1.cpp ex_branch_cut1.cpp ex_branch_cut2.cpp ex_adjacencymatrix.cpp generate_random_euclidean_graph.cpp
OBJEX = $(EX:.cpp=.o)

EXE = $(EX:.cpp=.e)

all: mylib.a $(OBJEX) $(EXE)

mylib.a: $(OBJLIB)
	libtool -o $@ $(OBJLIB) 

%.o: %.cpp 
	$(CC) -c $^ $(CFLAGS) -o $@  

%.e: %.o mylib.a
	$(CC) mylib.a $^ $(CFLAGS) $(LDFLAGS) $(LIBDIR) -o $@

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJ) $(OBJLIB) $(EXE) $(OBJEX) *~ core mylib.a
