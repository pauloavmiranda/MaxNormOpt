
INCLUDE=./include
SRC=./src
OBJ=./obj


all: OptMaxNorm OptMaxNorm_graph OptMaxNorm_rag

#Compiladores
CC=gcc
CXX=g++

FLAGS= -O3 -Wall
#FLAGS= -Wall -g

LFLAGS = -lpthread -lz

#Bibliotecas
GLIPLIB  = -L./lib/glip/lib -lglip
GLIPFLAGS  = -I./lib/glip/include

#Rules
libglip:
	$(MAKE) -C ./lib/glip

OptMaxNorm: \
OptMaxNorm.cpp libglip
	$(CXX) $(FLAGS) $(GLIPFLAGS) -I$(INCLUDE) OptMaxNorm.cpp \
	$(GLIPLIB) -o OptMaxNorm -lm $(LFLAGS)

OptMaxNorm_graph: \
OptMaxNorm_graph.cpp libglip
	$(CXX) $(FLAGS) $(GLIPFLAGS) -I$(INCLUDE) OptMaxNorm_graph.cpp \
	$(GLIPLIB) -o OptMaxNorm_graph -lm $(LFLAGS)

OptMaxNorm_rag: \
OptMaxNorm_rag.cpp libglip
	$(CXX) $(FLAGS) $(GLIPFLAGS) -I$(INCLUDE) OptMaxNorm_rag.cpp \
	$(GLIPLIB) -o OptMaxNorm_rag -lm $(LFLAGS)

clean:
	$(RM) *~ *.o $(OBJ)/*.o OptMaxNorm OptMaxNorm_graph OptMaxNorm_rag -f
