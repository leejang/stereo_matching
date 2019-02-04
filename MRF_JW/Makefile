VERSION = MRF_JW

SRC =  mrf.cpp GraphCuts.cpp graph.cpp maxflow.cpp \
       LinkedBlockList.cpp LoopyBP.cpp

CC = g++

WARN = -W -Wall
OPT ?= -O3
CPPFLAGS = $(OPT) $(WARN) -DUSE_64_BIT_PTR_CAST
#CPPFLAGS = $(OPT) $(WARN)   ### use this line instead to compile on 32-bit systems

OBJ = $(SRC:.cpp=.o)

all: libMRF.a

libMRF.a: $(OBJ)
	rm -f libMRF.a
	ar ruc libMRF.a $(OBJ)
	ranlib libMRF.a

clean: 
	rm -f $(OBJ) core core.* *.stackdump *.bak

allclean: clean
	rm -f libMRF.a example example.exe

depend:
	@makedepend -Y -- $(CPPFLAGS) -- $(SRC) 2>> /dev/null

# DO NOT DELETE THIS LINE -- make depend depends on it.

mrf.o: mrf.h
GraphCuts.o: energy.h graph.h block.h mrf.h GraphCuts.h
GraphCuts.o: LinkedBlockList.h
graph.o: graph.h block.h mrf.h
maxflow.o: graph.h block.h mrf.h
LoopyBP.o: LoopyBP.h mrf.h
