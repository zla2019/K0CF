OBJS = MyTree.o Config.o Hist.o generateHists.o
EXE = generateHists

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOTLIBS     += -lEG

INCFLAGS = -I$(ROOTSYS)/include
LDFLAGS = -L$(ROOTSYS)/lib

CXX = g++ -m64
FLAGS = -Wall -g $(INCFLAGS) $(LDFLAGS) -fPIC

COMPILE = $(CXX) $(FLAGS) -c -std=c++14

all: $(EXE)

$(EXE): $(OBJS)
	$(CXX) -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS)

%.o: %.cxx
	$(COMPILE)  $< $(ROOTFLAGS) $(ROOTLIBS)

%.o: %.C
	$(COMPILE)  $< $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm -f *.o *.so generateHists
