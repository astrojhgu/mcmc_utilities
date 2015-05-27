targets=fit_beta2d.out

all:$(targets)

CXX=g++
CXXFLAGS=-I ../ -O3 -std=c++14 -g -fopenmp
LIBFLAGS=-lSDL2 -lfio -lcfitsio -fopenmp


%.out:%.o
	$(CXX) $< -o $@ $(LIBFLAGS)

%.o:	%.CXX
	$(CXX) -c $(CXXFLAGS) $< 

clean:
	rm -f $(targets) *.o
