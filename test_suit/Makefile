targets=test_coal.out  test_eff.out  test_gibbs.out  test_pf.out test_eff_graph.out zero_trick.out test_poisson.out test_composed.out

all:$(targets)

CXX=clang++
CXXFLAGS=-I ../ -O2 -std=c++11 -g -I ../../east/include
LDFLAGS=-lSDL2 -lfio -lcfitsio -O3 -L ../../east/lib -least

%.out :	%.o
	$(CXX) $< -o $@ $(LDFLAGS)

%.o:	%.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)


clean:
	rm -f $(targets) *.o
