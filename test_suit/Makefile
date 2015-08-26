targets=test_coal.out  test_dual_err.out  test_eff.out  test_gibbs.out  test_pf.out test_eff_graph.out test_eff_graph2.out test_normal.out

all:$(targets)

CXX=clang++
CXXFLAGS=-I ../ -O2 -std=c++14 -g -I ../../arms/
LIBFLAGS=-lSDL2 -lfio -lcfitsio -O3


%.out:%.o
	$(CXX) $< -o $@ $(LIBFLAGS)

%.o:	%.CXX
	$(CXX) -c $(CXXFLAGS) $< 

clean:
	rm -f $(targets) *.o
