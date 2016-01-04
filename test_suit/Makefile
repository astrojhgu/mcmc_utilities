targets=test_arms.out test_coal.out  test_eff.out  test_gibbs.out  test_pf.out test_eff_graph.out zero_trick.out test_poisson.out test_composed.out test_t.out test_cosmos.out

all:$(targets)

CXX=g++
CXXFLAGS=-I ../ -std=c++11 -g -I ../../east/include -I ../../coscalcpp/include/ -I ../rng -O3
LDFLAGS=-L ../../east/lib -least -L ../../coscalcpp/lib/ -lcoscalcpp -O3 -fopenmp

%.out :	%.o
	$(CXX) $< -o $@ $(LDFLAGS)

%.o:	%.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

test_pf.out:test_pf.o
	$(CXX) $< -o $@  -lSDL2 $(LDFLAGS)

test_pf.o:test_pf.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -f $(targets) *.o
