targets=test_coal  test_eff  test_gibbs  test_pf test_eff_graph zero_trick test_poisson

all:$(targets)

CXX=clang++
CXXFLAGS=-I ../ -O2 -std=c++11 -g -I ../../arms/
LDFLAGS=-lSDL2 -lfio -lcfitsio -O3

% :	%.o
	$(CXX) $< -o $@ $(LDFLAGS)

%.o:	%.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)


clean:
	rm -f $(targets) *.o
