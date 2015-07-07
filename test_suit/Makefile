targets=test_arms1.out  test_coal.out  test_dual_err.out  test_dugongs.out  test_eff.out  test_gibbs.out  test_pf.out

all:$(targets)

CXX=clang++
CXXFLAGS=-I ../ -O3 -std=c++14 -g
LIBFLAGS=-lSDL2 -lfio -lcfitsio


%.out:%.o
	$(CXX) $< -o $@ $(LIBFLAGS)

%.o:	%.CXX
	$(CXX) -c $(CXXFLAGS) $< 

clean:
	rm -f $(targets) *.o
