targets=test_fit.out test_arms1.out test_arms2.out test_coal.out test_pf.out test_eff.out test_dual_err.out test_gibbs.out simulate_lm.out simulate_ddexp.out fit_lm.out simulate_dpar.out test_dugongs.out simulate_norm2d.out fit_norm2d.out fit_y_vs_x.out

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
