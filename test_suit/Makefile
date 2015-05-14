targets=test_fit test_arms1 test_arms2 test_coal test_pf test_eff test_dual_err test_gibbs

all:$(targets)

cpp=g++
cppflags=-I ../ -O2


test_fit:test_fit.o
	$(cpp) $< -o $@ 

test_fit.o:test_fit.cpp
	$(cpp) -c $(cppflags) $< 

test_arms1:test_arms1.o
	$(cpp) $< -o $@ 

test_arms1.o:test_arms1.cpp
	$(cpp) -c $(cppflags) $< 

test_arms2:test_arms2.o
	$(cpp) $< -o $@ 

test_arms2.o:test_arms2.cpp
	$(cpp) -c $(cppflags) $< 

test_coal:test_coal.o
	$(cpp) $< -o $@ 

test_coal.o:test_coal.cpp
	$(cpp) -c $(cppflags) $< 

test_eff:test_eff.o
	$(cpp) $< -o $@ 

test_eff.o:test_eff.cpp
	$(cpp) -c $(cppflags) $< 

test_dual_err:test_dual_err.o
	$(cpp) $< -o $@ 

test_dual_err.o:test_dual_err.cpp
	$(cpp) -c $(cppflags) $< -std=c++14

test_gibbs:test_gibbs.o
	$(cpp) $< -o $@ 

test_gibbs.o:test_gibbs.cpp
	$(cpp) -c $(cppflags) $< -std=c++14

test_pf:test_pf.o
	$(cpp) $< -o $@  -lSDL2

test_pf.o:test_pf.cpp
	$(cpp) -c $(cppflags) $< 


clean:
	rm -fv *.o
