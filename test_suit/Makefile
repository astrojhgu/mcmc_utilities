targets=test_fit test_arms1 test_arms2 test_coal test_pf

all:$(targets)

cppflags=-I ../

test_fit:test_fit.o
	clang++ $< -o $@ -O2

test_fit.o:test_fit.cpp
	clang++ -c $(cppflags) $< -O2

test_arms1:test_arms1.o
	clang++ $< -o $@ -O2

test_arms1.o:test_arms1.cpp
	clang++ -c $(cppflags) $< -O2

test_arms2:test_arms2.o
	clang++ $< -o $@ -O2

test_arms2.o:test_arms2.cpp
	clang++ -c $(cppflags) $< -O2

test_coal:test_coal.o
	clang++ $< -o $@ -O2

test_coal.o:test_coal.cpp
	clang++ -c $(cppflags) $< -O2

test_pf:test_pf.o
	clang++ $< -o $@ -O2 -lSDL2

test_pf.o:test_pf.cpp
	clang++ -c $(cppflags) $< -O2


clean:
	rm -fv *.o
