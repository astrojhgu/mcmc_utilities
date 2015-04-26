targets=test_fit test_arms1 test_arms2

all:$(targets)

test_fit:test_fit.o
	clang++ $< -o $@ -O2

test_fit.o:test_fit.cpp
	clang++ -c $< -O2

test_arms1:test_arms1.o
	clang++ $< -o $@ -O2

test_arms1.o:test_arms1.cpp
	clang++ -c $< -O2

test_arms2:test_arms2.o
	clang++ $< -o $@ -O2

test_arms2.o:test_arms2.cpp
	clang++ -c $< -O2

clean:
	rm -fv *.o
